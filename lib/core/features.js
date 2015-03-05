var features = exports,
    EventEmitter = require('events').EventEmitter,
    Location = require('./utils').Location,
    sprintf = require('sprintf'),
    ObjectID = require('mongodb').ObjectID;

features.BaseBaseInteraction = function(orientation, edge1, edge2, pos1, pos2) {
    this.orientation = orientation;
    this.edge1 = edge1;
    this.edge2 = edge2;
    this.location = new Location([pos1,pos1,pos2,pos2]);
};

features.SecondaryStructure = function(rna, basePairs) {
    this.name = "2D";
    this.rna = rna;
    this.helices = [];
    this.singleStrands = [];
    this.tertiaryInteractions = [];
    this.junctions = []
    this.source = "?:?:?";
    this._id = new ObjectID();
    
    var self = this,
        basePairs = basePairs || [];
    
    this.addHelix = function(name, start, end, length) {
        var helix = {
            name: name,
            location: new Location([start,start+length-1,end-length+1,end])
        };
        self.helices.push(helix);
        self.emit("new helix", helix);
        return helix;       
    };

    this.addSingleStrand = function(name, start, length) {
        var singleStrand = {
            name: name,
            location: new Location([start,start+length-1])
        };
        self.singleStrands.push(singleStrand);
        self.emit("new single strand", singleStrand);
        return singleStrand;  
    };

    this.addBasePair = function(orientation, edge1, edge2, pos1, pos2) {
        var isSecondaryInteraction = false;
        self.helices.forEach(function(helix) {
            var start = helix.location.getStart(),
                end = helix.location.getEnd(),
                length = helix.location.getLength()/2;

            if (pos1 >= start && pos1 <= start+length-1 ) {
                var diff = pos1 - start;
                if (end - pos2 == diff) {
                    //if not canonical (not AU, GC or GU, neither cWWW, we add it to the helix as a non-canonical secondary interaction
                    if (!(self.rna.sequence[pos1-1] == 'A' && self.rna.sequence[pos2-1] == 'U' ||
                          self.rna.sequence[pos1-1] == 'U' && self.rna.sequence[pos2-1] == 'A' ||
                          self.rna.sequence[pos1-1] == 'G' && self.rna.sequence[pos2-1] == 'C' ||
                          self.rna.sequence[pos1-1] == 'C' && self.rna.sequence[pos2-1] == 'G' ||
                          self.rna.sequence[pos1-1] == 'G' && self.rna.sequence[pos2-1] == 'U' ||
                          self.rna.sequence[pos1-1] == 'U' && self.rna.sequence[pos2-1] == 'G') ||
                          orientation != 'C' || edge1 != '(' || edge2 != ')') { //we have a non-canonical secondary-interaction
                        if (!helix.hasOwnProperty('interactions'))
                            helix['interactions'] = [];
                        helix['interactions'].push(new features.BaseBaseInteraction(orientation, edge1, edge2, pos1, pos2));
                    }
                    isSecondaryInteraction = true;
                }
            }  
        });
        if ( !isSecondaryInteraction) {
            //if we reach this point, its a tertiary interaction
            self.addTertiaryInteraction(orientation, edge1, edge2, pos1, pos2);
        }
    };

    this.addTertiaryInteraction = function(orientation, edge1, edge2, pos1, pos2) {
        var interaction = new features.BaseBaseInteraction(orientation, edge1, edge2, pos1, pos2);
        self.tertiaryInteractions.push(interaction);
        self.emit("new tertiary interaction", interaction);
    };
    
    this.toBn = function() {
        var bn = "";
        SEQUENCE: for (var pos = 1 ; pos <= self.rna.sequence.length ; pos++) {
            for (var i = 0 ; i < self.helices.length ; i++) {
                var helix = self.helices[i];
                var chunk = helix.location.ends[0];
                
                if (pos >= chunk[0] && pos <= chunk[1]) {
                    bn = bn.concat('(');
                    continue SEQUENCE;
                }
                
                chunk = helix.location.ends[1];
                
                if (pos >= chunk[0] && pos <= chunk[1]) {
                    bn = bn.concat(')');
                    continue SEQUENCE;
                }
            }
            for (var i = 0 ; i < self.tertiaryInteractions.length ; i++) {
                var interaction = self.tertiaryInteractions[i];
                
                if (pos == interaction.location.ends[0][0]) {
                    bn = bn.concat('(');
                    continue SEQUENCE;
                } else if (pos == interaction.location.ends[1][0]) {
                    bn = bn.concat(')');
                    continue SEQUENCE;
                }
            }
            bn = bn.concat('.');
        }
        return bn;
    };
    
    this.toVienna = function (single_line) {
        var toVienna = ">"+self.name+"\n",
            bn = self.toBn();
        if (single_line) {
            toVienna = toVienna.concat(self.rna.sequence+"\n");
            toVienna = toVienna.concat(bn+"\n");
        } else {
            var chunks = utils.chunks(self.rna.sequence.split(''),79);
            for (var i = 0; i < chunks.length ; i++) {
                toVienna = toVienna.concat(chunks[i].join("")+"\n");
            }
            chunks = utils.chunks(bn.split(''),79);
            for (var i = 0; i < chunks.length ; i++) {
                toVienna = toVienna.concat(chunks[i].join("")+"\n");
            }
        }
        return toVienna;
    };

    this.toBpseq = function() {
        var i = 1,
            toBpseq = "";
        self.rna.sequence.split('').forEach(function(residue) {
            toBpseq += i+"\t"+residue+"\t";
            var pairedPos = self.getPairedResidue(i);
            if (pairedPos)
                toBpseq += pairedPos+"\n";
            else  {
                var found = false; 
                for (var j = 0 ; i < self.tertiaryInteractions.length ; j++) {
                    var tertiaryInteraction = self.tertiaryInteractions[j]; 
                    if (i == tertiaryInteraction.location[0]) {
                        toBpseq += tertiaryInteraction.location[2]+"\n";
                        found = true;
                        break; 
                    }   
                }
                if (! found)
                    toBpseq +="0\n";   
            }
            i++;
        });
        return toBpseq;            
    };

    this.getPairedResidue = function(pos) {
        for (var i = 0 ; i < self.helices.length ; i++) {
            var helix = self.helices[i];
            if (pos >= helix.location.getStart() && pos <= (helix.location.getStart() + helix.location.getLength()/2-1)) {
                return helix.location.getEnd() - (pos - helix.location.getStart());
            } else if (pos <= helix.location.getEnd() && pos >= (helix.location.getEnd() - helix.location.getLength()/2+1)) {
                return helix.location.getStart() + (helix.location.getEnd() - pos);
            }
        }
        return undefined;
    };

    this.findJunctions = function() {
        self.junctions = [];
        for (var i = 0 ; i < self.singleStrands.length ; i++) {
            var singleStrand = self.singleStrands[i];
            if (singleStrand.location.getStart() == 1 || singleStrand.location.getEnd() == self.rna.sequence.length || self.junctions.filter(function(junction) { return junction.singleStrands.indexOf(singleStrand) != -1 }).length != 0 ) {// a single-strand cannot be in two different junctions
                continue;
            }
            var strands = [singleStrand],
                descr = rna.sequence.substring(singleStrand.location.getStart()-1,singleStrand.location.getEnd())+" ",
                currentPos = self.getPairedResidue(singleStrand.location.getEnd()+1)+1,
                crown = [[singleStrand.location.getStart()-1, singleStrand.location.getEnd()+1]];
            while (currentPos <= rna.sequence.length) {
                var nextSingleStrand = self.singleStrands.filter( function(singleStrand) {
                    return singleStrand.location.getStart() == currentPos ;
                }).pop();
                if (nextSingleStrand == singleStrand) {
                    break;
                } else if (nextSingleStrand) {
                    strands.push(nextSingleStrand);
                    crown.push([nextSingleStrand.location.getStart()-1, nextSingleStrand.location.getEnd()+1]);
                    descr = descr.concat(rna.sequence.substring(nextSingleStrand.location.getStart()-1,nextSingleStrand.location.getEnd())+" ");
                    currentPos = self.getPairedResidue(nextSingleStrand.location.getEnd()+1)+1;
                    continue;
                }
                var nextHelix = self.helices.filter( function(helix) {
                    return currentPos == helix.location.getStart() ||  helix.location.getEnd() -  helix.location.getLength()/2+1; 
                }).pop();
                if (nextHelix) {
                    descr = descr.concat("- ");
                    crown.push([currentPos-1, currentPos]);   
                    currentPos = self.getPairedResidue(currentPos)+1;
                }
            } 

            if (nextSingleStrand == singleStrand) { //we have a junction
                self.junctions.push({
                    singleStrands:strands, 
                    description:descr.trim(),
                    'crown': crown
                });
            }
        }
    };

    if (basePairs.length > 0) {

        //we sort the base-pairs with their first position in the molecule
        basePairs = basePairs.sort(function(basePair1,basePair2) {
            return basePair1[0]-basePair2[0];
        });

        var helixLength = 0, 
            start = -1, 
            end = -1, 
            helixNb = 1;
            
        //construction of helices
        var ssLocation = new Location([1,rna.sequence.length]);
        for (var i = 0 ; i < basePairs.length ; i++ ) {
            if (basePairs[i][0] == start+helixLength && basePairs[i][1] == end-helixLength) {
                helixLength++;
            } else {
                if (helixLength > 1) {
                    var helix = this.addHelix("H"+helixNb++, start, end, helixLength);
                    ssLocation.removeLocation(helix.location);
                } else if (start != -1 && end != -1) {
                    this.addTertiaryInteraction('c','(',')', start, end);
                }
                start = basePairs[i][0];
                end = basePairs[i][1];
                helixLength = 1;
            }   
        }
        //last helix
        if (helixLength > 1) {
            var helix = this.addHelix("H"+helixNb++, start, end, helixLength);
            ssLocation.removeLocation(helix.location);
        } else if (helixLength == 1) {
            this.addTertiaryInteraction('c','(',')', start, end);
        }

        for (var i = 0; i < ssLocation.ends.length ; i++) {
            var chunk = ssLocation.ends[i];
            self.addSingleStrand("SS"+(i+1),chunk[0], chunk[1]-chunk[0]+1); 
        }
    }
};

features.SecondaryStructure.prototype.__proto__ = EventEmitter.prototype;

features.TertiaryStructure = function(rna) {
    var self = this;
    this.name = "3D";
    this.rna = rna;
    this.residues = {}; //the keys are thee absPos for residues
    this.source = "?:?:?";
    this.title = null;
    this.resolution = null;
    this.publication = null;
    this._id = null,
    this.numbering_system = {}; // a dict of (label => abs_pos)
    
    this.addAtom = function(name, position, coords) {
        name = name.replace("*","'");
        if (name == 'OP1')
            name = 'O1P';
        else if (name == 'OP2')
            name = 'O2P'; 
        else if (name == 'OP3')
            name = 'O3P';    
        if (self.residues.hasOwnProperty(position)) {
            self.residues[position].atoms.push({
                'name': name,
                coords: coords
            });
        } else {
              self.residues[position] = {
                    atoms: [{
                        'name': name,
                        coords: coords
                    }]
              }  
        }
    }

    this.toPdb = function(exportNumberingSystem) {
        var toPdb = "",
            i = 1,
            keys =[];

        for (var key in self.residues)
            keys.push(key);

        keys.sort(function(a,b){return a - b}); //the positions are sorted numerically

        keys.forEach(function(key) {
            self.residues[key].atoms.forEach(function(atom) {
                if (exportNumberingSystem)
                    toPdb = toPdb.concat(sprintf.sprintf('%-6s%5u  %-4s%3s %s%4s    %8.3f%8.3f%8.3f\n',"ATOM", i++, atom.name, self.rna.sequence[parseInt(key)-1], self.rna.name[0], self.getLabelByAbsPos(parseInt(key)), atom.coords[0], atom.coords[1], atom.coords[2]));
                else
                    toPdb = toPdb.concat(sprintf.sprintf('%-6s%5u  %-4s%3s %s%4u    %8.3f%8.3f%8.3f\n',"ATOM", i++, atom.name, self.rna.sequence[parseInt(key)-1], self.rna.name[0], parseInt(key), atom.coords[0], atom.coords[1], atom.coords[2]));
            });    
        })
        

        self.emit('end pdb export', toPdb);

    }

   this.getLabelByAbsPos = function( absPos ) {
    for( var prop in self.numbering_system ) {
        if( self.numbering_system.hasOwnProperty( prop ) ) {
             if( self.numbering_system[ prop ] === absPos )
                 return prop;
        }
    }
}

};

features.TertiaryStructure.prototype.__proto__ = EventEmitter.prototype;