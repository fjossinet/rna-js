var parsers = exports,
    EventEmitter = require('events').EventEmitter,
    fs = require('fs'),
    path = require('path'),
    carrier = require('carrier'),
    RNA = require('../core/molecules').RNA,
    DNA = require('../core/molecules').DNA,
    SecondaryStructure = require('../core/features').SecondaryStructure,
    BaseBaseInteraction = require('../core/features').BaseBaseInteraction,
    TertiaryStructure = require('../core/features').TertiaryStructure,
    utils = require('../core/utils'),
    os = require("os");
    
parsers.StringParser = function() {
    var self = this;

    this.parseGenbank = function(genbankData) {
        var lines = genbankData.split('\n'),
            data = {};

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];
            if (line.trim().match(/^ORGANISM/)) {
               data.organism = line.trim().split("ORGANISM")[1].trim();
            }
        }

        self.emit('end', data);
    };

    this.parseStockholm = function(stockholmData) {
        var lines = stockholmData.split('\n'),
            alignedSequences = {},
            aligned2D = "",
            consensusStructure = [],
            rnas = [],
            rfamId = undefined;

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i].trim(),
                tokens = line.split(new RegExp(" +"));
            if (line.length != 0 && !line.match(/^#/) && tokens.length == 2) {
                if (tokens[0] in alignedSequences)
                    alignedSequences[tokens[0]] = alignedSequences[tokens[0]]+tokens[1];
                else
                    alignedSequences[tokens[0]] = tokens[1];
            } else if (line.length != 0 && line.match(/^#=GC SS_cons/))
                aligned2D += tokens[2].replace(/</g, '(').replace(/>/g, ')');
            else if (line.length != 0 && line.match(/^#=GF AC/))
                rfamId = tokens[2].trim();
        }

        var i = 0,
            lastPairedPos = [],
            basePairs = [];

        aligned2D.split('').forEach(function(pos) {
            i++;
            if (pos == '(')
                lastPairedPos.push(i);
            else if (pos == ')')
                basePairs.push(new BaseBaseInteraction('c', '(',')', lastPairedPos.pop(), i));
        });

        for (var sequenceName in alignedSequences) {
            var rna = new RNA(sequenceName, alignedSequences[sequenceName]);
            rna.source ="db:rfam:"+rfamId;
            rnas.push(rna);
        }

        self.emit('end', rnas, basePairs);
    };

    this.parseSam = function(samData) {
        var lines = samData.split('\n'),
            reads = [];

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i],
                tokens = line.split('\t');
            if (tokens.length >= 10) {
                var ncRNA = new RNA("RNA-Seq read");
                ncRNA.source = "file:toto:toto";
                ncRNA.genomicStrand = '+';
                ncRNA.sequence = tokens[9].trim();
                ncRNA.genomicPositions = [parseInt(tokens[3].trim()), parseInt(tokens[3].trim())+ncRNA.sequence.length-1];
                ncRNA.score = 0;
                ncRNA.organism = "Kuca";
                ncRNA.family = "RNA-Seq read";
                ncRNA.tool = "RNA-Seq read";
                ncRNA.genome = "506a6dfb45e3e6a41d000002@genomes";                    
                reads.push(ncRNA);
            }
        }

        self.emit('end', reads);
    };

    this.parseFasta = function(fastaData, molecularType) {
        var lines = fastaData.split('\n'),
            molecules = [],
            molecule = null,
            seq = "",
            type = molecularType || 'RNA';

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];
            if (line.charAt(0) == '>') {
                if (molecule) {
                    seq.toUpperCase().split('').forEach(function(residue) {
                        molecule.addResidue(residue);
                    });
                    molecules.push(molecule);   
                }
                if (type == 'RNA') {
                    molecule = new RNA(line.substr(1).trim());
                } else if (type == 'DNA') {
                    molecule = new DNA(line.substr(1).trim());
                }
                seq = "";
            } else {
                seq += line.trim();
            }
        }
        //the last molecule
        if (molecule) {
            seq.toUpperCase().split('').forEach(function(residue) {
                molecule.addResidue(residue);
            });
            molecules.push(molecule);   
        }
        self.emit('end',molecules);
    };
    
    this.parseVienna = function(viennaData) {
        var lines = viennaData.split('\n'),
            rna = null,
            seq = "",
            bn = "",
            regexp = /[.()]+/,
            basePairs = [],
            secondaryStructures = [];
            
        for (var i = 0; i < lines.length; i++) {
            var line = lines[i];
            if (line.charAt(0) == '>') {
                if (rna) {
                    seq.toUpperCase().split('').forEach(function(residue) {rna.addResidue(residue)});
                    secondaryStructures.push(new SecondaryStructure(rna, utils.bnToBps(bn)));   
                }
                rna = new RNA(line.substr(1).trim());
                sequence = "";
                bn = "";
                basePairs = [];
            } else if (line.match(regexp)) {
                bn = bn.concat(line.trim());
            } else {
                seq = seq.concat(line.trim());
            }
        }
        //the last rna
        if (rna) {
            seq.toUpperCase().split('').forEach(function(residue) {rna.addResidue(residue)});
            secondaryStructures.push(new SecondaryStructure(rna, utils.bnToBps(bn)));
        }
        self.emit('end',secondaryStructures);
    };

    this.parsePDB = function(pdbData) {
        var lines = pdbData.split('\n'),
            molecules = [],
            tertiaryStructures = [],
            currentChain = '', 
            currentResidue = '', 
            currentRNA = null, 
            currentTs = null,
            currentResiduePos = null, 
            absolutePosition = 0,
            title = ""

        for (var i = 0; i < lines.length; i++) {
            var line = lines[i],
                header = line.substring(0,6).trim();
                atomName = line.substring(12,16).trim(),
                residueName = line.substring(17,20).trim(),
                chainName = line.substring(21,22).trim(),
                residuePos = line.substring(22,27).trim();

            if ((header == "ATOM" || header == "HETATM") && ["FMN","PRF","HOH","MG","OHX","MN","ZN", "SO4", "CA", "UNK"].indexOf(residueName) == -1 && ["MG","K", "NA", "SR", "CL", "CD", "ACA"].indexOf(atomName) == -1 && chainName.length != 0) {
                if (currentChain != chainName) {
                    //new chain
                    currentResidue = residueName;
                    currentResiduePos = residuePos;
                    currentChain = chainName;
                    absolutePosition = 1;
                    currentRNA = new RNA(currentChain);
                    currentRNA.on('error', function(error) {
                        self.emit('error', error);
                    });
                    currentRNA.addResidue(currentResidue);
                    currentTs = new TertiaryStructure(currentRNA);
                    currentTs.title = title.replace( /  +/g, ' ' );
                    currentTs.numbering_system[currentResiduePos] = absolutePosition;
                } else if (currentResiduePos != residuePos) {
                    //new residue
                    currentResidue = residueName;
                    currentResiduePos = residuePos;
                    currentRNA.addResidue(currentResidue);
                    absolutePosition++;
                    currentTs.numbering_system[currentResiduePos] = absolutePosition;
                }
                var x = parseFloat(line.substring(30,38).trim()),
                    y = parseFloat(line.substring(38,46).trim()),
                    z = parseFloat(line.substring(46,54).trim());
                currentTs.addAtom(atomName,absolutePosition,[x,y,z]);
                if ((atomName == "O4'" || atomName == "O4*") && molecules.indexOf(currentRNA) == -1) {
                    molecules.push(currentRNA);
                    tertiaryStructures.push(currentTs);
                }
            } else if (header == "TITLE") {
                title = title.concat(line.substring(10));
            } else if (header == "JRNL") {

            }
        }
        self.emit('end',tertiaryStructures);
    };
    
}

parsers.StringParser.prototype.__proto__ = EventEmitter.prototype;

parsers.FileParser = function() {
    var self = this;

    this.parseCvs = function(cvsFile, sep) {
        var sep = sep || ';',
            inStream = fs.createReadStream(path.normalize(cvsFile), {flags:'r'}),
            lines = [];
        carrier.carry(inStream).
            on('line', function(line) {
                lines.push(line.split(sep));
            }).on('end',function() {
                self.emit('end', lines);    
            });   
    };

    this.parseSam = function(samFile) {
        var inStream = fs.createReadStream(path.normalize(samFile), {flags:'r'}),
        genomic_locations = [];

        carrier.carry(inStream).
            on('line', function(line) {
                if (!line.match(/^@/)) {
                    var tokens = line.split(/\t/);
                    if (tokens[2] != "*") {
                        var annotation = {};
                        annotation.id = tokens[0].trim();
                        annotation.targetSequence = tokens[2].trim();
                        annotation.start = parseInt(tokens[3].trim());
                        annotation.strand = tokens[1].trim() == "16" ? '-' : '+';
                        annotation.end = annotation.start+tokens[9].trim().length-1;
                        annotation.sequence = tokens[9].trim();
                        self.emit('got new annotation', annotation);    
                    }
                }
            }).on('end',function() {
                self.emit('end');    
            });
    };

    this.parseGff3 = function(gff3File) { //specifications here http://www.sequenceontology.org/gff3.shtml
        var inStream = fs.createReadStream(path.normalize(gff3File), {flags:'r'}),
            annotations = [];

        carrier.carry(inStream).
            on('line', function(line) {
                var tokens = line.split(/\t/);
                if (tokens.length > 1) {
                    var annotation = {};
                    annotation.targetSequence = tokens[0].trim();
                    annotation.type = tokens[2].trim();
                    annotation.start = parseInt(tokens[3].trim());
                    annotation.end = parseInt(tokens[4].trim());
                    if (!isNaN(tokens[5].trim()))
                        annotation.score = parseFloat(tokens[5].trim());
                    annotation.strand = tokens[6].trim();
                    annotation.phase = tokens[7].trim();
                    annotation.attributes = {};

                    tokens[8].trim().split(';').forEach(function(attribute) {
                        var pair = attribute.split('=');
                        annotation.attributes[pair[0].toLowerCase()] = pair[1];
                    });

                    annotations.push(annotation);
                }         
            }).
            on('end',function() {
                self.emit('end',annotations);
            });
    };
    
    this.parseFasta = function(fastaFile, molecularType) {
        var inStream = fs.createReadStream(path.normalize(fastaFile), {flags:'r'}),
            fastaData = "";

        carrier.carry(inStream).
            on('line', function(line) {
                fastaData += line+'\n';
            }).
            on('end',function() {
                var stringParser = new parsers.StringParser();
                stringParser.on('end', function(molecules) {
                    for (var i = 0 ; i < molecules.length ; i++) {
                        molecules[i].source = "file:"+os.hostname()+":"+path.basename(fastaFile);    
                    }
                    self.emit('end',molecules);
                });
                stringParser.parseFasta(fastaData, molecularType);
            });
    };

    this.parseBpseq = function(bpseqFile) {
        var molecules = [],
            secondaryStructures = [],
            basePairs = [],
            rna = new RNA(),
            inStream = fs.createReadStream(path.normalize(bpseqFile), {flags:'r'});

        carrier.carry(inStream).
        on('line', function(line) {
            var isDigit = /^\d+$/,
                tokens = line.split(/\s+/);
            if (tokens.length == 3 && isDigit.test(tokens[0])) {
                rna.addResidue(tokens[1]);
                if (parseInt(tokens[2]) != 0 && parseInt(tokens[0]) < parseInt(tokens[2])) {
                    basePairs.push([parseInt(tokens[0]),parseInt(tokens[2])]);
                }
            }   
        }).
        on('end',function() {
            rna.source = "file:"+os.hostname()+":"+bpseqFile;
            molecules.push(rna);
            var ss = new SecondaryStructure(rna, basePairs);
            ss.source = "file:"+os.hostname()+":"+path.basename(bpseqFile);
            secondaryStructures.push(ss);
            self.emit('end',molecules,secondaryStructures,[]);
        });
    };

    this.parseCt = function(ctFile) {
        var molecules = [],
            secondaryStructures = [],
            basePairs = [],
            rna = new RNA(),
            inStream = fs.createReadStream(path.normalize(ctFile), {flags:'r'});

        carrier.carry(inStream).
        on('line', function(line) {
            var isDigit = /^\d+$/,
                tokens = line.split(/\s+/);
            if (tokens.length == 6 && isDigit.test(tokens[0])) {
                rna.addResidue(tokens[1]);
                if (parseInt(tokens[4]) != 0 && parseInt(tokens[0]) < parseInt(tokens[4])) {
                    basePairs.push([parseInt(tokens[0]),parseInt(tokens[4])]);
                }
            }   
        }).
        on('end',function() {
            rna.source = "file:"+os.hostname()+":"+path.basename(ctFile);
            molecules.push(rna);
            var ss = new SecondaryStructure(rna, basePairs);
            ss.source = "file:"+os.hostname()+":"+ctFile;
            secondaryStructures.push(ss);
            self.emit('end',molecules,secondaryStructures,[]);
        });
    };


    this.parseVienna = function(viennaFile) {
        var inStream = fs.createReadStream(path.normalize(viennaFile), {flags:'r'}),
            viennaData = "";

        carrier.carry(inStream).
            on('line', function(line) {
                viennaData += line+'\n';   
            }).
            on('end',function() {
                var stringParser = new parsers.StringParser();
                stringParser.on('end', function(secondaryStructures) {
                    for (var i = 0 ; i < secondaryStructures.length ; i++) {
                        secondaryStructures[i].rna.molecules[i].source = "file:"+os.hostname()+":"+path.basename(viennaFile);
                        secondaryStructures[i].source = "file:"+os.hostname()+":"+path.basename(viennaFile);    
                    }
                    self.emit('end',secondaryStructures);
                });
                stringParser.parseVienna(viennaData);
            });
    };

    this.parsePdb = function(pdbFile) {
        var inStream = fs.createReadStream(path.normalize(pdbFile), {flags:'r'}),
            pdbData = "";

        carrier.carry(inStream).
            on('line', function(line) {
                pdbData += line+'\n';   
            }).
            on('end',function() {
                var stringParser = new parsers.StringParser();
                stringParser.on('end', function(tertiaryStructures) {
                    for (var i = 0 ; i < tertiaryStructures.length ; i++) {
                        tertiaryStructures[i].rna.source = "file:"+os.hostname()+":"+path.basename(pdbFile);
                        tertiaryStructures[i].source = "file:"+os.hostname()+":"+path.basename(pdbFile);    
                    }
                    self.emit('end',tertiaryStructures);
                });
                stringParser.on('error', function(error) {
                    self.emit('error', error);
                });
                stringParser.parsePDB(pdbData);
            });
    };
    
};

parsers.FileParser.prototype.__proto__ = EventEmitter.prototype;

parsers.FileWriter = function() {
    var self = this;
    
    this.writeFasta = function(molecules, file, single_line) {
        var outStream = fs.createWriteStream(path.normalize(file), {flags:'a'}),
        content = "";
        for (var i = 0 ; i < molecules.length ; i++ )
            content += molecules[i].toFasta(single_line)+'\n';
        outStream.write(content);
        self.emit('end');
    };
    
    this.writeVienna = function(secondaryStructure, file, single_line) {
        var outStream = fs.createWriteStream(path.normalize(file), {flags:'w'});
        outStream.write(secondaryStructure.toVienna(single_line));
        self.emit('end');
    };

    this.writeBpseq = function(secondaryStructure, file) {
        var outStream = fs.createWriteStream(path.normalize(file), {flags:'w'});
        outStream.write(secondaryStructure.toBpseq());
        self.emit('end');
    };

    this.writePdb = function(tertiaryStructure, file, exportNumberingSystem) {
        var outStream = fs.createWriteStream(path.normalize(file), {flags:'w'});
        tertiaryStructure.on('end pdb export', function(data) {
            fs.writeFile(file,data, function(err) {
                if(err) {
                    console.log(err);
                    self.emit('error');
                } else
                    self.emit('end');
            }); 
        });
        tertiaryStructure.toPdb(exportNumberingSystem || false);
        
    };

    this.writeRnaml = function(feature, file) {
        var builder = require('xmlbuilder'),
            doc = builder.create(),
            outStream = fs.createWriteStream(path.normalize(file), {flags:'w'});
        if (feature instanceof SecondaryStructure) {
            var str_annotation = doc.begin('rnaml')
                .ele('molecule')
                    .att('id','1')
                    .ele('sequence')
                        .ele('seq-data')
                            .txt(feature.rna.sequence)
                        .up()
                    .up()
                .up()
                .ele('structure')
                    .ele('model')
                        .att('id','1')
                        .ele('str-annotation');
            feature.helices.forEach(function(helix) {
              str_annotation
                .ele('helix')
                    .att('id',helix.name)
                    .ele('base-id-5p')
                        .ele('base-id')
                            .ele('position')
                                .txt(''+helix.location.ends[0][0])
                            .up()
                        .up()
                    .up()
                    .ele('base-id-3p')
                        .ele('base-id')
                            .ele('position')
                                .txt(''+helix.location.ends[1][0])
                            .up()
                        .up()
                    .up()
                    .ele('length')
                        .txt(''+(helix.location.ends[0][1]-helix.location.ends[0][0]+1)); 
            }); 
        }
        //outStream.write(doc.toString({ pretty: true }));
        console.log(doc.toString({ pretty: true }));
        self.emit('end');
    };
};

parsers.FileWriter.prototype.__proto__ = EventEmitter.prototype;
