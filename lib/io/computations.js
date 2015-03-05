var computations = exports,
    EventEmitter = require('events').EventEmitter,
    FileWriter = require('./parsers').FileWriter,
    StringParser = require('./parsers').StringParser,
    SecondaryStructure = require('../core/features').SecondaryStructure,
    TertiaryStructure = require('../core/features').TertiaryStructure,
    BaseBaseInteraction = require('../core/features').BaseBaseInteraction,
    RNA = require('../core/molecules').RNA,
    DNA = require('../core/molecules').DNA,
    db = require('./db'),
    exec = require('child_process'),
    fs = require('fs'),
    path = require('path'),
    carrier = require('carrier'),
    tmp = require('tmp'),
    xml2js = require('xml2js'),
    config = require('../default-config.js').config,
    ObjectID = require('mongodb').ObjectID;


computations.Samtools = function() {
    var self = this;

    this.view = function(bam_file, genomeId, start, end) { //the start and end parameters are optional.
        var samtools = exec.spawn("samtools", ["view", bam_file, start && end ? genomeId+":"+start+"-"+end : genomeId]),
            parser = new StringParser();
        parser.on('end', function(reads) {
            self.emit('got reads', reads);
        });
        samtools.stdout.setEncoding('utf8');
        samtools.stdout.on('data', function(data) { 
            parser.parseSam(data.trim());            
        });
        samtools.on('exit', function(code) {
            if (code != 0)
                self.emit('error', code);
            else 
                self.emit('end');
        });
    };
};

computations.Samtools.prototype.__proto__ = EventEmitter.prototype;

computations.SnoGPS = function() {
    var self = this,
        targetMolecules = undefined,
        targetMoleculesFilePath = undefined,
        installationPath = undefined;

    this.initialize = function(targetMolecules) {
        process.nextTick( function() {                    
            exec.exec("which pseudoU_test", function(err, stderr, stdout) {
                if (stderr.length == 0) 
                    self.emit('error', "pseudoU_test is not available in your PATH");
                else {
                    self.installationPath = path.normalize(stderr.split('pseudoU_test')[0] +"..");
                    tmp.file(function _tempFileCreated(err, filePath, fd) {
                        if (err)
                             self.emit('error', err);
                        else {
                            var fileWriter = new FileWriter();
                            fileWriter.on('end', function() {
                                self.targetMolecules = targetMolecules;
                                self.targetMoleculesFilePath = filePath;
                                self.emit('target molecules saved', filePath);              
                            });
                            fileWriter.writeFasta(targetMolecules, filePath);
                        }
                    });
                }
            });
        });
    };

    this.search = function(descriptor_file, scores_table_file, targets_file, threshold) {
        process.nextTick( function() {
            //console.log("pseudoU_test -L "+path.normalize(self.installationPath+"/scoretables/"+scores_table_file)+" -T "+path.normalize(self.installationPath+"/targs/"+targets_file)+" -F /dev/null "+self.targetMoleculesFilePath+" "+path.normalize(self.installationPath+"/desc/"+descriptor_file));
            var pseudo = exec.spawn("pseudoU_test", ["-L", path.normalize(self.installationPath+"/scoretables/"+scores_table_file), "-T", path.normalize(self.installationPath+"/targs/"+targets_file), "-F", "/dev/null", self.targetMoleculesFilePath, path.normalize(self.installationPath+"/desc/"+descriptor_file)]),
                content = "";
            pseudo.stdout.setEncoding('utf8');
            pseudo.stdout.on('data', function(data) { 
                content += data;
            });
            pseudo.on('exit', function(code) {
                if (code != 0)
                    self.emit('error', code);
                else {
                    var lines = content.split('\n'),
                        hits = [];
                    for (var i = 0 ; i < lines.length ; i++) {
                        var line = lines[i],
                            hit = {};
                        if (line.match(/^>/)) {                
                            line = line.replace(/\s+/g, ' ');
                            tokens = line.split(' ');  
                            for (var j = 0 ; j < tokens.length ; j++) { 
                                if (tokens[j].match(/^>/)) {
                                    hit.targetName = tokens[j].substring(1);
                                    hit.organism = targetMolecules[0].organism;
                                    hit.snoGPS_score = parseFloat(tokens[j+1]);
                                    if (hit.snoGPS_score >= (threshold || 25))
                                        hits.push(hit);
                                    hit.start_position = parseInt(tokens[j+2].split('-')[0].substring(1));
                                    hit.end_position = parseInt(tokens[j+2].split('-')[1].substring(0, tokens[j+2].trim().split('-')[1].length-1));
                                }
                                if (tokens[j].match(/^Cmpl:$/))
                                    hit.predicted_target_uridine = tokens[j+2];
                                if (tokens[j].match(/^Pairs:$/)) {
                                    hit.base_pairings = tokens[j+1],
                                    hit.snoRNA_length = tokens[j+2];
                                    if (tokens[j+4].match(/^\(W\)$/))
                                        hit.strand = "+";
                                    else
                                        hit.strand = "-";
                                }
                                    
                            }
                            var line2 = lines[i+1].replace(/\s+/g, ''),
                                line3 = lines[i+2].replace(/\s+/g, ' '),
                                line4 = lines[i+3].replace(/\s+/g, ' '),
                                line5 = lines[i+4].replace(/\s+/g, ' '),
                                line6 = lines[i+5].replace(/\s+/g, ' '),
                                line7 = lines[i+7].replace(/\s+/g, ' '),
                                line8 = lines[i+8].replace(/\s+/g, ' '),
                                line9 = lines[i+9].replace(/\s+/g, ' '),
                                line10 = lines[i+11].replace(/\s+/g, ' '),
                                tokens3 = line3.split(' '),
                                tokens4 = line4.split(' '),
                                tokens5 = line5.split(' '),
                                tokens6 = line6.split(' '),
                                tokens7 = line7.split(' '),
                                tokens8 = line8.split(' '),
                                tokens9 = line9.split(' '),
                                tokens10 = line10.split(' ');

                            hit.sequence = line2; //snoGPS is producing the correct sequence for the RNA

                            for (var k = 0 ; k < tokens3.length ; k++) { 
                                if (tokens3[k].match(/^Gap/)) {
                                    hit.gap_len = tokens3[k+2];
                                    hit.gap_sc = tokens3[k+4];
                                    hit.compl_len_score = tokens3[k+8];
                                    hit.intStem_len_sc = tokens3[k+12];
                                    hit.intStem_st_sc = tokens3[k+16];
                                }
                            }
                            for (var k = 0 ; k < tokens4.length ; k++) { 
                                if (tokens4[k].match(/^IntStem/)) {
                                    hit.intStem_dither = tokens4[k+2];
                                    hit.intStem_dither_score = tokens4[k+4];
                                    hit.extStem_dither = tokens4[k+7];
                                    hit.extStem_dither_score = tokens4[k+9];
                                    hit.HACA_RCstart_int = tokens4[k+12];
                                    hit.HACA_RCstart_sc = tokens4[k+14]; 
                                }
                            }
                            for (var k = 0 ; k < tokens5.length ; k++) { 
                                if (tokens5[k].match(/^HACA:/)) {
                                    hit.HACA_sequence = tokens5[k+1];
                                    hit.HACA_sequence_len = tokens5[k+2].substring(1, tokens5[k+2].length-1);
                                    hit.HACA_sequence_sc = tokens5[k+4]; 
                                }
                            }
                            for (var k = 0 ; k < tokens6.length ; k++) { 
                                if (tokens6[k].match(/^HACA-ISTEMRend/))
                                    hit.HACA_INSTEMRend_int = tokens6[k+2];  
                            }
                            for (var k = 0 ; k < tokens7.length ; k++) {
                                if (tokens7[k].match(/^LC/)) {
                                    hit.LC_sc = tokens7[k+2];
                                    hit.LC_sc_2 = tokens7[k+3].substring(1, tokens7[k+3].length-1);
                                    hit.RC_sc = tokens7[k+6];
                                    hit.RC_sc_2 = tokens7[k+7].substring(1, tokens7[k+7].length-1);
                                    hit.intStem_sc = tokens7[k+10];
                                    hit.intStem_sc_2 = tokens7[k+11].substring(1, tokens7[k+11].length-1);
                                    hit.extStem_sc = tokens7[k+14];
                                    hit.extStem_sc_2 = tokens7[k+15].substring(1, tokens7[k+15].length-1);
                                }
                            }
                            for (var k = 0 ; k < tokens8.length ; k++) {
                                if (tokens8[k].match(/^LCompl:/)) {
                                    hit.LCompl_sequence = tokens8[k+1].split('-')[0];
                                    hit.RComp_sequence = tokens8[k+3].split('-')[0];
                                    hit.intStem_sequence = tokens8[k+5].split("'")[1].substring(0, tokens8[k+5].split("'")[1].length-4);
                                    hit.extStem_sequence = tokens8[k+7].split("'")[1];
                                }
                            }
                            for (var k = 0 ; k < tokens9.length ; k++) {
                                if (tokens9[k].match(/^rRNA:/)) {
                                    hit.LrRNA_sequence = tokens9[k+1].split('-')[0].split('').reverse().join('');
                                    hit.RrRNA_sequence = tokens9[k+3].split('-')[0].split('').reverse().join(''); 
                                    hit.intStem_sequence_2 = tokens9[k+4].substring(0, tokens9[k+4].length-4).split('').reverse().join('');
                                    hit.extStem_sequence_2 = tokens9[k+6].split("'")[1].split('').reverse().join('');
                                }
                            }
                            hit.stem1_sequence = line10;
                        }
                    }                    
                    self.emit('hits',hits);
                }
            });
        });
    };
};

computations.SnoGPS.prototype.__proto__ = EventEmitter.prototype;

computations.Snoscan = function(targetMolecules) {
    var self = this,
        targetMolecules = undefined,
        targetMoleculesFilePath = undefined,
        installationPath = undefined;

    this.initialize = function(targetMolecules) {
        process.nextTick( function() {                    
            exec.exec("which snoscan", function(err, stderr, stdout) {
                if (stderr.length == 0) 
                    self.emit('error', "snoscan is not available in your PATH");
                else {
                    var tokens = stderr.split('/');
                    self.installationPath = path.normalize('/'+tokens.slice(1,tokens.length-1).join('/')+'/');
                    tmp.file(function _tempFileCreated(err, filePath, fd) {
                        if (err)
                             self.emit('error', err);
                        else {
                            var fileWriter = new FileWriter();
                            fileWriter.on('end', function() {
                                self.targetMolecules = targetMolecules;
                                self.targetMoleculesFilePath = filePath;
                                self.emit('target molecules saved', filePath);              
                            });
                            fileWriter.writeFasta(targetMolecules, filePath);
                        }
                    });
                }
            });
        });
    };

    this.search = function(methylation_sites_file, rRNA_file, threshold) {
        process.nextTick( function() {
            //console.log("snoscan "+["-m", path.normalize(self.installationPath+"/"+methylation_sites_file), path.normalize(self.installationPath+"/"+rRNA_file), self.targetMoleculesFilePath].join(' '));
            var snoscan = exec.spawn("snoscan", ["-m", path.normalize(self.installationPath+"/"+methylation_sites_file), path.normalize(self.installationPath+"/"+rRNA_file), self.targetMoleculesFilePath]),
                content = "";
            snoscan.stdout.setEncoding('utf8');
            snoscan.stdout.on('data', function(data) { 
                content += data;
            });
            snoscan.on('exit', function(code) {
                if (code != 0)
                    self.emit('error', code);
                else {
                    var lines = content.split('\n'),
                        hits = [];
                    for (var i = 0 ; i < lines.length ; i++) {
                        var line = lines[i],
                            tokens = line.split(' '),
                            hit = {};
                        if (line.match(/^>>/)) {               
                            for (var j = 0 ; j < tokens.length ; j++) {
                                if (tokens[j].match(/^>>/)) {
                                    hit.query_name = tokens[j+1];
                                    hit.organism = self.targetMolecules[0].organism;
                                    hit.score = parseFloat(tokens[j+3]);
                                    if (hit.score > (threshold || 25)) //min threshold to be determined
                                        hits.push(hit); 
                                    hit.start_position = parseInt(tokens[j+5].split('-')[0].substring(1));
                                    hit.end_position = parseInt(tokens[j+5].split('-')[1].substring(0, tokens[j+5].split('-')[1].length-1));     
                                    if (hit.start_position < hit.end_position)
                                        hit.strand = "+";
                                    else
                                        hit.strand = "-";
                                    //snoscan doesn't output the RNA sequence like snoGPS. We compute it...
                                    var nameRegexp = new RegExp("^"+hit.query_name);
                                    for (var n = 0 ; n < self.targetMolecules.length ; n++) {
                                        var targetMolecule = self.targetMolecules[n];
                                        if (nameRegexp.test(targetMolecule.name)) {
                                            if (hit.start_position < hit.end_position)
                                                hit.sequence = targetMolecule.sequence.substring(hit.start_position - 1, hit.end_position - 1);
                                            else
                                                hit.sequence = new DNA("fake_mol",targetMolecule.sequence.substring(hit.start_position - 1, hit.end_position - 1)).complement().split('').reverse().join(''); //we create a fake DNA molecule to compute the complement for what is REALLY necessary (needs less computing resources)
                                            break;
                                        }
                                    }
                                } 

                                if (tokens[j].match(/^Cmpl:/)) {
                                    var target = tokens[j+1].split('-');
                                    hit.target_name = target[0]+"-"+target[1];
                                    if (target[2]) {
                                        hit.methylation_site = target[2];
                                        hit.methylation_nucleotide = target[2].split('m')[0];
                                        hit.methylation_position = target[2].split('m')[1];
                                    }
                                }

                                if (tokens[j].match(/^Gs-DpBox:/)) {
                                    pairings = tokens[j-3].split('/');
                                    hit.pairings = pairings[0];
                                    hit.mismatches = pairings[1];
                                    hit.guide_region_adjacent_to = "Dprime Box"; 
                                    hit.guide_region_start_position_on_query = tokens[j+1];
                                    hit.guide_region_start_position_on_hit = tokens[j+2].substring(1, tokens[j+2].length-1);                                
                                }

                                if (tokens[j].match(/^box:/)) {
                                    pairings = tokens[j-4].split('/');
                                    hit.pairings = pairings[0];
                                    hit.mismatches = pairings[1];
                                    hit.guide_region_adjacent_to = "D Box";
                                    hit.guide_region_start_position_on_query = tokens[j+1];
                                    hit.guide_region_start_position_on_hit = tokens[j+2].substring(1, tokens[j+2].length-1);                                  
                                }

                                if (tokens[j].match(/^Len:/))
                                    hit.snoRNA_length = tokens[j+1];

                                if (lines[i].match(/TS$/))
                                    hit.terminal_stem = "present"
                                else
                                    hit.terminal_stem = "not present"
                            }
                                   
                            var line2 = lines[i+2],
                                line3 = lines[i+3],
                                line4 = lines[i+4],
                                line6 = lines[i+6],
                                line9 = lines[i+9],
                                line11 = lines[i+11],
                                line14 = lines[i+14],
                                line16 = lines[i+16],
                                line18 = lines[i+18],
                                tokens2 = line2.split(' '),
                                tokens3 = line3.split(' '),
                                tokens4 = line4.split(' '),
                                tokens6 = line6.split(' '),
                                tokens9 = line9.split(' '),
                                tokens11 = line11.split(' '),
                                tokens14 = line14.split(' '),
                                tokens16 = line16.split(' '),
                                tokens18 = line18.split(' ');

                            for (var l = 0 ; l < tokens2.length ; l++) {
                                hit.C_Box_sequence = tokens2[3];
                                hit.C_Box_score = tokens2[7];
                                hit.C_Box_start_position = tokens2[11].split('-')[0].substring(1);
                                hit.C_Box_end_position = tokens2[11].split('-')[1].substring(0, tokens2[11].trim().split('-')[1].length-1);
                                hit.C_D_box_dist = tokens2[16];
                            }

                            for (var l = 0 ; l < tokens3.length ; l++) {  
                                hit.D_Box_sequence = tokens3[3];
                                hit.D_Box_score = tokens3[10];
                            }

                            for (var l = 0 ; l < tokens4.length ; l++) {

                                if (line4.match(/[AUGC]+/)) {
                                    hit.Dprime_Box_sequence = tokens4[2];
                                    hit.Dprime_Box_score = tokens4[9];
                                }

                                if (line4.match(/None/)) {
                                    hit.Dprime_Box_sequence = "none";
                                    hit.Dprime_Box_score = "none";   
                                }

                                if (tokens4[l].match(/^Sc:/)) {
                                    if (line4.match(/None/)) {
                                        hit.D_Box_guide_transit_score = tokens4[l+1];
                                    } else {
                                        hit.Dprime_Box_guide_transit_score = tokens4[l+1];
                                    }
                                }
                            }

                            for (var l = 0 ; l < tokens6.length ; l++) {

                                if (line6.match(/^Meth/)) {
                                    hit.known_methylation_site_position = tokens6[3];
                                    hit.known_methylation_site_name = tokens6[4].substring(1, tokens6[4].length-1);
                                } else {
                                    hit.known_methylation_site_position = "none";
                                    hit.known_methylation_site_name = "none";
                                }

                                if (tokens6[l].match(/^Sc:/)) {
                                    hit.guide_sequence_score = tokens6[l+1];
                                    hit.guide_sequence_score_1 = tokens6[l+3].substring(1);
                                    hit.guide_sequence_score_2 = tokens6[l+4];
                                    hit.guide_sequence_score_3 = tokens6[l+5];
                                    hit.guide_sequence_score_4 = tokens6[l+6].substring(0, tokens6[l+6].length-1);
                                }
                            }

                            for (var l = 0 ; l < tokens9.length ; l++) {
                                hit.Db_seq = /[AUGC]{1,}/.exec(line9)[0];
                                hit.Db_seq_start_position = /\([0-9]{0,}-[0-9]{0,}\)/.exec(line9)[0].split('-')[0].substring(1);
                                hit.Db_seq_end_position = /\([0-9]{0,}-[0-9]{0,}\)/.exec(line9)[0].split('-')[1].substring(0, /\([0-9]{0,}-[0-9]{0,}\)/.exec(line9)[0].split('-')[1].length-1); 
                            }

                            for (var l = 0 ; l < tokens11.length ; l++) {
                                Qry_seq_reverse = /[AUGC]{1,}/.exec(line11)[0];
                                hit.Qry_seq = Qry_seq_reverse.split("").reverse().join("");
                                hit.Qry_seq_start_position = /\([0-9]{0,}-[0-9]{0,}\)/.exec(line11)[0].split('-')[1].substring(0, /\([0-9]{0,}-[0-9]{0,}\)/.exec(line11)[0].split('-')[1].length-1); 
                                hit.Qry_seq_end_position = /\([0-9]{0,}-[0-9]{0,}\)/.exec(line11)[0].split('-')[0].substring(1);
                            }

                            for (var l = 0 ; l < tokens14.length ; l++) {
                                hit.C_Box_gap_score = tokens14[6];
                                hit.C_Box_gap_length = tokens14[7].substring(1);
                                hit.D_Box_gap_score = tokens14[17];
                                if (D_Box_gap_length = tokens14[18]) {
                                    hit.D_Box_gap_length = D_Box_gap_length.substring(1);
                                }
                            }

                            for (var l = 0 ; l < tokens16.length ; l++) {
                                if (line16.match(/Sc:/)) {
                                    hit.stem_score = tokens16[l-2];
                                    if (stem_length = tokens16[l-1]) 
                                        hit.stem_length = stem_length.substring(1);
                                }
                            }

                            for (var l = 0 ; l < tokens18.length ; l++) {
                                if (line18.match(/Sc:/))
                                    hit.stem_transit_score = tokens18[l];
                            }
                        }  
                    }
                    self.emit('hits',hits);  
                }
            });
        });
    };
};

computations.Snoscan.prototype.__proto__ = EventEmitter.prototype;

computations.Rnaview = function() {
    var self = this;

    this.annotate = function(tertiaryStructure, rawData) {
        process.nextTick( function() {
            var fileWriter = new FileWriter(),
                filePath = config.data.root+"rnaview_input_"+(Math.floor(Math.random() * 90000) + 10000)+".pdb" ;
            fileWriter.on('error', function() {
                self.emit('error', "Problem to save the 3D in the file "+filePath);
            });
            fileWriter.on('end', function() {
                var rnaview = exec.spawn("rnaview", ["-p", filePath]),
                    content = "";
                rnaview.stdout.setEncoding('utf8');
                rnaview.stdout.on('data', function(data) { 
                    content += data;
                });
                rnaview.on('exit', function(code) {
                    if (code != 0) {
                        self.emit('error', code);
                    }
                    else {
                        var xmlParser = new xml2js.Parser();
                        path.exists(filePath+".xml", function(exists) {
                            if (exists) {
                                if (rawData || false) {
                                    fs.readFile(filePath+".xml", 'utf8', function(err, data) {
                                        self.emit('end',{'2D': data});
                                    });    
                                }
                                else {
                                    fs.readFile(filePath+".xml", function(err, data) {
                                        xmlParser.parseString(data, function (err, result) {

                                            var rna = new RNA(tertiaryStructure.rna.name, result.rnaml.molecule[0].sequence[0]['seq-data'][0].replace(/\s+/g,"")), // the rna molecule available from the XML output
                                                secondaryStructure = new SecondaryStructure(rna),
                                                newTertiaryStructure = undefined;
                                            secondaryStructure['source'] = "tool:rnaview:N.A.";
                                            
                                            if (rna.sequence.length != tertiaryStructure.rna.sequence.length) { //RNAVIEW had problems with some residues. Consequently, RNAVIEW has produced a new RNA molecule. We need to fit the 3D to this molecule.
                                                newTertiaryStructure = new TertiaryStructure(rna);
                                                newTertiaryStructure['source'] = "tool:rnaview:N.A.";
                                                var numbering_table = result.rnaml.molecule[0].sequence[0]['numbering-table'][0]['_'].replace(/\s{2,}/g," ").trim().split(' ');
                                                var previousAbsPos = 0, missingResidues = 0;
                                                //the strategy is the following:
                                                //- the numbering-table in the XML output stores the labels of the 3D residues used by RNAVIEW
                                                //- for each residue label, we recover its absolute position in the numbering system of the original 3D
                                                //- we're keeping track of the number of missing residues
                                                //- in the new 3D, the new absolute position = the originl absolute position - the missing residues
                                                //- in the new 3D, the numbering-system link the residue label to its new absolute position 
                                                numbering_table.forEach(function(residue_label) {
                                                    var absPos = parseInt(tertiaryStructure.numbering_system[residue_label]); //we get the absolute position in the original 3D according to the residue label in the numbering table
                                                    missingResidues += (absPos-previousAbsPos-1);
                                                    newTertiaryStructure.residues[absPos-missingResidues] = tertiaryStructure.residues[absPos]; //in the new tertiary structure, the new absPos is the previous one minus the missing residues
                                                    newTertiaryStructure.numbering_system[residue_label] = absPos-missingResidues;
                                                    previousAbsPos = absPos;
                                                });
                                            } else //no problem, then we can substitute the RNA of the 2D for the RNA of the 3D 
                                                secondaryStructure.rna = tertiaryStructure.rna;

                                            var helices = result.rnaml.molecule[0].structure[0].model[0]['str-annotation'][0].helix;
                                            helices.forEach(function(helix) {
                                                secondaryStructure.addHelix(helix['$'].id, parseInt(helix['base-id-5p'][0]['base-id'][0].position), parseInt(helix['base-id-3p'][0]['base-id'][0].position), parseInt(helix['length'][0]));    
                                            });
                                               
                                            var singleStrands = result.rnaml.molecule[0].structure[0].model[0]['str-annotation'][0]['single-strand'];
                                            singleStrands.forEach(function(singleStrand) {
                                                var end5 = parseInt(singleStrand.segment[0]['base-id-5p'][0]['base-id'][0].position),
                                                    end3 = parseInt(singleStrand.segment[0]['base-id-3p'][0]['base-id'][0].position);
                                                secondaryStructure.addSingleStrand(singleStrand.segment[0]['seg-name'][0], end5, end3-end5+1);    
                                            });
                                            
                                            var basePairs = result.rnaml.molecule[0].structure[0].model[0]['str-annotation'][0]['base-pair'];
                                            basePairs.forEach(function(basePair) {
                                                var edge1 = '(',
                                                    edge2 = ')';
                                                if (basePair['edge-5p'][0] == 'H')
                                                    edge1 = '[';
                                                else if (basePair['edge-5p'][0] == 'S')
                                                    edge1 = '{';
                                                else if (basePair['edge-5p'][0] == '!')
                                                    edge1 = '!';

                                                if (basePair['edge-3p'][0] == 'H')
                                                    edge2 = ']';
                                                else if (basePair['edge-3p'][0] == 'S')
                                                    edge2 = '}'; 
                                                else if (basePair['edge-3p'][0] == '!')
                                                    edge2 = '!';    
                                                secondaryStructure.addBasePair(basePair['bond-orientation'][0].toUpperCase(), edge1, edge2, parseInt(basePair['base-id-5p'][0]['base-id'][0].position), parseInt(basePair['base-id-3p'][0]['base-id'][0].position));    
                                            });
                                            if (newTertiaryStructure)
                                                self.emit('end',{'2D': secondaryStructure, '3D': newTertiaryStructure});    
                                            else
                                                self.emit('end',{'2D': secondaryStructure, '3D': tertiaryStructure});   
                                        });
                                        
                                    });
                                }
                            } else {
                                self.emit('error', "rnaview produced no xml for "+filePath+" "+tertiaryStructure.source);
                            }                            
                        });    
                    }    
                });
            });
            fileWriter.writePdb(tertiaryStructure, filePath, true); 
        });
    };

};

computations.Rnaview.prototype.__proto__ = EventEmitter.prototype;
    
computations.Contrafold = function() {
    var self = this;

    this.fold = function(rna, output) {
        tmp.file(function _tempFileCreated(err, file_path, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    exec.exec("cd "+path.dirname(file_path)+" ; contrafold predict " + file_path, function(err, stdout) {
                        if (err)
                            console.log(err);
                        var viennaData = "";
                        stdout.split("\n").forEach(function(line){
                            if (line != ">structure") {
                                viennaData = viennaData.concat(line+'\n');
                            }
                        });
                        if (output == "bpseq") {
                            var stringParser = new StringParser();
                            stringParser.on('end',function(secondaryStructures) {
                                var ss = secondaryStructures.pop();
                                self.emit('end',ss.toBpseq());
                            });
                            stringParser.parseVienna(viennaData);    
                        } else if (output == "vienna") {
                            self.emit('end',viennaData);
                        } else {
                            var stringParser = new StringParser();
                            stringParser.on('end',function(secondaryStructures) {
                                var ss = secondaryStructures.pop();
                                ss.rna = rna;
                                ss.source = "tool:contrafold:N.A.";
                                self.emit('end',ss);
                            });
                            stringParser.parseVienna(viennaData);
                        }
                    });
                });
            });
            fileWriter.writeFasta([rna],file_path,true);
        });
    }
};

computations.Contrafold.prototype.__proto__ = EventEmitter.prototype;

computations.Rnafold = function() {
    var self = this;
    
    this.fold = function(rna, output) {
        tmp.file(function _tempFileCreated(err, file_path, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    exec.exec("cd "+path.dirname(file_path)+" ; RNAfold < " + file_path, function(err, stdout) {
                        if (stdout.length != 0) {
                            var viennaData = "";
                            stdout.split("\n").forEach(function(line) {
                                var tokens = line.split(" ");
                                if (tokens.length >= 2) {
                                    viennaData = viennaData.concat(tokens[0]+'\n');
                                }
                                else {
                                    viennaData = viennaData.concat(line+'\n');
                                }
                            });
                            if (output == "bpseq") {
                                var stringParser = new StringParser();
                                stringParser.on('end',function(secondaryStructures) {
                                    var ss = secondaryStructures.pop();
                                    self.emit('end', ss.toBpseq());
                                });
                                stringParser.parseVienna(viennaData);    
                            } else if (output == "vienna") {
                                self.emit('end', viennaData);
                            } else {
                                var stringParser = new StringParser();
                                stringParser.on('end',function(secondaryStructures) {
                                    var ss = secondaryStructures.pop();
                                    ss.rna = rna;
                                    ss.source = "tool:rnafold:N.A.";
                                    self.emit('end', ss);
                                });
                                stringParser.parseVienna(viennaData);
                            }
                        }
                    });
                });
            });
            fileWriter.writeFasta([rna],file_path,true);
        });
    }
};

computations.Rnafold.prototype.__proto__ = EventEmitter.prototype;

computations.Rnaplot = function() {
    var self = this;
    
    this.plot = function(secondaryStructure, output) {
        tmp.file(function _tempFileCreated(err, file_path, fd) {
            var previousName = secondaryStructure.name;
            secondaryStructure.name = "2D";
            var fileWriter = new FileWriter(),
                coords2D = [];
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    exec.exec("cd "+path.dirname(file_path)+" ; RNAplot -o svg < " + file_path, function(err, stdout) {
                        if (stdout.length != 0) {
                            var inStream = fs.createReadStream(path.normalize(path.dirname(file_path)+"/2D_ss.svg"), {flags:'r'}),
                                svgOutput = "";
                            
                            carrier.carry(inStream).
                                on('line', function(line) {
                                    if (output == 'svg')
                                        svgOutput += line;
                                    else {   
                                        var regex = /<text x="(.+)" y="(.+)">.+/,
                                            coords = regex.exec(line.trim());
                                        if (coords)
                                            coords2D.push([parseFloat(coords[1]),parseFloat(coords[2])]);
                                    }
                                }).
                                on('end',function() {
                                    if (output == 'svg')
                                        self.emit('end',svgOutput);
                                    else {    
                                        secondaryStructure.plot = coords2D;
                                        self.emit('end',secondaryStructure);
                                    }
                                });
                        }
                    });
                });
            });
            fileWriter.writeVienna(secondaryStructure,file_path,true);
            secondaryStructure.name = previousName;
        });
    }
};

computations.Rnaplot.prototype.__proto__ = EventEmitter.prototype;

computations.Rnadistance = function() {
    var self = this;
    
    this.compute = function(bracketNotations) {
        tmp.file(function _tempFileCreated(err, file_path, fd) {
            var input="";
            bracketNotations.forEach(function(bn) {
                input+=(bn+"\n");
            });
            console.log(input);
            fs.writeFileSync(file_path, input); 
            process.nextTick( function() {
                exec.exec("RNAdistance -Xf -B < " + file_path, function(err, stdout) {
                    if (err)
                        self.emit('error', err);
                    else
                        self.emit('end', stdout);  
                });
            });
        });
    }
};

computations.Rnadistance.prototype.__proto__ = EventEmitter.prototype;

computations.RnaAlifold = function() {
    var self = this;
    
    this.getMostInformativeSequence = function(alignment) {
        tmp.file(function _tempFileCreated(err, file_path, fd) {
            fs.writeFileSync(file_path, alignment); 
            process.nextTick( function() {
                exec.exec("RNAalifold -mis " + file_path, function(err, stdout) {
                    if (err)
                        self.emit('error', err);
                    else
                        self.emit('end', stdout);  
                });
            });
        });
    };
};

computations.RnaAlifold.prototype.__proto__ = EventEmitter.prototype;

computations.Blast = function(targetMolecules) {
    var self = this;
    this.db = null;
    this.targetMolecules = targetMolecules || [];

    this.formatDb = function(isNucleotide) {
        tmp.dir(function _tempDirCreated(err, dirPath) {
            tmp.file(function _tempFileCreated(err, fastaFile, fd) {
                var fileWriter = new FileWriter();
                fileWriter.on('end', function() {
                    process.nextTick( function() {
                        exec.exec("cd "+path.dirname(fastaFile)+" ; formatdb -i "+fastaFile+" -p "+((isNucleotide || true) ? "F" : "T")+" -o", function(err, stderr, stdout) {
                            if (err)
                                console.log(err);
                            self.db = fastaFile;
                            self.emit("database formatted", fastaFile);
                        });
                    });
                });
                fileWriter.writeFasta(self.targetMolecules, fastaFile, false);
            });
        });
    };

    this.blastclust = function(molecules, isNucleotide) {
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    exec.exec("blastclust -p "+((isNucleotide || true) ? "F" : "T")+" -i " + fastaFile, function(err, stderr, stdout) {
                        if (err)
                            console.log(err);
                        var lines = stderr.split('\n'),
                            clusters = [];
                        for (var i = 0 ; i < lines.length ; i++) {
                            var line = lines[i].trim();
                            if (line.length != 0 && !line.match(/Start clustering/)) {
                                var sequence_ids = line.split(' '),
                                    cluster = [];
                                molecules.forEach(function(molecule) {
                                    for (var j = 0 ; j < sequence_ids.length ; j++) {
                                        var sequence_id = sequence_ids[j];
                                        if (new RegExp("^"+sequence_id).test(molecule.name)) {
                                            cluster.push(molecule);
                                            break;
                                        }
                                    };
                                });
                                clusters.push(cluster);
                            }
                        }
                        self.emit('end', clusters);
                    });
                });
            });
            fileWriter.writeFasta(molecules || [], fastaFile, false);
        });
    };

    this.blastn = function(queryMolecule) {
        if (this.db == null) {
            console.log("You have to format a database before to be able to blast");
            return;
        }
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    exec.exec("blastall -p blastn -d "+self.db+" -i " + fastaFile, function(err, stderr, stdout) {
                        if (err)
                            console.log(err);
                        //console.log(stderr);
                        var hits = self.parseOutput(stderr);
                        self.emit("hits", hits, queryMolecule);
                    });
                });
            });
            fileWriter.writeFasta([queryMolecule],fastaFile,false);
        });
    };

    this.parseOutput = function(output) {
        var lines = output.split('\n'),
            source = null,
            queryName = null,
            queryPositions = [],
            targetPositions = [],
            targetName = null,
            queryName = null,
            queryStrand = null,
            targetStrand = null,
            hits = [];

        for (var i = 0 ; i < lines.length ; i++) {
            line = lines[i].trim();
            if (line.match(/^BLAST/))
                source = line.toLowerCase();  
            else if (line.match(/^Query= /))
                queryName = line.split("Query= ")[1];
            else if (line.match(/^Score =/)) {
                if (queryPositions.length > 0 && targetPositions.length > 0) { // we have a previous hit to store
                    var hit = {};
                    hit.source = source;
                    if (targetStrand == '-') {
                        targetPositions = targetPositions.reverse();
                        hit.targetPositions = [targetPositions[0][0],targetPositions.pop()[1]];
                    }
                    else
                        hit.targetPositions = [targetPositions[0][0],targetPositions.pop()[1]];
                    
                    if (queryStrand == '-') {
                        queryPositions = queryPositions.reverse();
                        hit.queryPositions = [queryPositions[0][0],queryPositions.pop()[1]];
                    }
                    else
                        hit.queryPositions = [queryPositions[0][0],queryPositions.pop()[1]];    
                    hit.eValue = eValue;
                    hit.targetStrand = targetStrand;
                    hit.queryStrand = queryStrand;
                    hit.targetName = targetName;
                    hit.queryName = queryName;
                    for (var j = 0 ; j < self.targetMolecules.length ; j++) {
                        var targetMolecule = self.targetMolecules[j];
                        if (targetMolecule.name.match(new RegExp(hit.targetName))) {
                            if (targetStrand == '+')
                                hit.sequence = targetMolecule.sequence.substring(hit.targetPositions[0]-1,hit.targetPositions[1]);  
                            else
                                hit.sequence = targetMolecule.complement().substring(hit.targetPositions[0]-1,hit.targetPositions[1]).split('').reverse().join('');
                            break;
                        }
                    }
                    hits.push(hit);
                }
                queryPositions = [];
                targetPositions = [];
                eValue = parseFloat(line.split("Expect =")[1].trim());
		
                if (lines[i-3].trim().length > 0 ) {//otherwise its a hit for the same sequence
                    targetName = lines[i-3].trim().substring(1);
                }
            } else if (line.match(/^Strand =/)) {
                var strandOrientations = line.split("Strand = ")[1].trim().split('/');
                queryStrand = strandOrientations[0].trim() == "Plus" ? '+' : '-';
                targetStrand = strandOrientations[1].trim() == "Plus" ? '+' : '-';
            } else if (line.match(/^Query:/)) {
                var tokens = line.split(' '),
                    fragmentStart = parseInt(tokens[1]),
                    fragmentEnd = parseInt(tokens.pop());
                if (fragmentStart < fragmentEnd)
                    queryPositions.push([fragmentStart,fragmentEnd]);
                else
                    queryPositions.push([fragmentEnd,fragmentStart]);
            } else if (line.match(/^Sbjct:/)) {
                var tokens = line.split(' '),
                    fragmentStart = parseInt(tokens[1]),
                    fragmentEnd = parseInt(tokens.pop());
                if (fragmentStart < fragmentEnd)
                    targetPositions.push([fragmentStart,fragmentEnd]);
                else
                    targetPositions.push([fragmentEnd,fragmentStart]);
            }
        }
        if (queryPositions.length > 0 && targetPositions.length > 0) { // the last hit to store
            var hit = {};
            hit.source = source;
            if (targetStrand == '-') {
                targetPositions = targetPositions.reverse();
                hit.targetPositions = [targetPositions[0][0],targetPositions.pop()[1]];
            }
            else
                hit.targetPositions = [targetPositions[0][0],targetPositions.pop()[1]];
            
            if (queryStrand == '-') {
                queryPositions = queryPositions.reverse();
                hit.queryPositions = [queryPositions[0][0],queryPositions.pop()[1]];
            }
            else
                hit.queryPositions = [queryPositions[0][0],queryPositions.pop()[1]]; 
            hit.eValue = eValue;
            hit.targetStrand = targetStrand;
            hit.queryStrand = queryStrand;
            hit.targetName = targetName;
            hit.queryName= queryName;
            for (var j = 0 ; j < self.targetMolecules.length ; j++) {
                var targetMolecule = self.targetMolecules[j];
                if (targetMolecule.name.match(new RegExp(hit.targetName))) {
                    if (targetStrand == '+')
                        hit.sequence = targetMolecule.sequence.substring(hit.targetPositions[0]-1,hit.targetPositions[1]);  
                    else
                        hit.sequence = targetMolecule.complement().substring(hit.targetPositions[0]-1,hit.targetPositions[1]).split('').reverse().join('');
                    break;
                }
            }
            hits.push(hit);
        }
        return hits;
    };

};

computations.Blast.prototype.__proto__ = EventEmitter.prototype;

computations.Clustalw = function() {

    var self = this;

    this.align = function(molecules) {
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    tmp.file(function _tempFileCreated(err, outputFile, fd) {
                        exec.exec("clustalw2 -infile="+fastaFile+" -outfile="+outputFile, function(err, stderr, stdout) {
                            if (err)
                                self.emit('error', err);
                            else {
                                var inStream = fs.createReadStream(outputFile, {flags:'r'}),
                                    alignment = "";
                                carrier.carry(inStream).
                                    on('line', function(line) {
                                        alignment += line+"\n";
                                    }).
                                    on('end',function() {
                                        self.emit('end', alignment);
                                    });
                            }
                        });
                    });
                });
            });
            fileWriter.writeFasta(molecules || [], fastaFile, false);
        });        
    };

};

computations.Clustalw.prototype.__proto__ = EventEmitter.prototype;

computations.Cmsearch = function(targetMolecules) {
    var self = this;
    this.targetMolecules = targetMolecules || [];

    this.search = function(rfamId, rfamRoot) {
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    exec.exec("cmsearch --ga " + path.normalize(rfamRoot+"/CMs/"+rfamId+".cm") + " "+fastaFile, function(err, stderr, stdout) {
                        if (err)
                            console.log(err);
                        self.emit('hits',self.parseOutput(stderr), rfamId);
                    });
                });
            });
            fileWriter.writeFasta(self.targetMolecules,fastaFile,false);
        });
    };

    this.parseOutput = function(output) {
        var lines = output.split('\n'),
            line = null,
            plusStrand = true,
            hit = null,
            source = null,
            cm = null,
            rfamId = null,
            targetStrand = null,
            targetSequence = null,
            targetName = null,
            hits = [];
        for (var i = 0 ; i < lines.length ; i++) {
            line = lines[i].trim();
            if (line.match(/^# INFERNAL/))
                source = "tool:cmsearch:N.A.";
            else if (line.match(/^# command:/))
                cm = line.split("cmsearch --ga")[1].split(' ')[1].split('/').pop();
            else if (line.match(/^CM:/))
                rfamId = line.split("CM:")[1].trim();
            else if (line.match(/^Plus strand results:/)) {
                targetStrand = '+';
                if (lines[i-2].trim().length != 0) { //otherwise its Plus Strand results after Minus Strand results for the same sequence
                    targetName = lines[i-2].substring(1).trim();
                }
            } else if (line.match(/^Minus strand results:/)) {
                targetStrand = '-';
                if (lines[i-2].trim().length != 0) { //otherwise its Minus Strand results after Plus Strand results for the same sequence
                    targetName = lines[i-2].substring(1).trim();
                }
            } else if (line.match(/^Query =/)) {
                var tokens = line.split(", Target = ")[0].trim().split("Query = ")[1].split(" - "),
                    queryStart = parseInt(tokens[0]),
                    queryEnd = parseInt(tokens[1]);

                tokens = line.split(", Target = ")[1].trim().split(" - ");

                var targetStart = parseInt(tokens[0]),
                    targetEnd = parseInt(tokens[1]);

                i++;

                var scores = lines[i].split(", "),
                    score = null,
                    pValue = null,
                    eValue = null;

                scores.forEach(function(s) {
                    if (s.trim().match(/^Score/)) {
                        score = parseFloat(s.trim().split(" = ")[1].trim());
                    } else if (s.match(/^E/)) {
                        eValue = parseFloat(s.trim().split(" = ")[1].trim());  
                    } else if (s.match(/^P/)) {
                        pValue = parseFloat(s.trim().split(" = ")[1].trim());  
                    }
                });
                    
                i+=3;

                var queryPositions = [],
                    targetPositions = [],
                    currentQueryStart = queryStart,
                    currentTargetStart = targetStart,
                    querySequence = "",
                    targetSequence = "",
                    targetTurn = false ; //we can have situation where the target and query lines have the same start position

                while (true) {
                    if (!targetTurn) {
                        tokens = lines[i].trim().split(' ');
                        querySequence += tokens.splice(1,tokens.length-2).join('');
                        targetTurn = true;
                        i+=2;
                    } else {
                        tokens = lines[i].trim().split(' ');
                        targetSequence += tokens.splice(1,tokens.length-2).join('');
                        targetTurn = false;
                        if (lines[i].trim().match(new RegExp(targetEnd+'$')) || lines[i+2].trim().length == 0) {
                            break;
                        }
                        i += 3;
                    }
                }

                //console.log("querySequence: "+querySequence);
                //console.log("targetSequence: "+targetSequence);

                currentQueryStart = queryStart;
                currentTargetStart = targetStart;
                subsequences = querySequence.split('*');

                if (subsequences.length == 1) {
                    var residues = subsequences[0].length-(subsequences[0].split('.').length - 1)-(subsequences[0].split('-').length - 1);
                    queryPositions.push([currentQueryStart,currentQueryStart + residues-1]);
                } else {
                    subsequences.forEach(function(subsequence) {
                        if (subsequence.match(/^\[/))
                            currentQueryStart += parseInt(subsequence.substring(1,subsequence.length-1).trim());
                        else {
                            residues = subsequence.length-(subsequence.split('.').length - 1)-(subsequence.split('-').length - 1);
                            queryPositions.push([currentQueryStart,currentQueryStart+residues-1]);
                            currentQueryStart += residues
                        }
                    });
                }

                subsequences = targetSequence.split('*');

                if (subsequences.length == 1) {
                    var residues = subsequences[0].length-(subsequences[0].split('.').length - 1)-(subsequences[0].split('-').length - 1);
                    if (targetStrand == '+')
                        targetPositions.push([currentTargetStart, currentTargetStart + residues-1]);
                    else
                        targetPositions.push([currentTargetStart, currentTargetStart - residues+1]);
                } else {
                    subsequences.forEach(function(subsequence) {
                        if (subsequence.match(/^\[/)) {
                            if (targetStrand == '+')
                                currentTargetStart += parseInt(subsequence.substring(1,subsequence.length-1).trim());
                            else
                                currentTargetStart -= parseInt(subsequence.substring(1,subsequence.length-1).trim());
                        }
                        else {
                            residues = subsequence.length-(subsequence.split('.').length - 1)-(subsequence.split('-').length - 1);
                            if (targetStrand == '+') {
                                targetPositions.push([currentTargetStart,currentTargetStart+residues-1]);
                                currentTargetStart += residues;
                            }
                            else {
                                targetPositions.push([currentTargetStart,currentTargetStart-residues+1])
                                currentTargetStart -= residues;
                            }
                        }
                    });
                }

                if (eValue != null) {
                    hit = {};
                    hit.source = source;
                    hit.rfamId = rfamId;
                    hit.cm = cm; 
                    hit.targetStrand = targetStrand;
                    hit.targetSequence = targetSequence;
                    hit.targetName = targetName;
                    if (score != null)
                        hit.score = score;
                    if (eValue != null)
                        hit.eValue = eValue;
                    if (pValue != null)
                        hit.pValue = pValue;
                    hit.targetPositions = targetPositions[0][0] < targetPositions[targetPositions.length-1][1] ? [targetPositions[0][0],targetPositions[targetPositions.length-1][1]] : [targetPositions[targetPositions.length-1][1], targetPositions[0][0]];
                    hit.queryPositions = queryPositions[0][0] < queryPositions[queryPositions.length-1][1] ? [queryPositions[0][0],queryPositions[queryPositions.length-1][1]] : [queryPositions[queryPositions.length-1][1], queryPositions[0][0]];
                    for (var j = 0 ; j < self.targetMolecules.length ; j++) {
                        var targetMolecule = self.targetMolecules[j];
                        if (targetMolecule.name.match(new RegExp(hit.targetName))) {
                            if (hit.targetStrand == '+')
                                hit.sequence = targetMolecule.sequence.substring(hit.targetPositions[0]-1,hit.targetPositions[1]);  
                            else
                                hit.sequence = targetMolecule.complement().substring(hit.targetPositions[0]-1,hit.targetPositions[1]).split('').reverse().join('');
                            break;
                        }
                        
                    }
                    hits.push(hit);
                } else {
                    console.log("No e-value found. The covariance model "+cm+" was probably not calibrated");
                    break; //the parsing is stopped 
                }

            }
        }
        return hits;
    };
};

computations.Cmsearch.prototype.__proto__ = EventEmitter.prototype;

computations.Mlocarna = function() {

    var self = this;

    this.align = function(molecules) {
        var alignedMolecules = {};
        molecules.forEach(function(molecule) {
            alignedMolecules[molecule.name] = "";
        });
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    tmp.file(function _tempFileCreated(err, outputFile, fd) {
                        exec.exec("mlocarna "+fastaFile, function(err, stderr, stdout) {
                            if (err)
                                self.emit('error', err);
                            else {
                                var basePairs = [];
                                stderr.split('\n').forEach(function(line) {
                                    var tokens = line.split(/\s+/);

                                    if (tokens.length == 2 && alignedMolecules.hasOwnProperty(tokens[0]))
                                        alignedMolecules[tokens[0]] =  alignedMolecules[tokens[0]]+tokens[1];
                                    else if (tokens[0] == 'alifold') {
                                        var bn = tokens[1],
                                            i = 0,
                                            lastPairedPos = [];

                                        bn.split('').forEach(function(pos) {
                                            i++;
                                            if (pos == '(')
                                                lastPairedPos.push(i);
                                            else if (pos == ')')
                                                basePairs.push(new BaseBaseInteraction('c', '(',')', lastPairedPos.pop(), i));
                                        });
                                    }

                                });
                                
                                var _alignedMolecules = [];
                                for (var key in alignedMolecules)
                                    _alignedMolecules.push(new RNA(key, alignedMolecules[key]));

                                self.emit('end', _alignedMolecules, basePairs);
                            }
                        });
                    });
                });
            });
            fileWriter.writeFasta(molecules || [], fastaFile, false);
        });    
    };

};

computations.Mlocarna.prototype.__proto__ = EventEmitter.prototype;

computations.Foldalign = function() {

    var self = this;

    this.align = function(molecules) {
        var alignedMolecules = {};
        molecules.forEach(function(molecule) {
            alignedMolecules[molecule.name] = "";
        });
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    tmp.file(function _tempFileCreated(err, outputFile, fd) {
                        exec.exec("foldalign "+fastaFile, function(err, stderr, stdout) {
                            if (err)
                                self.emit('error', err);
                            else {
                                console.log(stderr);
                            }
                        });
                    });
                });
            });
            fileWriter.writeFasta(molecules || [], fastaFile, false);
        });    
    };

};

computations.Foldalign.prototype.__proto__ = EventEmitter.prototype;

computations.Cmalign= function() {

    var self = this;

    this.align = function(molecules, rfamid) {
        tmp.file(function _tempFileCreated(err, fastaFile, fd) {
            var fileWriter = new FileWriter();
            fileWriter.on('end', function() {
                process.nextTick( function() {
                    tmp.file(function _tempFileCreated(err, outputFile, fd) {
                        exec.exec("cmalign.py -f "+fastaFile+" -rfamid "+rfamid, function(err, stderr, stdout) {
                            if (err)
                                self.emit('error', err);
                            else
                                self.emit('end', stdout);
                        });
                    });
                });
            });
            fileWriter.writeFasta(molecules || [], fastaFile, false);    
        });
    };
};

computations.Cmalign.prototype.__proto__ = EventEmitter.prototype;


