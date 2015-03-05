/**
A module giving access to several public databases (RFAM, NCBI, Ensembl,...)
@module db
**/

var db = exports,
    parsers = require('./parsers'),
    EventEmitter = require('events').EventEmitter
    exec = require('child_process'),
    tmp = require('tmp'),
    path = require('path'),
    url = require('url');
    carrier = require('carrier'),
    fs = require('fs'),
    config = require('../default-config.js').config,
    mongo = require('mongoskin'),
    jsdom = require("jsdom"),
    RNA = require('../core/molecules').RNA,
    ftpClient = require('ftp'),
    rna = require('../rna.js');

db.Ensembl = function(ensemblRoot) {
    var self = this;
    this.ensemblRoot = path.normalize(ensemblRoot || config.data.root+"/ensembl");

    this.getAllNcRNAs = function() {
        exec.exec(path.normalize(module['filename']+"/../../../files/scripts/get_ensembl_ncrnas.sh")+ " " + self.ensemblRoot, function(err, stderr, stdout) {
            if (stdout.length != 0) {
                console.log(stdout);
            }
            if (stderr.length != 0) {
                console.log(stderr);
            }
            if (err) {
                console.log(err);
                self.emit('error');
            }
            fs.readdir(self.ensemblRoot, function(err, files) {
                var ncRNAs_files = [],
                    parsed_Files = 0;
                for (var h = 0 ; h < files.length ; h++) {
                    if (files[h].match(/ncrna.fa$/))
                        ncRNAs_files.push(files[h]);
                }
                for (var h = 0 ; h < ncRNAs_files.length ; h++) {
                    var file = files[h],
                        parser = new parsers.FileParser();
                    parser.on('end', function(ncRNAs) {
                        parsed_Files++;
                        for (var i = 0 ; i < ncRNAs.length ; i++) {
                            var ncRNA = ncRNAs[i],
                                tokens = ncRNA.name.split(' ');
                            ncRNA.organism = ncRNA.source.split(':')[2].split('.')[0].split('_').join(' ');
                            for (var j = 0 ; j < tokens.length ; j++) {
                                var token = tokens[j];
                                if (token.match("^gene_biotype"))
                                    ncRNA.family = token.split(":")[1];
                                else if (token.match('^gene'))
                                    ncRNA.source = 'db:ensembl:'+token.split(":")[1];
                                else if (token.match("^(chromosome|scaffold|ultracontig|reftig)")) {
                                    var _tokens = token.split(':'),
                                        strand = _tokens.pop();
                                    if (strand == '-1')
                                      ncRNA.genomicStrand = '-';
                                    else if (strand == '1')
                                      ncRNA.genomicStrand = '+';
                                    ncRNA.genomicPositions = [parseInt(_tokens.pop()),parseInt(_tokens.pop())].reverse();
                                }
                                else if (token.match("^transcript_biotype")) {
                                    //do nothing
                                }
                                else if (token.match("^ncrna")) {
                                    //do nothing
                                }
                                /*else
                                    console.log(token)*/
                            }
                        }
                        self.emit('entry parsed', ncRNAs);
                        if (parsed_Files == ncRNAs_files.length)
                            self.emit('end');
                    });
                    parser.parseFasta(path.normalize(self.ensemblRoot+'/'+file, 'RNA'));
                }
            });
            
        });      
    };

};

db.Ensembl.prototype.__proto__ = EventEmitter.prototype;

db.PDB = function() {
    var self = this;

    this.getEntry = function(pdbId) {
        var options = {
                host: 'www.rcsb.org',
                port: 80,
                path: '/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId='+pdbId
            },
            pdbContent = "";

        require('http').get(options, function(res) {
            res.on('data', function(data) {
                pdbContent = pdbContent.concat(data.toString());    
            }).on('end', function() {
                var pdbParser = new parsers.StringParser();
                pdbParser.on('end', function(tertiaryStructures) {
                    for (var i = 0 ; i < tertiaryStructures.length ; i++) {
                        tertiaryStructures[i].source = "db:pdb:"+pdbId;
                        tertiaryStructures[i].rna.source = "db:pdb:"+pdbId;  
                    }
                    self.emit('entry parsed',tertiaryStructures);
                });
                pdbParser.on('error', function(error) {
                    self.emit('error', pdbId+": "+error);
                });
                pdbParser.parsePDB(pdbContent);    
            });
        }).on('error', function(e) {
            console.log("Got error: " + e.message);
        });
    };

    this.getHeader = function(pdbId) {
        var options = {
                host: 'www.rcsb.org',
                port: 80,
                path: '/pdb/files/'+pdbId+'.pdb?headerOnly=YES'
            },
            headerData  = "",
            header= {'pdbId':pdbId};

        require('http').get(options, function(res) {
            res.on('data', function(data) {
                headerData  = headerData.concat(data.toString());    
            }).on('end', function() {
                var lines = headerData.split('\n'),
                    title = "",
                    authors = "",
                    date = "";
                for (var i = 0 ; i < lines.length ; i++) {
                    var line = lines[i];
                    if (line.substring(0,6).trim() == "HEADER") {
                        date = line.substring(50,59).trim()    
                    } else if (line.substring(0,6).trim() == "TITLE") {
                        title = title.concat(line.substring(10,70)).trim()+" ";
                    } else if (line.substring(0,6).trim() == "AUTHOR") {
                        authors = authors.concat(line.substring(10,70)).trim()+" ";
                    }
                }
                header['title'] = title.trim();
                header['date'] = date.trim();
                header['authors'] = authors.trim();
                self.emit('header parsed',header);
            });
        }).on('error', function(e) {
            console.log("Got error: " + e.message);
        });
    }

    this.query = function(query) {
        var minRes = query.minRes || null,
            maxRes = query.maxRes || null,
            minDate = query.minDate || null,
            maxDate = query.maxDate || null,
            keywords = query.keywords || [],
            authors = query.authors || [],
            pdbIds = query.pdbIds || [],
            titleContains = query.titleContains || [],
            containsRNA = query.containsRNA || '?',
            containsProtein = query.containsProtein || '?',
            containsDNA = query.containsDNA || '?',
            containsHybrid = query.containsHybrid || '?',
            experimentalMethod = query.experimentalMethod || null,
            post_data = '<orgPdbCompositeQuery version="1.0">',
            refinementLevel = 0,
            ids = "";

        if (maxRes != null || minRes != null) {
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.ResolutionQuery</queryType>\
<description>Resolution query</description>\
<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>');
            if (minRes != null) {
                post_data = post_data.concat('\
<refine.ls_d_res_high.min>'+minRes+'</refine.ls_d_res_high.min>');                   
            }
            if (maxRes != null) {
                post_data = post_data.concat('\
<refine.ls_d_res_high.max>'+maxRes+'</refine.ls_d_res_high.max>');
            }
            post_data = post_data.concat('</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        if (maxDate != null || minDate != null) {
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.ReleaseDateQuery</queryType>\
<description>Release Date query</description>\
<refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>');
            if (minDate != null) {
                post_data = post_data.concat('\
<database_PDB_rev.date.min>'+minDate+'</database_PDB_rev.date.min>');                   
            }
            if (maxDate != null) {
                post_data = post_data.concat('\
<database_PDB_rev.date.max>'+maxDate+'</database_PDB_rev.date.max>');
            }
            post_data = post_data.concat('</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        for (var i = 0 ; i < titleContains.length ; i++) {
            var titleContain = titleContains[i];
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.StructTitleQuery</queryType>\
<description>StructTitleQuery: struct.title.comparator=contains struct.title.value='+titleContain+'</description>\
<struct.title.comparator>contains</struct.title.comparator>\
<struct.title.value>'+titleContain+'</struct.title.value>\
</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        if (keywords.length != 0) {
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.AdvancedKeywordQuery</queryType>\
<description>Text Search for: '+keywords.join(" ")+'</description>\
<keywords>'+keywords.join(" ")+'</keywords>\
</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        if (pdbIds.length != 0) {
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.StructureIdQuery</queryType>\
<description>Simple query for a list of PDB IDs ('+pdbIds.length+' IDs) :'+pdbIds.join(", ")+'</description>\
<structureIdList>'+pdbIds.join(", ")+'</structureIdList>\
</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        if (experimentalMethod != null) {
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.ExpTypeQuery</queryType>\
<description>Experimental Method is '+experimentalMethod+'</description>\
<mvStructure.expMethod.value>'+experimentalMethod+'</mvStructure.expMethod.value>\
</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        for (var i = 0 ; i < authors.length ; i++) {
            var author = authors[i];
            if (refinementLevel > 0) {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
            } else {
                post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
            }
            post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.AdvancedAuthorQuery</queryType>\
<description>Author Search: Author Search: audit_author.name='+author+' OR (citation_author.name='+author+' AND citation_author.citation_id=primary)</description>\
<exactMatch>false</exactMatch>\
<audit_author.name>'+author+'</audit_author.name>\
</orgPdbQuery></queryRefinement>');
            refinementLevel++;
        }

        //chain type
        if (refinementLevel > 0) {
            post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel><conjunctionType>and</conjunctionType>');
        } else {
            post_data = post_data.concat('<queryRefinement><queryRefinementLevel>'+refinementLevel+'</queryRefinementLevel>');   
        }
        post_data = post_data.concat('\
<orgPdbQuery>\
<version>head</version>\
<queryType>org.pdb.query.simple.ChainTypeQuery</queryType>\
<description>Chain Type</description>\
<containsProtein>'+containsProtein+'</containsProtein>\
<containsDna>'+containsDNA+'</containsDna>\
<containsRna>'+containsRNA+'</containsRna>\
<containsHybrid>'+containsHybrid+'</containsHybrid>\
</orgPdbQuery></queryRefinement>');
        refinementLevel++;

        post_data = post_data.concat('</orgPdbCompositeQuery>');

        var options = {
                host: 'www.rcsb.org',
                port: 80,
                method: 'POST',
                path: '/pdb/rest/search',
                headers: {
                    'Content-Type': 'application/x-www-form-urlencoded',
                    'Content-Length': post_data.length
                }
            };

        var post_req = require('http').request(options, function(res) {
                res.setEncoding('utf8');
                res.on('data', function (data) {
                    ids = ids.concat(data);
                }).on('end', function() {
                    query["hits"] = ids.split('\n').slice(0,-1)
                    self.emit('end search',query);
                });
            });

        // post the data
        post_req.write(post_data);
        post_req.end();

    };
};

db.PDB.prototype.__proto__ = EventEmitter.prototype;

db.NCBI = function() {
    var self = this;

    var listFungalGenomes = function() {
        var conn = new ftpClient({host:'ftp.ncbi.nlm.nih.gov'}),
            fungalGenomes = [];
        conn.on('connect', function() {
            conn.auth(function(e) {
                if (e)
                    console.log(e);
                conn.list('genomes/Fungi/', function(e, eventEmiter) {
                    if (e)
                        console.log(e);
                    eventEmiter.on('entry', function(entry) {
                        if (entry.type === 'd')
                            fungalGenomes.push(entry.name); 
                    });
                    eventEmiter.on('end', function() {
                        self.emit('end genomes list',fungalGenomes);
                        conn.end();
                    });
                });
            });
        });
        conn.connect();
    };

    var listFungalGenomeFiles = function(fungalGenome) {
        var conn = new ftpClient({host:'ftp.ncbi.nlm.nih.gov'}),
            files = [];
        conn.on('connect', function() {
            conn.auth(function(e) {
                if (e)
                    console.log(e);
                conn.list('genomes/Fungi/'+fungalGenome, function(e, eventEmiter) {
                    if (e)
                        console.log(e);
                    eventEmiter.on('entry', function(entry) {
                        if (entry.type === '-')
                            files.push('genomes/Fungi/'+fungalGenome+'/'+entry.name);
                    });
                    eventEmiter.on('end', function() {
                        self.emit('end files list',files);
                        conn.end();
                    });
                });
            });
        });
        conn.connect();
    };

    self.getFungalGenomicFiles = function() {
        var genomes = {};

        self.on('end genomes list', function(fungalGenomes) {
            var i = 0;

            self.on('end files list', function(files) {
                genomes[fungalGenomes[i-1]] = files;
                if (i < fungalGenomes.length)
                    listFungalGenomeFiles(fungalGenomes[i++]);
                else 
                    self.emit('got fungal genomic files', genomes);          
            });

            listFungalGenomeFiles(fungalGenomes[i++]);
        });

        listFungalGenomes();
    };

    self.getFtpFileContent = function(filePath) {
        console.log(filePath);
        var conn = new ftpClient({host:'ftp.ncbi.nlm.nih.gov'}),
            content = "";
        conn.on('connect', function() {
            conn.auth(function(e) {
                if (e)
                    console.log(e);
                conn.get(filePath, function(e, stream) {
                    if (e)
                        console.log(e);
                    stream.setEncoding('utf8');
                    stream.on('success', function() {
                        self.emit('got ftp file content',content);
                        conn.end();
                    });
                    stream.on('error', function(e) {
                        console.log(e);
                        conn.end();
                    });
                    stream.on('data', function(data) {
                        content += data;
                    });
                });
            });
        });
        conn.connect();
    };
};

db.NCBI.prototype.__proto__ = EventEmitter.prototype;

db.Eutils = function() {
    var self = this;
    this.options = {
        host: 'eutils.ncbi.nlm.nih.gov',
        port: 80
    };

    self.fetch = function(options) {
        var db = options.db,
            id = options.id,
            rettype = options.rettype || "asn1",
            retmode = options.retmode|| "text",
            content = "";

        if (db == null) {
           console.log('Entrez database is missing ("db" key in params).');
           return;  
        }

        if (id == null) {
           console.log('Entrez UIDs are missing ("id" key in params).');
           return;  
        }

        self.options.path = "/entrez/eutils/efetch.fcgi?db="+db+"&id="+id.join(',')+"&rettype="+rettype+"&retmode="+retmode;

        require('http').get(self.options, function(res) {
            res.on('data', function(data) {
                content += data;  
            }).on('end', function() {
                if (content.match(/Temporarily Unavailable/)) {
                    self.emit('error',"Service Temporarily Unavailable");
                }
                else if (rettype == "fasta") {
                    var parser = new parsers.StringParser();
                    parser.on("end", function(molecules) {
                        self.emit('end fetch',molecules);    
                    });
                    if (db == "nucleotide") {
                        parser.parseFasta(content,'DNA');
                    }
                }
                else if (rettype == "gb") {
                    var parser = new parsers.StringParser();
                    parser.on("end", function(data) {
                        if (!data.organism)
                            console.log(content); 
                        self.emit('end fetch',data);    
                    });
                    parser.parseGenbank(content);
                }
            });
        });

    };
};

db.Eutils.prototype.__proto__ = EventEmitter.prototype;

db.RnajsDb = function(database, host, port) {
    var self = this;
    this.db = mongo.db((host || "localhost") +':'+ (port || 27017) +'/'+database, {safe:true});

    //emit all the computation objects ids recorded in MongoDB associated to an executable name. 
    this.getExecutables = function(organismName) {
        self.db.collection('genomes').find().toArray(function(err, genomes) {
            var organismNames2GenomesIds = {};
            if (err)
                console.log(err);
            else if (genomes) {
                genomes.forEach(function(genome) {
                    if (organismNames2GenomesIds.hasOwnProperty(genome.organism))
                        organismNames2GenomesIds[genome.organism].push(genome._id.toString());
                    else
                        organismNames2GenomesIds[genome.organism] = [genome._id.toString()];       
                });
                self.db.collection('computations').find().toArray(function(err, computations) {
                    if (err)
                        console.log(err);
                    else if (computations) {
                        var executables = {}, computationsProcessed = 0 ;
                        for (var j=0 ; j < computations.length ; j++) {
                            var computation = computations[j],
                                inputs = computation['inputs'];
                            for (var i=0 ; i < inputs.length ; i++) {
                                if (inputs[i].match(/@genomes/) && organismNames2GenomesIds[organismName].indexOf(inputs[i].split('@genomes')[0]) != -1) {
                                    if (!executables.hasOwnProperty(computation.executable))
                                        executables[computation.executable] = [computation._id];
                                    else
                                       executables[computation.executable].push(computation._id);
                                    break;
                                }
                            }
                            computationsProcessed++;
                            if ( computationsProcessed == computations.length) {
                                self.emit("got executables", executables); 
                                break;
                            }
                        };
                    }
                });   
            }
        });      
    };

    this.getComputations = function(executableName, organismName) {
        self.db.collection('genomes').find().toArray(function(err, genomes) {
            var organismNames2GenomesIds = {};
            if (err)
                console.log(err);
            else if (genomes) {
                genomes.forEach(function(genome) {
                    if (organismNames2GenomesIds.hasOwnProperty(genome.organism))
                        organismNames2GenomesIds[genome.organism].push(genome._id.toString());
                    else
                        organismNames2GenomesIds[genome.organism] = [genome._id.toString()];       
                });
                self.db.collection('computations').find({'executable':executableName}).toArray(function(err, computations) {
                    if (err)
                        console.log(err);
                    else if (computations) {
                        var executables = {}, computationsProcessed = 0, hits = [];
                        for (var j=0 ; j < computations.length ; j++) {
                            var computation = computations[j],
                                inputs = computation['inputs'];
                            for (var i=0 ; i < inputs.length ; i++) {
                                if (inputs[i].match(/@genomes/) && organismNames2GenomesIds[organismName].indexOf(inputs[i].split('@genomes')[0]) != -1) {
                                    hits.push(computation);
                                    break;
                                }
                            }
                            computationsProcessed++;
                            if ( computationsProcessed == computations.length) {
                                self.emit("got computations", hits);
                                break;
                            }
                        };
                    }
                });   
            }
        });       
    };

    this.getncRNAs = function(computationId) {
         self.db.collection('computations').findOne({_id:computationId}, function(err, computation) {
            if (err)
                console.log(err);
            else if (computation) {
                var outputs = computation['outputs'], ncRNAs = [], ncRNAs_ids = [], ncRNAsProcessed = 0;
                for (var i=0 ; i < outputs.length ; i++) {
                    if (outputs[i].match(/@ncRNAs/))
                        ncRNAs_ids.push(outputs[i].split('@ncRNAs')[0])
                }
                ncRNAs_ids.forEach(function(id) {
                    self.db.collection('ncRNAs').findOne({_id:id}, function(err, ncRNA) {
                        ncRNAsProcessed++;
                        if (err)
                            console.log(err)
                        else if (ncRNA) {
                            ncRNAs.push(ncRNA);
                            if (ncRNAsProcessed == ncRNAs_ids.length) {
                                self.emit("got ncRNAs per computation", ncRNAs);
                            }
                        }
                    });
                });  
            }
         });
    };

    this.removeExecutable = function (executableName, organismName) {
        self.db.collection('computations').find({"executable":executableName}).toArray(function(err, computations) {
            if (err)
                console.log(err);
            else if (computations) {
               computations.forEach(function(computation) {
                    var inputs = computation['inputs'];
                    for (var i=0 ; i < inputs.length ; i++) {
                        if (inputs[i].match(/@genomes/)) {
                            self.db.collection('genomes').findOne({_id:inputs[i].split('@genomes')[0]}, function(err, genome) {
                                if (err)
                                    console.log(err)
                                else if (genome && genome['organism'] == organismName) {
                                    computation.outputs.forEach(function(output_id) {
                                        var tokens = output_id.split('@');
                                        self.db.collection(tokens[1]).removeById(tokens[0], function(err, reply) {
                                            if (err)
                                                console.log(err);
                                            else if (reply)
                                                console.log(reply); 
                                        });
                                    });
                                    self.db.collection('computations').removeById(computation._id, function(err, reply) {
                                        if (err)
                                            console.log(err);
                                        else if (reply)
                                            console.log(reply); 
                                    });
                                }
                            });
                        }
                    }
               });
               self.emit("executable removed", executableName);  
            }      
        });
    };

    this.removeComputation = function (computationId) {
        self.db.collection('computations').findOne({_id:computationId}, function(err, computation) {
            if (err)
                console.log(err);
            else if (computation) {
                computation.outputs.forEach(function(output_id) {
                    var tokens = output_id.split('@');
                    self.db.collection(tokens[1]).removeById(tokens[0], function(err, reply) {
                        if (err)
                            console.log(err);
                        else if (reply)
                            console.log(reply); 
                    });
                });
                self.db.collection('computations').removeById(computation._id, function(err, reply) {
                    if (err)
                        console.log(err);
                    else if (reply)
                        console.log(reply);
                    self.emit("computation removed", computation); 
                });
            }      
        });
    };

    this.removeNcRNA = function (ncRNAId) {
        self.db.collection('ncRNAs').removeById(ncRNAId, function(err, reply) {
            if (err)
                console.log(err);
            else if (reply)
                console.log(reply);
            else {
                //we need to update the array of outputs for the corresponding computation entry
                self.db.collection('computations').find().toArray(function(err, computations) {
                    if (err)
                        console.log(err);
                    else if (computations) {
                        for (var i = 0 ; i < computations.length ; i++) {
                            var computation = computations[i], found = false;
                            computation.outputs.forEach(function(output_id) {
                                if (output_id.split('@')[0] == ncRNAId && output_id.split('@')[1] == 'ncRNAs') {
                                    computation.outputs.splice(computation.outputs.indexOf(ncRNAId),1);
                                    console.log(computation._id);
                                    self.db.collection('computations').update({_id:computation._id}, {$set:{'outputs': computation.outputs}}, {safe:true}, function(error, record) {
                                        if (error)
                                            console.log(error);
                                    });
                                    found = true;   
                                }
                            })
                            if (found)
                                break;
                        };
                    }
                });
                self.emit("ncRNA removed", ncRNAId); 
            } 
        });
    }; 
};

db.RnajsDb.prototype.__proto__ = EventEmitter.prototype;

var args = process.argv.splice(2);
