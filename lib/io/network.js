var network = exports,
    http = require('http'),
    express = require("express"),
    path = require("path"),
    configFile = undefined,
    config = undefined,
    mongo = require("mongodb"),
    rnajs = require('../rna.js'),
    EventEmitter = require('events').EventEmitter,
    ObjectID = require('mongodb').ObjectID;

//A minimal REST server giving access to computational tools and to a MongoDB running instance
network.RnajsServer = function(rootDir, config) {

    var self = this;
    this.config = config;
    this.ark = config.mongodb ? new rnajs.db.RnajsDb(config.mongodb.ark.current, config.mongodb.host, config.mongodb.port) : undefined;
    this.mongodbpath = config.mongodb && config.mongodb.dbpath ? path.normalize(config.mongodb.dbpath) : undefined;
    this.app = express(config.webserver.host);

    this.app.configure(function () {
      self.app.use(express.limit('10mb')); //here we change the upload limit size
      self.app.use(express.bodyParser());
      self.app.use(express.methodOverride());
      self.app.use(self.app.router);
      self.app.use(express.static(path.normalize(path.join(rootDir, "public"))));
      self.app.use(express.errorHandler({ dumpExceptions: true, showStack: true }));
    });

    //################### GET WEBSERVICES ##############################//

    // MONGODB QUERIES

    this.app.get('/api/db/:db?/:collection?/:id?', function(req, res) { //adapted from https://github.com/tdegrunt/mongodb-rest/blob/master/lib/rest.js

        var query = {}, options = {};

        if (req.params.id) 
            query = {'_id': req.params.id};

        if (!req.params.db) { //if no db precised, we send the list of annotations projects //USED BY BACARA 
            res.send( Object.keys(config.mongodb.genomic_annotation_projects));    
        }
        else {
      
            var db = new mongo.Db(req.params.db, new mongo.Server(config.mongodb.host, config.mongodb.port, {'auto_reconnect':true}), {'safe':true});
          
            db.open(function(err,db) {
                db.collection(req.params.collection, function(err, collection) {
                    collection.find(query, options, function(err, cursor) {
                        cursor.toArray(function(err, docs){
                            if (req.params.id && req.params.collection == 'ncRNAs' && docs.length == 1) { //USED BY BACARA 
                                ncRNA = docs.pop()
                                if (ncRNA['ortholog'].split(':')[1] ==  'rfam') {
                                    tmp.file(function _tempFileCreated(err, file_path, fd) {
                                        var outStream = fs.createWriteStream(file_path, {flags:'w'});
                                        if (ncRNA['genomicStrand'] == '+')
                                            outStream.write(">"+req.params.id+"/"+ncRNA['genomicPositions'][0]+"-"+ncRNA['genomicPositions'][1]+"\n"+ncRNA['sequence']);
                                        else
                                            outStream.write(">"+req.params.id+"/"+ncRNA['genomicPositions'][1]+"-"+ncRNA['genomicPositions'][0]+"\n"+ncRNA['sequence']);     
                                        exec.exec("cmalign.py -f " + file_path +" -id " + ncRNA['ortholog'].split(':')[2], function(err, stderr, stdout) {
                                            if (err)
                                                console.log(err);
                                            res.send(stderr);
                                        });
                                    });    
                                }

                            }
                            else if (!req.params.id && req.params.collection == 'genomes' ) { //USED BY BACARA 
                                genomes = []
                                for (var i = 0 ; i < docs.length ; i++) //Genomes can be huged, so in this case, only some informations are sent
                                    genomes.push({
                                        '_id': docs[i]._id,
                                        'name': docs[i].name 
                                    })
                                res.set('Content-Type', 'text/plain');
                                res.send(JSON.stringify(genomes));   
                            }
                            else {
                                res.header('Content-Type', 'application/json');
                                res.send(docs);
                            }
                            db.close();
                        });
                    });
                });
            });
        }         

    });

    /*this.app.get('/api/ark/:output/:source/:feature/:organism/:class/:id', function(req, res) {
        var query = {}, 
            options = {},
            output = req.params.output,
            db_name = req.params.source,
            collection_name = req.params.feature;

        if (req.params.id != '*')
            query = {'_id': req.params.id};
        else {
            if (req.params.organism != '*')
                query['organism'] = req.params.organism;
            if (req.params.class != '*')
                query['class'] = req.params.class;   
        }    

        var db = new mongo.Db(db_name, new mongo.Server(config.mongodb.host, config.mongodb.port, {'auto_reconnect':true}), {'safe':true});

        db.open(function(err,db) {
            db.collection(collection_name, function(err, collection) {
                collection.find(query, options, function(err, cursor) {
                    cursor.toArray(function(err, docs){
                        if (output = 'json') {
                            res.header('Content-Type', 'application/json');
                            res.send(docs);
                        }
                        db.close();
                    });
                });
            });
        });

    });*/

    /*this.app.get('/api/db/:db/:query', function (req, res) {
        if (req.query.coll && req.query.query)
            self.ark.db.collection(req.query.coll).find(JSON.parse(req.query.query)).toArray(function(err, ncRNAs) {
                if (err)
                    console.log(err);
                res.send(ncRNAs);
            });
        else if (req.query.coll) {//the client get all the data from a collection
            self.ark.db.collection(req.query.coll).find().toArray(function(err, entries) {
                if (err)
                    console.log(err);
                if (req.query.coll == 'ncRNAs' && req.query.output && req.query.output == "gff3") {
                    var output = "";
                    self.ark.db.collection('genomes').find().toArray(function(err, genomes) {
                        if (err)
                            console.log(err);
                        var genomes_idsToName = {}
                        genomes.forEach(function(genome) {
                            genomes_idsToName[""+genome._id] = genome.name;
                        });
                        entries.forEach(function(ncRNA) {
                            output += genomes_idsToName[ncRNA.genome.split('@genomes')[0]]+"\t"+ncRNA.tool.split(' ')[0]+"\tncRNA\t"+ncRNA.genomicPositions[0]+"\t"+ncRNA.genomicPositions[1]+"\t"+ncRNA.score+"\t"+ncRNA.genomicStrand+"\t.\tID="+ncRNA._id+"\n";
                        });
                        res.set('Content-Type', 'text/plain');
                        res.send(output);
                    });
                } else {
                    res.set('Content-Type', 'text/plain');
                    res.send(entries);
                }
            });
        } else
            res.redirect("https://bitbucket.org/fjossinet/rna-js-clients/raw/default/README.md");
    });*/

    // COMPUTATIONS

    this.app.get('/api/compute/2d', function (req, res) {
        if (req.query.seq && req.query.tool && (req.query.tool == 'contrafold' || req.query.tool == 'rnafold')) {
            var parser = new rnajs.parsers.StringParser();
            parser.on('end', function(molecules) {
                if (molecules.length > 0) {
                    if (req.query.tool == 'contrafold') {
                        var contrafold = new rnajs.computations.Contrafold();
                        contrafold.on('end', function(secondaryStructure) {
                            res.send(secondaryStructure);
                        });
                        contrafold.fold(molecules[0], req.query.output);    
                    } else if (req.query.tool == 'rnafold') {
                        var rnafold = new rnajs.computations.Rnafold();
                        rnafold.on('end', function(secondaryStructure) {
                            res.send(secondaryStructure);
                        });
                        rnafold.fold(molecules[0], req.query.output);    
                    }        
                }
            });
            parser.parseFasta('>'+(req.query.name || 'rna')+"\n"+req.query.seq);
        } else if (req.query.pdbid && req.query.tool == 'rnaview') { 
            var pdb = new rnajs.db.PDB(),
                rnaview = new rnajs.computations.Rnaview(),
                result= [],
                total = 0,
                parsed = 0;
            
            pdb.on('entry parsed', function(tertiaryStructures) {
                total = tertiaryStructures.length;
                rnaview.on('end', function(annotated2D) {
                    if (req.query.output != 'rnaml')
                        annotated2D['2D'].findJunctions();
                    result.push(annotated2D);
                    parsed++;
                    if (parsed == total) {
                        res.set('Content-Type', 'text/plain');
                        res.send(JSON.stringify(result));
                    } else
                        rnaview.annotate(tertiaryStructures[parsed], req.query.output == 'rnaml');

                }).on('error', function(error) {
                    console.log(error);
                    result.push({"name":"2D","rna":tertiaryStructures[parsed].rna});
                    parsed++;                                   
                    if (parsed == total) {
                        res.set('Content-Type', 'text/plain');
                        res.send(JSON.stringify(result));
                    } else 
                        rnaview.annotate(tertiaryStructures[parsed], req.query.output == 'rnaml');
                });
                rnaview.annotate(tertiaryStructures[parsed], req.query.output == 'rnaml');
            });
            pdb.getEntry(req.query.pdbid);
        } else
            res.redirect("https://bitbucket.org/fjossinet/rna-js-clients/raw/default/README.md");
    });

    //################### POST WEBSERVICES ##############################//

    // MONGODB QUERIES

    this.app.post('/api/ark', function (req, res) { //USED BY ASSEMBLE2 TO GET THE JUNCTIONS
        if (req.body.coll && req.body.query) {
            self.ark.db.collection(req.body.coll).find(JSON.parse(req.body.query)).toArray(function(err, ncRNAs) {
                if (err)
                    console.log(err);
                res.send(ncRNAs);
            });
        } else if (req.body.coll && req.body.id) {
            self.ark.db.collection(req.body.coll).findOne({_id:new ObjectID(req.body.id)},function(err, result) { //in the first version of ARK, the id is stored as an ObjectId
                if (err)
                    console.log(err);
                if (result)
                    res.send(result);
                else {
                    self.ark.db.collection(req.body.coll).findOne({_id:req.body.id},function(err, result) { //in the first version of ARK, the id is stored as a String
                        if (err)
                            console.log(err);
                        if (result)
                            res.send(result);
                    });    
                }                   
                    
            });
        } else if (req.body.coll) {//the client get all the data from a collection
            self.ark.db.collection(req.body.coll).find().toArray(function(err, entries) {
                if (err)
                    console.log(err);
                if (req.body.coll == 'ncRNAs' && req.body.output && req.body.output == "gff3") {
                    var output = "";
                    userDb.db.collection('genomes').find().toArray(function(err, genomes) {
                        if (err)
                            console.log(err);
                        var genomes_idsToName = {}
                        genomes.forEach(function(genome) {
                            genomes_idsToName[""+genome._id] = genome.name;
                        });
                        entries.forEach(function(ncRNA) {
                            output += genomes_idsToName[ncRNA.genome.split('@genomes')[0]]+"\t"+ncRNA.tool.split(' ')[0]+"\tncRNA\t"+ncRNA.genomicPositions[0]+"\t"+ncRNA.genomicPositions[1]+"\t"+ncRNA.score+"\t"+ncRNA.genomicStrand+"\t.\tID="+ncRNA._id+"\n";
                        });
                        res.set('Content-Type', 'text/plain');
                        res.send(output);
                    });
                }
                else {
                    res.set('Content-Type', 'text/plain');
                    res.send(entries);
                }
            });
        } else
            res.redirect("https://bitbucket.org/fjossinet/rna-js-clients/raw/default/README.md");
    }); 

    // COMPUTATIONS
    this.app.post('/api/compute/2d', function (req, res) {
        if (req.body.data && req.body.tool) {
            var parser = new rnajs.parsers.StringParser();
            if (req.body.data.match(/^>/)) { 
                parser.on('end', function(molecules) {
                    if (molecules.length >= 2) { //structural alignment
                        if (req.body.tool == 'foldalign') {      
                        } else if (req.body.tool == 'mlocarna') {
                            var mlocarna = new rnajs.computations.Mlocarna();
                            mlocarna.on('end', function(alignedMolecules, consensus2D) {
                                res.send(JSON.stringify({'sequences':alignedMolecules, 'consensus2D':consensus2D}));    
                            });
                            mlocarna.align(molecules);     
                        }       
                    } else { //USED BY ASSEMBLE2 TO PREDICT A 2D
                        if (req.body.tool == 'contrafold') {
                            var contrafold = new rnajs.computations.Contrafold();
                            contrafold.on('end', function(secondaryStructure) {
                                res.send(secondaryStructure);
                            });
                            contrafold.fold(molecules[0], req.body.output);    
                        } else if (req.body.tool == 'rnafold') {
                            var rnafold = new rnajs.computations.Rnafold();
                            rnafold.on('end', function(secondaryStructure) {
                                res.send(secondaryStructure);
                            });
                            rnafold.fold(molecules[0], req.body.output);    
                        } 
                    }
                });
                parser.parseFasta(req.body.data);
            } else if (req.body.tool == 'rnaview') { //USED BY ASSEMBLE2 TO ANNOTATE A 3D
                var rnaview = new rnajs.computations.Rnaview(),
                    total = 0, 
                    parsed = 0,
                    result = [];
                parser.on('end', function(tertiaryStructures) {
                    total = tertiaryStructures.length;
                    rnaview.on('end', function(annotated2D) {
                        if (req.body.output != 'rnaml')
                            annotated2D['2D'].findJunctions();
                        result.push(annotated2D);
                        parsed++;
                        if (parsed == total) {
                            res.set('Content-Type', 'text/plain');
                            res.send(JSON.stringify(result));
                        } else
                            rnaview.annotate(tertiaryStructures[parsed], req.body.output == 'rnaml');

                    }).on('error', function(error) {
                        console.log(error);
                        result.push({"name":"2D","rna":tertiaryStructures[parsed].rna});
                        parsed++;                                   
                        if (parsed == total) {
                            res.set('Content-Type', 'text/plain');
                            res.send(JSON.stringify(result));
                        } else 
                            rnaview.annotate(tertiaryStructures[parsed], req.body.output == 'rnaml');
                    });
                    rnaview.annotate(tertiaryStructures[parsed], req.body.output == 'rnaml');
                });
                parser.parsePDB(req.body.data);    
            }  else
                res.redirect("https://bitbucket.org/fjossinet/rna-js-clients/raw/default/README.md");
        }    
    });

    this.app.post('/api/compute/2dplot', function (req, res) { //USED BY ASSEMBLE2 TO PLOT A 2D
        if (req.body.data) {
            tmp.file(function _tempFileCreated(err, file_path, fd) {
                var outStream = fs.createWriteStream(file_path, {flags:'w'});
                outStream.write(req.body.data);
                exec.exec("plot_2d.py -f " + file_path, function(err, stderr, stdout) {
                    if (err)
                        console.log(err);
                    res.send(stderr);
                });
            });
        }  else
            res.redirect("https://bitbucket.org/fjossinet/rna-js-clients/raw/default/README.md");  
    });

    this.start = function () {

        if (self.mongodbpath) {

            if (!path.existsSync(self.mongodbpath))
                fs.mkdirSync(self.mongodbpath);

            console.log("Launch MongoDB server with database located at "+self.mongodbpath);

            exec.exec("mongod -dbpath "+self.mongodbpath, function(err, stderr, stdout) {
                if (err)
                    console.log(stderr);
            });
        }

        server = http.createServer(self.app);
        server.listen(self.config.webserver.port);
        console.log('Rnajs server ready!! The Web services are accessible at http://'+self.config.webserver.host+":"+self.config.webserver.port+"/api/");
    };
};

network.RnajsServer.prototype.__proto__ = EventEmitter.prototype;