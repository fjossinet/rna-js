//########################################################################################
// All the directories precised in this file have to exist before to use rna-js
//########################################################################################

var config = {};

//the webserver configuration (launched by typing ./lib/server.js)
config.webserver = {};
config.webserver.host = "localhost";
config.webserver.port = 8000;

//configuration for the external data
config.data = {};
config.data.root = "/tmp/"

//the mongodb instance to use 
config.mongodb = {};
config.mongodb.host = 'localhost';
config.mongodb.port = 27017;
config.mongodb.dbpath = "/Users/fjossinet/mongo-db"

//the name and details of the ARK database used by the script ./lib/ark.js
config.mongodb.ark = {};
config.mongodb.ark.current = 'ark_1_0';

//uncomment the lines for the public databases you want to import
//config.mongodb.ark.genolevures = {};
//config.mongodb.ark.rfam = {};
//config.mongodb.ark.pdb = {'assemble2':false, 'non_redundant_set':false};
//config.mongodb.ark.ncbi = {};
//config.mongodb.ark.ensembl = {};

exports.config = config;
