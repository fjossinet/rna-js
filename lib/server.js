#!/usr/bin/env node

//this is a minimal REST server giving access to computational tools and to a MongoDB running instance

var args = process.argv.splice(2),
	application_root = __dirname,
	rnajs = require('./rna.js');

if (args.indexOf("help") != -1) {
	console.log("Usage: ./server.js [-conf configuration_file]\n");
	console.log('This script starts the rna-js server.\n');
	return;
}

var rnajsServer = undefined;	

if (args.indexOf("-conf") != -1) {
	configFile = path.normalize(path.resolve(args[args.indexOf("-conf")+1]));
	rnajsServer = new rnajs.network.RnajsServer(application_root, require(configFile.split('.js')[0]).config);
}
else { //we use the configuration file provided with this module
	configFile = path.normalize(path.join(application_root, '/default-config.js'));
	rnajsServer = new rnajs.network.RnajsServer(application_root, require(configFile.split('.js')[0]).config);
}

rnajsServer.start();



