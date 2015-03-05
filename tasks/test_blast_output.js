#!/usr/bin/env node

var rna = require('../lib/rna');

var fastaParser1 = new rna.parsers.FileParser(),
	fastaParser2 = new rna.parsers.FileParser(),
    blast = null;

 fastaParser2.on("end", function(queryMolecules) {
    blast.on("database formatted", function() {
            console.log("Search for "+queryMolecules[0].name+" in "+blast.targetMolecules[0].name);
            console.log("Can take a while...");
            blast.blastn(queryMolecules[0]);
        })
        .on("hits", function(hits) {
            hits.forEach(function(hit) {
                console.log("New hit!!");
                console.log(hit);
            });
        });
    blast.formatDb();        
});

fastaParser1.on("end", function(targetMolecules) {
	blast = new rna.computations.Blast(targetMolecules);
    fastaParser2.parseFasta("/Users/fjossinet/tmp/test/DEHA2D16852r.fasta", 'DNA'); //the query     
});

fastaParser1.parseFasta("/Users/fjossinet/tmp/test/AQF_AB_Kodamaea_ohmeri_AQF_AB.scaffolds.fa", 'DNA'); 
        