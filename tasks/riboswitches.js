#!/usr/bin/env node 

var args = process.argv.splice(2),
    rna = require('../lib/rna');

var rfam =  new rna.db.Rfam();

rfam.on('entry parsed', function(molecules, consensus2D) {
	consensus2D.forEach(function(basebaseInteraction) {

		var currentResidue = undefined,
			currentPairedResidue = undefined,
			previousResidue = undefined,
			previousPairedResidue = undefined,
			covariationScore = 0;

		molecules.forEach(function(molecule) {
			console.log(molecule);
			currentResidue = molecule.sequence[basebaseInteraction.location.getStart()-1];
			currentPairedResidue = molecule.sequence[basebaseInteraction.location.getEnd()-1];
			if (currentResidue == 'A' && currentPairedResidue == 'U' ||
				currentResidue == 'U' && currentPairedResidue == 'A' ||
				currentResidue == 'G' && currentPairedResidue == 'C' ||
				currentResidue == 'C' && currentPairedResidue == 'G' ||
				currentResidue == 'G' && currentPairedResidue == 'U' ||
				currentResidue == 'U' && currentPairedResidue == 'G' ) {

				if (!previousResidue && !previousPairedResidue) {//start of counting
					previousResidue = currentResidue;
					previousPairedResidue = currentPairedResidue;
				} else if (previousResidue != currentResidue || previousPairedResidue != currentPairedResidue) {
					covariationScore++;
					previousResidue = currentResidue;
					previousPairedResidue = currentPairedResidue;
				}	
			}
		});

		var percent= (covariationScore/molecules.length)*100;

		console.log(percent.toFixed(2));

		//console.log(JSON.stringify(basebaseInteraction));

		covariationScore = 0;

	});
});

rfam.getEntry("RF00059", 'seed');