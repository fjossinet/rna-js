rna-js -- Produce and Analyse RNA Data with Node.js
===================================================

To use rna-js in your node.js project, do the following:

	var rna = require('rna-js/lib/rna.js');
	var parser = new rna.parsers.StringParser();
	parser.on('end', function(tertiaryStructures) {
		//...
	});
	parser.parsePDB(myPDBdata);

* ./core/molecules.js: RNA and DNA concepts,
* ./core/features.js: RNA 2D and 3D concepts,
* ./core/utils.js: some utils functions,
* ./io/computations.js: wrappers for algorithms (Blast, RNAfold,...),
* ./io/db.js: wrappers for databases (RFAM, NCBI,...),
* ./io/parsers.js: parsers for bioinformatics format (PDB, FASTA,...).




