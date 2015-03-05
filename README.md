rna-js -- Web Services to compute and query RNA data
====================================================

The module rna-js allows you to deploy Web Services:

* to compute RNA data using dozens of algorithms: 2D prediction from single molecule (RNAfold, Contrafold,...), 2D prediction/structural alignment from a set of aligned/unaligned molecules (mlocarna, FoldAlign,...), annotation of tertiary structures into extended secondary structures (RNAVIEW, MCAnnotate,...),... All the computational results are transformed to JSON. 
* to query a Data Warehouse storing ncRNA data and annotations extracted from public databases (PDB, Rfam, NCBI, Ensembl,...) and transformed to JSON.

#Step 1: installation of rna-js

To be able to run rna-js, you will need to install:

* [node.js](http://nodejs.org),
* [MongoDB](http://www.mongodb.org/).

Then create a new directory and type:

	$ npm install rna-js

#Step 2: installation of third-parties algorithms

This step is done automatically using the script 'install_algorithms.sh':
	
	$ ./node_modules/rna-js/files/scripts/install_algorithms.sh <absolute_path_for_your_directory>

Once done, add all the binaries in your PATH by adding the following line in your file .bashrc (or any configuration file for your shell): 

	$ source <your_directory>/setmyenv

# Step 3: feeding the Data Warehouse

rna-js allows you to:

* import data from public databases or from GFF3 files using the script 'ark.js'. For more details, type:

		$ ./node_modules/rna-js/lib/ark.js help

* to compute and store your own ncRNAs annotations using the script 'annotate.js'. For more details, type:

		$ ./node_modules/rna-js/lib/annotate.js help

## the script 'ark.js'

This script imports data from public databases or from GFF3 files and store them in a Data Warehouse named ARK (for "All RNA Knowledge").

To import data from public databases, type:

	$ ./node_modules/rna-js/lib/ark.js -conf configuration_file

The configuration file stores all the details that will be needed during the import. Copy the ./node_modules/rna-js/lib/config.js file in a safe place and edit it to fit your own purposes.

To import data stored in GFF3 files:

	$ ./node_modules/rna-js/lib/ark.js -gff3 file_path

##the script 'annotate.js':

This script allows you to compute and to import your own genome-wide annotation of ncRNAs into the Data Warehouse. Some annotation steps will need the availability of an ARK database in the running MongoDB instance.

The first step is to import and annotate genomic sequences stored in a FASTA file by typing:

	$ ./node_modules/rna-js/lib/annotate.js -conf configuration_file -org "My organism name" -seqs fasta_file -db mongod_database_name

Once the genomic sequences stored in the database, further annotations can be done by typing:

	$ ./node_modules/rna-js/lib/annotate.js -conf configuration_file -org "My organism name" -db mongod_database_name

Annotations stored in gff3 files can also be imported using:

	$ ./node_modules/rna-js/lib/annotate.js -conf configuration_file -org "My organism name" -gff3 gff3_file -db mongod_database_name

The configuration file stores all the details that will be needed during the annotation process. Copy the ./lib/default-config.js file in a safe place and edit it to fit your own purposes. If no configuration file is precised in the command line, the script uses the file ./lib/default-config.

In the configuration file, you will find a glite section allowing to deploy the annotation processes on a grid using the glite middleware. 

# Final Step: launching the server

The script 'server.js' launches the rna-js server. It will deploy REST web services to use the algorithms in a transparent way and to query the Data Warehouse.

The server can be launched by typing:

	$ ./node_modules/rna-js/lib/server.js -conf configuration_file

The configuration file stores all the details needed by the server process. Copy the ./lib/default-config.js file in a safe place and edit it to fit your own purposes. If no configuration file is precised in the command line, the script uses the file ./lib/default-config. Two sections of this file are important:

* the config.webserver section: define the host and port for the running web server,
* the config.mongodb section: define the host and port for the running MongoDb instance.

Once the server launched, open a browser at http://{config.webserver.host}:{config.webserver.port}.

#Related projects

* [rna-js clients](https://bitbucket.org/fjossinet/rna-js-clients): documentation and code samples for basic rna-js clients,
* [PyRNA](https://bitbucket.org/fjossinet/pyrna) a rna-js client written with Python. It transforms the JSON data recovered from the rna-js Web Services into [pandas data structures](http://pandas.pydata.org/) more adapted to data analysis and modeling,
* [Assemble2](http://www.bioinformatics.org/assemble/) a Java GUI using rna-js Web Services to construct RNA 3D models.


