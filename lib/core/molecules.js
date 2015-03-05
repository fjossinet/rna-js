var molecules = exports,
    EventEmitter = require('events').EventEmitter,
    utils = require('./utils'),
    ObjectID = require('mongodb').ObjectID;

molecules.DNA = function(name, sequence) {
    var self = this;
    this.name = name;
    this.sequence = "";
    this._id = new ObjectID();

    this.toFasta = function (single_line) {
        var toFasta = ">"+self.name+"\n";
        if (single_line) {
            toFasta = toFasta.concat(self.sequence+"\n");
        } else {
            var chunks = utils.chunks(self.sequence.split(''),79).map(function(chunk) {
                return chunk.join("");
            });
            toFasta = toFasta.concat(chunks.join('\n'));
        }
        return toFasta;
    };

    this.addResidue = function(residue) {
        self.sequence += residue; 
    };

    this.complement = function() {
        var basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'},
            letters = self.sequence.split(''),
            complement = "";
        letters.forEach(function(letter) {
            complement +=  basecomplement[letter] || '?'; 
        });
        return complement;
    }

    this.reversecomplement = function() { 
        var basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'},
            letters = self.sequence.split('').reverse(),
            complement = "";
        letters.forEach(function(letter) {
            complement +=  basecomplement[letter] || '?'; 
        });
        return complement;
    }

    if (sequence) {
        sequence.toUpperCase().split('').forEach(function(residue) {
            self.addResidue(residue);
        });
    }
};

molecules.RNA = function(name, sequence) {
    var self = this;
    this.name = name || "RNA";
    this.sequence = "";
    this.source = "?:?:?";
    this._id = new ObjectID();
    this.family = undefined;
    this.organism = undefined;
    this.modified_ribonucleotides = [];

    this.addResidue = function(residue) {
        if (regular_nucleotides[residue]) {
            residue = regular_nucleotides[residue]; 
        } else if (modified_ribonucleotides[residue]) {
            self.modified_ribonucleotides.push([residue,self.sequence.length+1]);   
            residue = modified_ribonucleotides[residue];
        } else if (regular_aminoacids[residue]) {
            return;
        }
        switch (residue) {
            case 'A' : self.sequence += residue; break;
            case 'U' : self.sequence += residue; break;
            case 'G' : self.sequence += residue; break;
            case 'C' : self.sequence += residue; break;
            case '\.': 
            case '-' :
            case '_' : self.sequence += '-'; break; //an aligned RNA
            default: console.log("Unknown residue "+residue) ; self.sequence += residue;
        }
    };
    
    this.getComplement = function() {
        var complementSeq = self.sequence.split('').reverse().map(function(residue){
            switch (residue) {
                case 'A' : return 'U'; break;
                case 'U' : return 'A'; break;
                case 'G' : return 'C'; break;
                case 'C' : return 'G'; break;
                case '-' : return '-'; break; //an aligned RNA
                default: return '?'; //for all the other ones (IUPAC symbols,...), no real rule
            }
        }).join('');
        return new RNA(self.name+" complement", complementSeq);
    }
    
    this.getGcContent = function() {
        var gcCount = self.sequence.split('').filter(function(residue) {
            return residue == 'G' || residue == 'C';
        }).length;
        return gcCount/self.sequence.length;
    }
    
    this.toFasta = function (single_line) {
        var toFasta = ">"+self.name+"\n";
        if (single_line) {
            toFasta = toFasta.concat(self.sequence+"\n");
        } else {
            var chunks = utils.chunks(self.sequence.split(''),79).map(function(chunk) {
                return chunk.join("");
            });
            toFasta = toFasta.concat(chunks.join('\n'));
        }
        return toFasta;
    }

    if (sequence) {
        sequence.toUpperCase().split('').forEach(function(residue) {
            self.addResidue(residue);
        });
    }
};

molecules.RNA.prototype.__proto__ = EventEmitter.prototype;

molecules.makeRandomRNA = function(size, name) {
    var sequence = "",
        residues = ['A','U','G','C'];
        
    for( var i = 0; i < size; i++ )
        sequence += residues[Math.floor(Math.random() * residues.length)];

    return new RNA(name || "RNA", sequence);
};


regular_nucleotides = {
    "ADE": "A",
    "URA": "U",
    "URI": "U",
    "GUA": "G",
    "CYT": "C",
};

regular_aminoacids = {
    "CYS": "C",
    "HIS": "H",
    "ILE": "I",
    "MET": "M",
    "SER": "S",
    "VAL": "V",
    "ALA": "A",
    "GLY": "G",
    "LEU": "L",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ARG": "R",
    "TYR": "Y",
    "TRP": "W",
    "ASP": "D",
    "ASN": "N",
    "GLU": "E",
    "GLN": "Q",
    "LYS": "K",
};

modified_ribonucleotides = {
    "T": "U",
    "PSU": "U",
    "I": "A",
    "N": "U",
    "S": "U",
    "+A": "A",
    "+C": "C",
    "+G": "G",
    "+I": "I",
    "+T": "U",
    "+U": "U",
    "PU": "A",
    "YG": "G",
    "1AP": "G",
    "1MA": "A",
    "1MG": "G",
    "2DA": "A",
    "2DT": "U",
    "2MA": "A",
    "2MG": "G",
    "4SC": "C",
    "4SU": "U",
    "5IU": "U",
    "5MC": "C",
    "5MU": "U",
    "5NC": "C",
    "6MP": "A",
    "7MG": "G",
    "A23": "A",
    "AD2": "A",
    "AET": "A",
    "AMD": "A",
    "AMP": "A",
    "APN": "A",
    "ATP": "A",
    "AZT": "U",
    "CCC": "C",
    "CMP": "A",
    "CPN": "C",
    "DAD": "A",
    "DCT": "C",
    "DDG": "G",
    "DG3": "G",
    "DHU": "U",
    "DOC": "C",
    "EDA": "A",
    "G7M": "G",
    "GDP": "G",
    "GNP": "G",
    "GPN": "G",
    "GTP": "G",
    "GUN": "G",
    "H2U": "U",
    "HPA": "A",
    "IPN": "U",
    "M2G": "G",
    "MGT": "G",
    "MIA": "A",
    "OMC": "C",
    "OMG": "G",
    "OMU": "U",
    "ONE": "U",
    "P2U": "P",
    "PGP": "G",
    "PPU": "A",
    "PRN": "A",
    "PST": "U",
    "QSI": "A",
    "QUO": "G",
    "RIA": "A",
    "SAH": "A",
    "SAM": "A",
    "T23": "U",
    "T6A": "A",
    "TAF": "U",
    "TLC": "U",
    "TPN": "U",
    "TSP": "U",
    "TTP": "U",
    "UCP": "U",
    "VAA": "A",
    "YYG": "G",
    "70U": "U",
    "12A": "A",
    "2MU": "U",
    "127": "U",
    "125": "U",
    "126": "U",
    "MEP": "U",
    "TLN": "U",
    "ADP": "A",
    "TTE": "U",
    "PYO": "U",
    "SUR": "U",
    "PSD": "A",
    "S4U": "U",
    "CP1": "C",
    "TP1": "U",
    "NEA": "A",
    "GCK": "C",
    "CH": "C",
    "EDC": "G",
    "DFC": "C",
    "DFG": "G",
    "DRT": "U",
    "2AR": "A",
    "8OG": "G",
    "IG": "G",
    "IC": "C",
    "IGU": "G",
    "IMC": "C",
    "GAO": "G",
    "UAR": "U",
    "CAR": "C",
    "PPZ": "A",
    "M1G": "G",
    "ABR": "A",
    "ABS": "A",
    "S6G": "G",
    "HEU": "U",
    "P": "G",
    "DNR": "C",
    "MCY": "C",
    "TCP": "U",
    "LGP": "G",
    "GSR": "G",
    "X": "G",
    "R": "A",
    "Y": "A",
    "E": "A",
    "GSS": "G",
    "THX": "U",
    "6CT": "U",
    "TEP": "G",
    "GN7": "G",
    "FAG": "G",
    "PDU": "U",
    "MA6": "A",
    "UMP": "U",
    "SC": "C",
    "GS": "G",
    "TS": "U",
    "AS": "A",
    "ATD": "U",
    "T3P": "U",
    "5AT": "U",
    "MMT": "U",
    "SRA": "A",
    "6HG": "G",
    "6HC": "C",
    "6HT": "U",
    "6HA": "A",
    "55C": "C",
    "U8U": "U",
    "BRO": "U",
    "BRU": "U",
    "5IT": "U",
    "ADI": "A",
    "5CM": "C",
    "IMP": "G",
    "THM": "U",
    "URI": "U",
    "AMO": "A",
    "FHU": "P",
    "TSB": "A",
    "CMR": "C",
    "RMP": "A",
    "SMP": "A",
    "5HT": "U",
    "RT": "U",
    "MAD": "A",
    "OXG": "G",
    "UDP": "U",
    "6MA": "A",
    "5IC": "C",
    "SPT": "U",
    "TGP": "G",
    "BLS": "A",
    "64T": "U",
    "CB2": "C",
    "DCP": "C",
    "ANG": "G",
    "BRG": "G",
    "Z": "A",
    "AVC": "A",
    "5CG": "G",
    "UDP": "U",
    "UMS": "U",
    "BGM": "G",
    "SMT": "U",
    "DU": "U",
    "CH1": "C",
    "GH3": "G",
    "GNG": "G",
    "TFT": "U",
    "U3H": "U",
    "MRG": "G",
    "ATM": "U",
    "GOM": "A",
    "UBB": "U",
    "A66": "A",
    "T66": "U",
    "C66": "C",
    "3ME": "A",
    "A3P": "A",
    "ANP": "A",
    "FA2": "A",
    "9DG": "G",
    "GMU": "U",
    "UTP": "U",
    "5BU": "U",
    "APC": "A",
    "DI": "I",
    "UR3": "U",
    "3DA": "A",
    "DDY": "C",
    "TTD": "U",
    "TFO": "U",
    "TNV": "U",
    "MTU": "U",
    "6OG": "G",
    "E1X": "A",
    "FOX": "A",
    "CTP": "C",
    "D3T": "U",
    "TPC": "C",
    "7DA": "A",
    "7GU": "U",
    "2PR": "A",
    "CBR": "C",
    "I5C": "C",
    "5FC": "C",
    "GMS": "G",
    "2BT": "U",
    "8FG": "G",
    "MNU": "U",
    "AGS": "A",
    "NMT": "U",
    "NMS": "U",
    "UPG": "U",
    "G2P": "G",
    "2NT": "U",
    "EIT": "U",
    "TFE": "U",
    "P2T": "U",
    "2AT": "U",
    "2GT": "U",
    "2OT": "U",
    "BOE": "U",
    "SFG": "G",
    "CSL": "I",
    "PPW": "G",
    "IU": "U",
    "D5M": "A",
    "ZDU": "U",
    "DGT": "U",
    "UD5": "U",
    "S4C": "C",
    "DTP": "A",
    "5AA": "A",
    "2OP": "A",
    "PO2": "A",
    "DC": "C",
    "DA": "A",
    "LOF": "A",
    "ACA": "A",
    "BTN": "A",
    "PAE": "A",
    "SPS": "A",
    "TSE": "A",
    "A2M": "A",
    "NCO": "A",
    "A5M": "C",
    "M5M": "C",
    "S2M": "U",
    "MSP": "A",
    "P1P": "A",
    "N6G": "G",
    "MA7": "A",
    "FE2": "G",
    "AKG": "G",
    "SIN": "G",
    "PR5": "G",
    "GOL": "G",
    "XCY": "G",
    "5HU": "U",
    "CME": "C",
    "EGL": "G",
    "LC": "C",
    "LHU": "U",
    "LG": "G",
    "PUY": "U",
    "PO4": "U",
    "PQ1": "U",
    "ROB": "U",
    "O2C": "C",
    "C30": "C",
    "C31": "C",
    "C32": "C",
    "C33": "C",
    "C34": "C",
    "C35": "C",
    "C36": "C",
    "C37": "C",
    "C38": "C",
    "C39": "C",
    "C40": "C",
    "C41": "C",
    "C42": "C",
    "C43": "C",
    "C44": "C",
    "C45": "C",
    "C46": "C",
    "C47": "C",
    "C48": "C",
    "C49": "C",
    "C50": "C",
    "A30": "A",
    "A31": "A",
    "A32": "A",
    "A33": "A",
    "A34": "A",
    "A35": "A",
    "A36": "A",
    "A37": "A",
    "A38": "A",
    "A39": "A",
    "A40": "A",
    "A41": "A",
    "A42": "A",
    "A43": "A",
    "A44": "A",
    "A45": "A",
    "A46": "A",
    "A47": "A",
    "A48": "A",
    "A49": "A",
    "A50": "A",
    "G30": "G",
    "G31": "G",
    "G32": "G",
    "G33": "G",
    "G34": "G",
    "G35": "G",
    "G36": "G",
    "G37": "G",
    "G38": "G",
    "G39": "G",
    "G40": "G",
    "G41": "G",
    "G42": "G",
    "G43": "G",
    "G44": "G",
    "G45": "G",
    "G46": "G",
    "G47": "G",
    "G48": "G",
    "G49": "G",
    "G50": "G",
    "T30": "U",
    "T31": "U",
    "T32": "U",
    "T33": "U",
    "T34": "U",
    "T35": "U",
    "T36": "U",
    "T37": "U",
    "T38": "U",
    "T39": "U",
    "T40": "U",
    "T41": "U",
    "T42": "U",
    "T43": "U",
    "T44": "U",
    "T45": "U",
    "T46": "U",
    "T47": "U",
    "T48": "U",
    "T49": "U",
    "T50": "U",
    "U30": "U",
    "U31": "U",
    "U32": "U",
    "U33": "U",
    "U34": "U",
    "U35": "U",
    "U36": "U",
    "U37": "U",
    "U38": "U",
    "U39": "U",
    "U40": "U",
    "U41": "U",
    "U42": "U",
    "U43": "U",
    "U44": "U",
    "U45": "U",
    "U46": "U",
    "U47": "U",
    "U48": "U",
    "U49": "U",
    "U50": "U",
    "UFP": "U",
    "UFR": "U",
    "UCL": "U",
    "3DR": "U",
    "CBV": "C",
    "HFA": "A",
    "MMA": "A",
    "DCZ": "C",
    "GNE": "C",
    "A1P": "A",
    "6IA": "A",
    "CTG": "G",
    "5FU": "U",
    "2AD": "A",
    "T2T": "U",
    "XUG": "G",
    "2ST": "U",
    "5PY": "U",
    "4PC": "C",
    "US1": "U",
    "M5C": "C",
    "DG": "G",
    "DA": "A",
    "DT": "U",
    "DC": "C",
    "P5P": "A",
    "FMU": "U"
};