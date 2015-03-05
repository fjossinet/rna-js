var utils = exports;

utils.range = function (begin, end) {
    var range = [];
    for (var i = begin; i < end; i++) {
        range.push(i);
    }
    return range;
};

utils.bnToBps = function(bn) {
    var bps = [],
        chars = bn.split(''),
        leftPositions = [];
    for (var i = 0; i < chars.length; i++) {
        if (chars[i] == '(') {
            leftPositions.push(i+1);
        } else if (chars[i] == ')') {
            bps.push([leftPositions.pop(),i+1]);
        }
    }
    return bps;
};

utils.chunks = function(array,size) {
    var chunks = [];
    for (var i = 0 ; i < array.length ; i+=size) {
        chunks.push(array.slice(i,i+size));
    }
    return chunks;
};

utils.Location = function(ends) {
    var self = this;
    this.ends = utils.chunks(ends,2);

    this.getLength = function() {
        var length = 0;
        self.ends.forEach(function(chunk) {
            length += chunk[1] - chunk[0] + 1;
        });
        return length;    
    }

    this.getStart = function() {
        return self.ends[0][0];
    };

    this.getEnd = function() {
        return self.ends[self.ends.length-1][1];
    };
    
    this.hasPosition = function(position) {
        return self.ends.some(function(chunk) {
            return position >= chunk[0] && position <= chunk[1];
        });
    };
    
    this.getPositions = function() {
        var positions = [];
        self.ends.forEach(function(chunk) {
            for (var i = chunk[0]; i <= chunk[1]; i++ ) {
                positions.push(i);
            }
        });
        return positions;
    }
    
    this.removeLocation = function(location) {
        location.getPositions().forEach(function(pos) {
            self.removePosition(pos)
        });
    }
    
    this.removePosition = function(position) {
        for (var i = 0; i < self.ends.length; i++) {
            var chunk = self.ends[i];
            if (position >= chunk[0] && position <= chunk[1]) {
                if (chunk[0] == chunk[1]) {
                    self.ends.splice(i,1);
                } else if (position == chunk[0]) {
                    self.ends.splice(i,1, [position+1,chunk[1]]);
                } else if (position == chunk[1]) {
                    self.ends.splice(i,1, [chunk[0],chunk[1]-1]);
                } else {
                    self.ends.splice(i,1, [chunk[0],position-1], [position+1,chunk[1]]);
                }
                break;
            }
        }
    };
};