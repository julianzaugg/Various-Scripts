#!/usr/bin/python

import argparse
import sys

from collections import Counter

class Alphabet(object):
    """ Defines an immutable biological alphabet (e.g. the alphabet for DNA is AGCT) 
    that can be used to create sequences (see sequence.py).
    We use alphabets to define "tuple" tables, where entries are keyed by combinations
    of symbols of an alphabet (see class TupleStore below). 
    Alphabets are used to define probability distributions for stochastic events
    (see prob.py). """
    
    def __init__(self, symbolString):
        """ Construct an alphabet from a string of symbols. Lower case characters 
        will be converted to upper case, repeated characters are ignored.
        Example of constructing the DNA alphabet:
        >>> alpha = Alphabet('ACGTttga')
        >>> alpha.symbols
        ('A', 'C', 'G', 'T') """
        # Add each symbol to the symbols list, one at a time, and ignore doubles (could use "set" here...)
        _symbols = [] # create a temporary list
        for s in symbolString:
            if not str(s).upper()[0] in _symbols:
                _symbols.append(str(s).upper()[0])
        _symbols.sort() # we put them in alphabetical (one canonical) order
        # OK done extracting, put them in place
        self.symbols = tuple(_symbols); # create the immutable tuple from the extracted list
        self.length = len(self.symbols)
        self.annotations = {}

    def __str__(self):
        return str(self.symbols)
    
    def __len__(self):
        return len(self.symbols)
    
    def __iter__(self):
        return self.symbols.__iter__()
    
    def __getitem__(self, ndx):
        """ Retrieve the symbol(s) at the specified index (or slice of indices) """
        return self.symbols[ndx]
    
    def __contains__(self, sym):
        """ Check if the given symbol is a member of the alphabet. """
        return sym in self.symbols
    
    def index(self, sym):
        """ Retrieve the index of the given symbol in the alphabet. """
        # If the symbol is valid, use the tuple's index function
        if sym in self.symbols:
            syms = self.symbols
            return syms.index(sym)
        else:
            raise RuntimeError('Symbol %s is not indexed by alphabet %s' % (sym, str(self.symbols)))
        
    def __eq__(self, rhs):
        """ Test if the rhs alphabet is equal to ours. """
        if rhs == None:
            return False
        if len(rhs) != len(self):
            return False
        # OK we know they're same size...
        for sym in self.symbols:
            if not sym in rhs:
                return False
        return True

    def isSubsetOf(self, alpha2):
        """ Test if this alphabet is a subset of alpha2. """
        for sym in self.symbols:
            if not alpha2.isValidSymbol(sym):
                return False
        return True
    
    def isSupersetOf(self, alpha2):
        """ Test if this alphabet is a superset of alpha2. """
        return alpha2.isSubsetOf(self)
    
    def annotateSym(self, label, sym, value):
        try:
            lookup = self.annotations[label]
        except KeyError:
            lookup = self.annotations[label] = {}
        lookup[sym] = value
            
    def annotateAll(self, label, symdictOrFilename):
        if isinstance(symdictOrFilename, str): # we assume it is a filename
            fh = open(symdictOrFilename)
            string = fh.read()
            d = {}
            for line in string.splitlines():
                if len(line.strip()) == 0:
                    continue
                sections = line.split()
                symstr, value = sections[0:2]
                for sym in symstr:
                    d[sym] = value
            fh.close()
        else: # we assume it is a dictionary 
            d = symdictOrFilename
        for sym in d:
            self.annotateSym(label, sym, d[sym])
        
    def getAnnotation(self, label, sym):
        try:
            lookup = self.annotations[label]
            return lookup[sym]
        except KeyError:
            return None


""" Below we declare alphabets that are going to be available when 
this module is imported """
Bool_Alphabet = Alphabet('TF')
DNA_Alphabet = Alphabet('ACGT')
DNA_Alphabet_wN = Alphabet('ACGTN')
RNA_Alphabet = Alphabet('ACGU')
Protein_Alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
Protein_Alphabet_Gappy = Alphabet('ACDEFGHIKLMNPQRSTVWY-')
Protein_Alphabet_wX = Protein_wX = Alphabet('ACDEFGHIKLMNPQRSTVWYX')
Protein_Alphabet_wSTOP = Alphabet('ACDEFGHIKLMNPQRSTVWY*')
DSSP_Alphabet = Alphabet('GHITEBSC')
DSSP3_Alphabet = Alphabet('HEC')

predefAlphabets = {'DNA': DNA_Alphabet,
                   'RNA': RNA_Alphabet,
                   'DNAwN': Alphabet('ACGTN'),
                   'RNAwN': Alphabet('ACGUN'),
                   'Protein': Protein_Alphabet,
                   'Protein_Gappy': Protein_Alphabet_Gappy,
                   'ProteinwX': Protein_wX}

class Sequence(object):

    def __init__(self, sequence, alphabet = None, name = '', info = '', **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.sequence = sequence
        self.name = name
        self.alphabet = alphabet

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        out = "%s: %s" % (self.name, self.sequence)
        return out

    def __contains__(self, item):
        return item in self.sequence

    def __getitem__(self, ndx):
        return self.sequence[ndx]

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.sequence == other.sequence
        return False

class Alignment(object):

    def __init__(self, sequences):
        self.seqs = [s for s in sequences]
        self.seqs_dict = dict([(s.name, s) for s in self.seqs])
        self.alignlen = len(sequences[0])
        self.alphabet = self.seqs[0].alphabet

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, ndx):
        return self.seqs[ndx]

    def __contains__(self, item):
        return item.sequence in [s.sequence for s in self.seqs]

    def __str__(self):
        output = ""
        for seq in self.seqs:
            output += "%s\t%s\n" % (seq.name, seq.sequence)
        return output

    def get_column(self, position):
        return [s[position] for s in self.seqs]

    def get_column_probabilities(self, position, pseudo = 0, normalise = True):
        """
        Returns probabilities of each symbol in the alphabet at a position
        """
        colstr = "".join(self.get_column(position))
        # print colstr
        # print position
        # sys.exit()
        if "-" in self.alphabet:
        	cnts = Counter(colstr)
        else:
        	cnts = Counter(colstr.replace("-", ''))
        for sym in self.alphabet:
            if sym not in cnts and pseudo:
                cnts[sym] = pseudo
            else:
                cnts[sym] += pseudo
        total = sum(cnts.values())
        probs = dict()
        if not normalise: return cnts
        for sym, cnt in cnts.items():
            probs[sym] = float(cnt)/total
        return probs

    def get_sequence(self, seq_name):
        try:
            return self.seqs_dict[seq_name]
        except KeyError:
            raise KeyError("Sequence %s was not found in the alignment" % seq_name)

    def write_clustal_file(self, filename):
        """
        Save a Alignment in CLUSTAL format.
        """
        symbolsPerLine = 60
        max_name_length = max(len(seq.name) for seq in self.seqs)
        namelen = 0
        string = ''
        for seq in self.seqs:
            namelen = max(len(seq.name), namelen)
        wholeRows = self.alignlen / symbolsPerLine
        for i in range(wholeRows):
            for j in range(len(self.seqs)):
                string += self.seqs[j].name.ljust(max_name_length) + ' '
                string += self.seqs[j][i * symbolsPerLine:(i + 1) * symbolsPerLine] + '\n'
            string += '\n'
        # Possible last row
        last_row_length = self.alignlen - wholeRows * symbolsPerLine
        if last_row_length > 0:
            for j in range(len(self.seqs)):
                if max_name_length > 0:
                    string += self.seqs[j].name.ljust(max_name_length) + ' '
                string += self.seqs[j][-last_row_length:] + '\n'
        if filename:
            fh = open(filename, 'w')
            # fake header so that clustal believes it
            fh.write('CLUSTAL O(1.2.0) multiple sequence alignment\n\n\n')
            fh.write(string)
            fh.close()
            return
        return string

    def get_ungapped(self):
        """
        Return new alignment with gappy columns removed
        """
        gappy_columns = set()
        for seq in self:
            for i in range(self.alignlen):
                if i not in gappy_columns and seq[i] == "-":
                    gappy_columns.add(i)
        new_seqs = []
        for seq in self:
            content = "".join([seq[i] for i in range(self.alignlen) if i not in gappy_columns])
            new_seqs.append(Sequence(sequence=content, name=seq.name))
        return Alignment(new_seqs)

    def get_ungapped_using_reference(self, seq_name):
        """
        Return a new alignment where gappy columns have been removed using in respect to
        a user specified reference sequence
        :param name: Name of template sequence
        :return:
        """
        template_seq = self.get_sequence(seq_name)
        #Find gapped column indices
        gappy_columns = [i for i in xrange(len(template_seq)) if template_seq[i] == "-"]
        new_seqs = []
        for seq in self:
            content = "".join([seq[i] for i in range(self.alignlen) if i not in gappy_columns])
            new_seqs.append(Sequence(sequence=content, name=seq.name))
        return Alignment(new_seqs)


def read_clustal_file(filename, alpha):
    """
    Read a CLUSTAL Alignment file and return a dictionary of sequence data
    """
    fh = open(filename, 'r')
    names = []
    seqdata = dict()
    data = [line.strip('\n') for line in fh.readlines() if line is not None]
    for line in data:
        if line.startswith('CLUSTAL') or line.startswith('#'):
            continue
        if len(line) == 0:
            continue
        if line[0] == ' ' or '*' in line or ':' in line:
            continue
        sections = line.split()
        name, seqstr = sections[0], "".join(sections[1:])
        names.append(name)
        if seqdata.has_key(name):
            seqdata[name] += seqstr
        else:
            seqdata[name] = seqstr
    sequences = [Sequence(seqstr, name=seqname, alphabet = alpha) for seqname, seqstr in sorted(seqdata.items(),
                                                                            key=lambda x: names.index(x[0]))]
    return Alignment(sequences)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get probabilities for a column')
    parser.add_argument('-i', '--input', help='Input alignment file', required=True)
    # parser.add_argument('-o', '--output', help='Output file', required=False, default="./")
    parser.add_argument('-a', '--alphabet', help='Alphabet, e.g Protein_Alphabet', required=False, default = "Protein")
    parser.add_argument('-c', '--columns', help='List of columns, sep = ",", indexing starts at 1', type = str, required=True)
    parser.add_argument('-p', '--pseudo', help='Pseudo counts', required=False, type = int, default = 0)
    args = parser.parse_args()
    columns = [int(c) for c in args.columns.split(",")]

    the_aln = read_clustal_file(args.input, predefAlphabets[args.alphabet])

    for c in columns:
    	probs = the_aln.get_column_probabilities(position = c-1, pseudo = args.pseudo, normalise = True)
    	counts =the_aln.get_column_probabilities(position = c-1, pseudo = args.pseudo, normalise = False)
    	for sym in the_aln.alphabet:
    		print "{}\t{}\t{}\t{:0.4f}".format(c, sym,counts[sym], probs[sym])

