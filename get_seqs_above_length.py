#!/usr/bin/python

import argparse

class Sequence(object):

    def __init__(self, sequence, name = '', info = '', **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.sequence = sequence
        self.name = name
        self.info = info

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

    def write_fasta(self):
        """ Write one sequence in FASTA format to a string and return it. """
        fasta = '>' + self.name + ' ' + self.info + '\n'
        data = self.sequence
        nlines = (len(self.sequence) - 1) / 60 + 1
        for i in range(nlines):
            lineofseq = ''.join(data[i*60 : (i+1)*60]) + '\n'
            fasta += lineofseq
        return fasta

class Alignment(object):

    def __init__(self, sequences):
        self.seqs = [s for s in sequences]
        self.seqs_dict = dict([(s.name, s) for s in self.seqs])
        self.alignlen = len(sequences[0])

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


def read_clustal_file(filename):
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
    sequences = [Sequence(seqstr, name=seqname) for seqname, seqstr in sorted(seqdata.items(),
                                                                            key=lambda x: names.index(x[0]))]
    return Alignment(sequences)

def write_fasta_file(filename, seqs):
    """ Write the specified sequences to a FASTA file. """
    fh = open(filename, 'w')
    for seq in seqs:
        fh.write(seq.write_fasta())
    fh.close()


def read_fasta_file(filename, alphabet):
    """
    Read a Fasta file and return a set of Sequence
    {name: seq_string}
    """
    fh = open(filename, 'r')
    seqdata = dict()
    order = []
    # temp = [line for line in fh.readlines() if line is not None]
    data = [line.strip() for line in fh.readlines() if line is not None]

    for line in data:
        if not line: continue
        if line[0] == '>':
            name = line.split()[0][1:]
            order.append(name)
            seqdata[name] = ''
        else:
            seqdata[name] += line
    output = []
    for sname in order:
        sseq = seqdata[sname]

        output.append(Sequence(name = sname, sequence = sseq.strip(), alphabet = alphabet))
    # output = [Sequence(name = sname, sequence = sseq.strip(), alphabet = alphabet) for sname,
    #                     sseq in seqdata.items()]
    return output

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-i', '--input', help='Input fasta file', required=True)
    parser.add_argument('-l', '--length', help='Minimum length', required=True, default = 0, type = int)
    parser.add_argument('-o', '--output', help='Output file', required=False, default="./")
    args = parser.parse_args()
    the_aln = read_fasta_file(args.input, Protein_Alphabet)
    my_seqs = write_fasta_file(args.output, the_aln.seqs)
    out = []
    for seq in my_seqs:
        if len(seq) >= parser.length:
            out.append(seq)
    write_fasta_file(out)


