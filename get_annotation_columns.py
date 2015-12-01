#!/usr/bin/python

import argparse


class Annotation:

    def __init__(self, delimiter = "\t", filename=None):
        self.header = []
        self.entries = dict()
        self.delimiter = delimiter
        self.number_of_entries = 0
        if filename:
            self.load_annotations(filename)

    def load_annotations(self, filename):
        """
        Load annotations from a file
        :param filename: filename
        """
        with open(filename, 'r') as fh:
            data = [line.strip().split(self.delimiter) for line in fh.readlines()]
            self.header = data[0]
            self.entries = dict([(c, []) for c in self.header])
            for line in data[1:]:
                for h_idx in range(len(self.header)):
                    self.entries[self.header[h_idx]].append(line[h_idx])
                self.number_of_entries += 1

    def get_column(self, column_name):
        """
        Return entries for a specific column
        :param column_name: Name of header
        :return:
        """
        assert column_name in self.entries, "No column with the header %s in annotation" % column_name
        return self.entries[column_name]

    def add_column(self, column_name, entries):
        if self.number_of_entries != 0:
            assert len(entries) == self.number_of_entries, "Number of entry rows is less than existing annotation"
        elif self.number_of_entries == 0:
            self.number_of_entries = len(entries)
        self.entries[column_name] = entries
        self.header.append(column_name)

    def __str__(self):
        out = "\t".join(self.header) + "\n"
        for i in range(self.number_of_entries):
            line_str = "\t".join([self.entries[column_name][i] for column_name in self.header])
            out += line_str + '\n'
        return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Gets columns from an delimited text file \
                                     and writes them to a new file')
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-o', '--output', help='Output file', required=True)
    parser.add_argument('-l', '--list', help='list of columns to extract',nargs='+', required=True)
    parser.add_argument('-d', '--delimiter', help='delimiter symbol to use, default is "\\t"', required=False, default = "\t")
    parser.add_argument('-v', '--verbose', help='print the table to command line', required=False, action='store_true')
    args = parser.parse_args()

    my_annotation = Annotation(delimiter = args.delimiter, filename = args.input)
    columns_to_grab = []
    for identifier in args.list:
        if identifier.isdigit() and int(identifier) <= len(my_annotation.header):
            columns_to_grab.append(my_annotation.header[int(identifier)])
            continue
        elif identifier in my_annotation.header:
            columns_to_grab.append(identifier)
            continue
        print "Identifier '%s' could not be found. If column number, it is out of range (start at 0)." % (identifier)
    new_annotation = Annotation(delimiter = args.delimiter)
    if args.verbose: print new_annotation
    for column_name in columns_to_grab:
        new_annotation.add_column(column_name, my_annotation.get_column(column_name))
    with open(args.output, 'w') as fh:
        print >> fh, str(new_annotation)
    




