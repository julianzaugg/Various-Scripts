import sys,os
import argparse



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Take sequences in fasta format and remove annotations")
    
    parser.add_argument('-i', '--input', help='Input file', required=True)
    parser.add_argument('-a', '--annotations_to_remove', help='List of comma-separated annotation indices to remove. E.g\
        "1,2,3" will remove the 1st to 3rd annotations', required=True, type = str)
    parser.add_argument('-d', '--delete_from', help='Remove all annotations from the first index provided (inclusive)', required=False, action='store_true')
    parser.add_argument('-o', '--output', help='Output filename', required=False, default = "processed_out_seqs.txt")

    # temp = ["-i","Blast_RpoH_20000-ed-2.txt", "-a", "2,5", "-d", "-o","test.txt"]
    args = parser.parse_args()
    
    indices = map(lambda x: int(x) - 1, args.annotations_to_remove.split(","))

    #Load sequences
    all_seqs = dict()
    with open(args.input, 'r') as fh:
        data = [line.strip() for line in fh.readlines()]
        for line in data:
            if not line: continue
            if line.startswith(">"):
                cur_name = line[1:]
                cur_name_parts = cur_name.split(" ")
                # print cur_name_parts
                updated_name = []
                for i in xrange(len(cur_name_parts)):
                    if i not in indices:
                        updated_name.append(cur_name_parts[i])
                if args.delete_from:
                    updated_name = updated_name[:indices[0]]
                updated_name = ' '.join(updated_name)
                # print cur_name
                # print
                # print updated_name
                # sys.exit()
                seq_string = ''
                all_seqs[updated_name] = seq_string
                continue
            else:
                all_seqs[updated_name] += line

    # Write out corrected file
    with open(args.output, 'w') as fh:
        for seq_name, seq_string in all_seqs.iteritems():
            print >> fh, ">" + seq_name
            chunks = [seq_string[i:i+60] for i in range(0, len(seq_string), 60)]
            for line in chunks:
                print >> fh, line

