#!/usr/bin/env python3
import re, sys

def rename_fasta(infile, outfile, resp, regex=None, split_char=None, group=None, n_groups=None):
    # Build RE
    if resp == 1:
        assert split_char is not None and group is not None and n_groups is not None, "Something went wrong"
        grps = [i+1 for i in range(n_groups)]
        re_string = '>{}(.*)$'.format("(.*?)[{}]".format(split_char)*(n_groups-1))
        re_compiled = re.compile(re_string, re.M)
        grps.remove(group)
    else:
        re_compiled = re.compile(regex)
        grps = [1,2]
        grps.remove(group)

    # Rename Seqs
    with open(infile, 'r') as f:
        with open(outfile, 'w+') as out:
            for line in f:
                if line.startswith('>'):
                    match = re_compiled.match(line)
                    new_string = '>{}_{}\n'.format(match.group(group), "_".join([match.group(i) for i in grps]))
                    out.write(new_string)
                else:
                    out.write(line)

def rename_dialogue(file_name):
    # User Interaction
    print("Reformat FASTA file for use with TreeFix-VP")
    print("-"*43)

    print("In order to properly calculate transmission cost, each sequence name must be formatted as 'HOST_UNIQUE-STRAIN-ID'. Select how to find the host name within the current sequence names:")
    print("1. Split on character\n2. Use regular expression")
    while True:
        resp = input()
        try:
            resp = int(resp)
        except ValueError:
            print("Invalid input, please enter 1 or 2")
            continue
        if resp == 1 or resp == 2:
            break
        else:
            print("Invalid input, please enter 1 or 2")
            continue

    if resp == 1:
        regex = None
        special_chars = ['.', '^', '$', '*', '+', '?', '{', '}', '[', ']']
        while True:
            split_char = input("Enter split character: ")
            if len(split_char) > 1:
                print("Invalid input, please enter a single character")
                continue
            if split_char in special_chars:
                split_char = '\\{}'.format(split_char)
            break
        while True:
            n_groups = input("Number of groups in sequence name: ")
            try:
                n_groups = int(n_groups)
            except ValueError:
                print("Invalid input")
                continue
            break
        while True:
            group = input("Host name appears in group number: ")
            try:
                group = int(group)
            except ValueError:
                print("Invalid input")
                continue
            if group > n_groups or group <= 0:
                print("Please enter a valid group number")
                continue
            break

    if resp == 2:
        split_char = None
        group = None
        while True:
            regex = input("Enter regular expression with groups for host name and unique ID: ")
            try:
                re.compile(regex.strip(), re.M)
                break
            except:
                print("Invalid input, please enter a valid regular expression")
                continue
        while True:
            group = input("Host name is group 1 or 2?: ")
            try:
                group = int(group)
            except ValueError:
                print("Invalid input")
                continue
            if group == 1 or group == 2:
                break
            else:
                print("Please enter 1 or 2")
                continue
        n_groups = 2

    rename_fasta(file_name, '{}.formatted.seqs'.format(file_name), resp, regex=regex, split_char=split_char, group=group, n_groups=n_groups)
    return "Done"

def short_help():
    print("Usage: ./tree_utils [operation] [sequence file path]")
    print("Operations: '{}'".format("', '".join(ops_dict)))
    sys.exit(1)

def long_help():
    print("Usage: ./tree_utils [operation] [sequence file path]\n")
    print("rename: formats multiple sequence alignment in FASTA format for use with TreeFix-VP")
    print("help: displays this message")
    sys.exit(1)

def main():

    try:
        op = sys.argv[1]
    except:
        short_help()
    if op == "help":
        long_help()
    try:
        tf = sys.argv[2]
    except:
        short_help()

    print("Operation = {}".format(op))
    print("Tree File = {}".format(tf))

    print(ops_dict.get(op, lambda x: 'Invalid Operation')(tf))

if __name__ == "__main__":
    ops_dict = {
        "rename": rename_dialogue,
        "help": long_help
      }
    main()

