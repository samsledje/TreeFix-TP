import re, sys

# Parse Arguments
#file_name = sys.argv[1]
#divider = sys.argv[2]
#grp = sys.argv[3]

def replace(file, resp, regex=None, split_char=None, group=None, n_groups=None):
    # Build RE
    if resp == 1:
        assert split_char is not None and group is not None and n_groups is not None, "Something went wrong"
        grps = [i+1 for i in range(n_groups)]
        re_string = '>{}{}$'.format('(.*)[{}]'.format(split_char)*(n_groups-1), '(.*)')
        re_compiled = re.compile(re_string, re.M)
        grps.remove(group)
    else:
        re_compiled = re.compile(regex)
        n_groups = 3
        group = 2
        l.remove(2)

    # Read file
    with open(file, 'r') as f:
        file_string = f.read()

    name_iter = re_compiled.finditer(file_string)
    for item in name_iter:
        print(item.group())
        new_string = '{}_{}_{}'.format(item.group(group),item.group(grps[0]), item.group(grps[1]))
        print(new_string)
        #file_string = re.sub(item.group(), '{}_{}_{}'.format(item.group(group),item.group(grps[0]), item.group(grps[1])), file_string)

    #Rewrite file
    with open('{}.formatted'.format(file), 'w+') as f:
        f.write(file_string)

def main():
    # User Interaction
    print("Reformat FASTA file for use with TreeFix-VP")
    print("-"*43)

    while True:
        file_name = input("Enter the name of your FASTA file: ")
        try:
            file = open(file_name, 'r')
            file.close()
            break
        except FileNotFoundError:
            print("Invalid file")
            continue

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
                split_char = '\{}'.format(split_char)    
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
            regex = input("Enter regular expression in the form '(PRE_HOST)(HOST)(POST_HOST)': ")
            try:
                re.compile(regex_string.strip(), re.M)
                break
            except:
                print("Invalid input, please enter a valid regular expression")
                continue
    
    replace(file_name, resp, regex=regex, split_char=split_char, group=group, n_groups=n_groups)

main()
