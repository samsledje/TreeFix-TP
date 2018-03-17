import sys
from rasmus import treelib, util

#=============================
# input

def add_common_options(parser,
                       infiles=False ,reroot=False,
                       stree=False, smap=False,
                       treeext=False, alignext=False,
                       clade=False):
    """Add common options to parser"""

    if infiles:
        parser.add_option("-i", "--input", dest="input",
                          metavar="<input file>",
                          help="file containing a viral phylogeny")
    if reroot:
        parser.add_option("-r", "--reroot", dest="reroot",
                          action="store_true", default=False,
                          metavar="<reroot tree>",
                          help="set to reroot the input tree")
    if stree:
        parser.add_option("-s", "--stree", dest="stree",
                          metavar="<species tree>",
                          help="gene tree file in newick format")
    if smap:
        parser.add_option("-S", "--smap", dest="smap",
                          metavar="<strain map>",
                          help="strain to host map")
    if treeext:
        parser.add_option("-T","--treeext", dest="treeext",
                          metavar="<tree file extension>",
                          default=".tree",
                          help="tree file extension (default: \".tree\")")
    if alignext:
        parser.add_option("-A","--alignext", dest="alignext",
                          metavar="<alignment file extension>",
                          default=".align",
                          help="alignment file extension (default: \".align\")")
    if clade:
        parser.add_option("-c", "--clade", dest="clade",
                          metavar="<clade file>",
                          help="clade file")

def move_option(parser, opt_str, opt_grp):
    """Move option 'opt_str' from 'parser' to 'opt_grp'"""

    if parser.has_option(opt_str):
        opt = parser.get_option(opt_str)
        parser.remove_option(opt_str)
        opt_grp.add_option(opt)

def check_req_options(parser, options,
                      species=True, clade=True):
    """Check if required options are present"""
    if species and ((not options.stree) or (not options.smap)):
        parser.error("--stree and --smap are required")
    if clade and (not options.clade):
        parser.error("--clade is required")

def get_input_files(parser, options, args):
    """Determine input files from options"""

    infiles = []
    if options.input:
        for arg in options.input:
            if arg == "-":
                infiles.append(sys.stdin)
            else:
                infiles.append(open(arg))

    # determine all input lines
    files = args
    for infile in infiles:
        files.extend(map(lambda fn: fn.rstrip("\n"), infile.readlines()))
    if len(files) == 0:
        parser.error("must specify input file(s)")

    return files

#=============================
# clades

def get_clade(names, stree):
    """Get clade of given species"""

    head = treelib.lca([stree.nodes[name] for name in names])
    assert sorted(names) == sorted(node.name for node in head.leaves())
    return [head] + head.descendants()

def read_clades(filename, stree):
    """Read a clade file"""

    clades = {}
    for toks in util.DelimReader(filename):
        name, sps = toks[0], toks[1].split(',')
        assert all(sp in stree.nodes and stree.nodes[sp].is_leaf() for sp in sps)
        clades[name] = get_clade(sps, stree)
    return clades

def label_clades(gtree, recon, clades):
    """
    Find all subclades by finding all branches separating clade,
    e.g. node is in subtree of clade but parent node is not.
    """

    nodes = {}
    for cname, clade in clades.iteritems():
        for node in gtree.preorder():
            if (node.parent) and \
               (recon[node] in clade) and (recon[node.parent] not in clade):
                nodes[node] = cname
    return nodes
