#!/usr/bin/env python

# python libraries
import os, sys
import optparse

# treefix_raxml library
try:
    import treefix_raxml as raxml
except:
    from os.path import realpath, dirname, join

    sys.path.append(join(realpath(dirname(dirname(__file__))), "python"))
    import treefix_raxml as raxml

# rasmus and compbio libraries
from rasmus import treelib, util
from compbio import phylo

# =============================
# parser

usage = "usage: %prog [options] <gene tree>"
parser = optparse.OptionParser(usage=usage)
parser.add_option(
    "-T",
    "--treeext",
    dest="treeext",
    metavar="<tree file extension>",
    default=".tree",
    help='tree file extension (default: ".tree")',
)
parser.add_option(
    "-A",
    "--alignext",
    dest="alignext",
    metavar="<alignment file extension>",
    default=".align.phylip",
    help='alignment file extension (default: ".align")',
)
parser.add_option(
    "--niter",
    dest="niter",
    metavar="<# iterations>",
    default=5,
    type="int",
    help="number of iterations (default: 5)",
)
parser.add_option(
    "-e",
    "--extra",
    dest="extra",
    metavar="<extra arguments to initialize RAxML>",
    default="-m GTRGAMMA -n test -e 2.0",
    help="extra arguments to pass to program",
)
parser.add_option(
    "-p",
    "--pval",
    dest="pval",
    metavar="<p-value>",
    default=0.05,
    type="float",
    help="p-value threshold (default: 0.05)",
)
parser.add_option("-o", "--output", dest="output", metavar="<output file>", default="-")

options, args = parser.parse_args()

# =============================
# check arguments

if options.niter < 1:
    parser.error("--niter must be >= 1: %d" % options.niter)

if len(args) != 1:
    parser.error("must specify input file")

# =============================
# main file

treefile = args[0]
seqfile = util.replace_ext(treefile, options.treeext, options.alignext)
out = util.open_stream(options.output, "w")

util.tic("Initializing RAXML and optimizing...")
module = raxml.RAxML()
module.optimize_model(treefile, seqfile, options.extra)
util.toc()

tree = treelib.read_tree(treefile)
for node in tree:
    node.dist = 0
    if "boot" in node.data:
        del node.data["boot"]
treehash = phylo.hash_tree(treelib.unroot(tree, newCopy=True))
treehashes = set([treehash])

for i in xrange(options.niter):
    while treehash in treehashes:
        util.log("random spr")
        node1, node2 = phylo.propose_random_spr(tree)
        phylo.perform_spr(tree, node1, node2)
        treehash = phylo.hash_tree(treelib.unroot(tree, newCopy=True))

    treehashes.add(treehash)
    tree.write(out, oneline=True)
    out.write("\n")
    out.flush()

    util.tic("Computing LH...")
    p, Dlnl = module.compute_lik_test(tree)
    util.log("pvalue: %.3f, Dlnl: %.3f" % (p, Dlnl))
    util.toc()

    if Dlnl <= 0:
        util.log("worse likelihood?: %s" % False)  # better topology (higher likelihood)
    else:
        util.log("worse likelihood?: %s" % True)  # worse topology (higher likelihood)

        if p < options.pval:
            util.log(
                "significant?: %s" % True
            )  # statistically significant (significantly worse likelihood)
        else:
            util.log(
                "significant?: %s" % False
            )  # statistically insignificant (equivalent likelihood)

    util.log("\n")

# cleanup
out.close()
