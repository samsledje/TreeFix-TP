#!/usr/bin/env python

# this program is deprecated - use test_mod.py instead

# python libraries
import os, sys
import optparse

# raxml library
from treefix_raxml import raxml

# scipy library
from scipy.stats import norm

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
    default=".align",
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


def draw_raxml_tree(tr, adef):
    util.tic("Tree to string...")
    treestr = raxml.tree_to_string(tr, adef)
    util.toc()

    util.tic("Drawing tree...")
    T = treelib.parse_newick(treestr)
    T2 = treelib.unroot(T)
    treelib.draw_tree(T2, out=sys.stdout, minlen=5, maxlen=5)
    util.toc()


treefile = args[0]
seqfile = util.replace_ext(treefile, options.treeext, options.alignext)
out = util.open_stream(options.output, "w")

adef = raxml.new_analdef()
raxml.init_adef(adef)
tr = raxml.new_tree()
cmd = "raxmlHPC -t %s -s %s %s" % (treefile, seqfile, options.extra)
raxml.init_program(adef, tr, cmd.split(" "))

util.tic("Optimizing model...")
raxml.optimize_model(adef, tr)
util.toc()

# draw_raxml_tree(tr, adef)

util.tic("Getting parameters for LH...")
bestVector, bestLH, weightSum = raxml.compute_best_LH(tr)
util.log("bestLH: %.3f" % bestLH)
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

    r, w = os.pipe()
    fr, fw = os.fdopen(r, "r"), os.fdopen(w, "w")

    tree.write(out, oneline=True)
    out.write("\n")
    out.flush()
    tree.write(fw, oneline=True)
    fw.write("\n")
    fw.close()

    raxml.read_tree(fr, tr, adef)
    fr.close()

    # draw_raxml_tree(tr, adef)

    util.tic("Computing LH...")
    zscore, Dlnl = raxml.compute_LH(adef, tr, bestLH, weightSum, bestVector)
    util.log("zscore: %.3f, Dlnl: %.3f" % (zscore, Dlnl))
    util.toc()

    if Dlnl <= 0:
        util.log("worse likelihood?: %s" % False)  # better topology (higher likelihood)
    else:
        util.log("worse likelihood?: %s" % True)  # worse topology (higher likelihood)

        # really should just use norm.sf(zscore) < pval, but raxml uses two-sided test
        if norm.sf(zscore) * 2 < options.pval:
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
raxml.delete_best_vector(bestVector)
raxml.delete_analdef(adef)
raxml.delete_tree(tr)
