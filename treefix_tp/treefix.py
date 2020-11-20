#!/usr/bin/env python

#
# This code requires underlying modules for
# (1) testing likelihood equivalence (inherits treefix.models.StatModel), and
# (2) computing the species tree aware cost (inherits treefix.models.CostModel).
#

# python libraries
import os, sys, optparse
from io import StringIO

# import random in order to seed since random is used in phylo
# import numpy.random since this is faster than the native random
import random

try:
    from numpy import random as nprnd

    NUMPY = True
except:
    NUMPY = False

# treefix libraries
import treefix_tp
from treefix_tp import common

# rasmus and compbio libraries
from treefix_tp.rasmus import treelib, util, timer
from treefix_tp.compbio import phylo, alignlib

# ==========================================================
# parser

VERSION = treefix_tp.PROGRAM_VERSION_TEXT


def parse_args():
    """parse input arguments"""

    parser = optparse.OptionParser(
        usage="usage: %prog [options] <gene tree> ...",
        version="%prog " + VERSION,
        description="TreeFix is a phylogenetic program for improving gene tree reconstructions using "
        + "a test statistic for likelihood equivalence and a species tree aware cost function. "
        + "See http://compbio.mit.edu/treefix for details.",
        epilog="Written by Yi-Chieh Wu (yjw@mit.edu), Massachusetts Institute of Technology. "
        + "(c) 2011. Released under the terms of the GNU General Public License.",
    )

    grp_io = optparse.OptionGroup(parser, "Input/Output")
    common.add_common_options(
        grp_io, infiles=True, reroot=True, stree=False, smap=False, alignext=True
    )
    grp_io.add_option(
        "-U",
        "--usertreeext",
        dest="usertreeext",
        metavar="<user tree file extension>",
        help="check if user tree is visited in search",
    )
    grp_io.add_option(
        "-o",
        "--oldext",
        dest="oldext",
        metavar="<old tree file extension>",
        default=".tree",
        help='old tree file extension (default: ".tree")',
    )
    grp_io.add_option(
        "-n",
        "--newext",
        dest="newext",
        metavar="<new tree file extension>",
        default=".treefix.tree",
        help='new tree file extension (default: ".treefix.tree")',
    )
    grp_io.add_option(
        "-B",
        "--boottreeext",
        dest="boottreeext",
        metavar="<bootstrap trees file extension>",
        default=".treefix.boot.trees",
        help='bootstrap trees file extension (default: ".treefix.boot.trees")',
    )
    parser.add_option_group(grp_io)

    default_module = "treefix.models.raxmlmodel.RAxMLModel"
    grp_model = optparse.OptionGroup(parser, "Likelihood Model")
    grp_model.add_option(
        "-m",
        "--module",
        dest="module",
        metavar="<module for likelihood calculations>",
        default=default_module,
        help="module for likelihood calculations " + '(default: "%s")' % default_module,
    )
    grp_model.add_option(
        "-e",
        "--extra",
        dest="extra",
        metavar="<extra arguments to module>",
        help="extra arguments to pass to program",
    )
    parser.add_option_group(grp_model)

    grp_test = optparse.OptionGroup(parser, "Likelihood Test")
    grp_test.add_option(
        "-t",
        "--test",
        dest="test",
        metavar="<test statistic>",
        choices=["AU", "NP", "BP", "KH", "SH", "WKH", "WSH"],
        default="SH",
        help='test statistic for likelihood equivalence (default: "SH")',
    )
    grp_test.add_option(
        "--alpha",
        dest="alpha",
        metavar="<alpha>",
        default=0.05,
        type="float",
        help="alpha threshold (default: 0.05)",
    )
    grp_test.add_option(
        "-p", "--pval", dest="alpha", metavar="<alpha>", type="float", help="same as --alpha"
    )
    parser.add_option_group(grp_test)

    default_smodule = "treefix.models.duplossmodel.DupLossModel"
    grp_smodel = optparse.OptionGroup(parser, "Reconciliation Cost Model")
    grp_smodel.add_option(
        "-M",
        "--smodule",
        dest="smodule",
        metavar="<module for species tree aware cost calculations>",
        default=default_smodule,
        help="module for species tree aware cost calculations "
        + '(default: "%s")' % default_smodule,
    )
    grp_smodel.add_option(
        "-E",
        "--sextra",
        dest="sextra",
        metavar="<extra arguments to module>",
        help="extra arguments to pass to program",
    )
    parser.add_option_group(grp_smodel)

    grp_search = optparse.OptionGroup(parser, "Search Options")
    grp_search.add_option(
        "-b",
        "--boot",
        dest="nboot",
        metavar="<# bootstraps>",
        default=1,
        type="int",
        help="number of bootstraps to perform (default: 1)",
    )
    grp_search.add_option(
        "-x",
        "--seed",
        dest="seed",
        metavar="<seed>",
        type="int",
        help="seed value for random generator",
    )
    grp_search.add_option(
        "--niter",
        dest="niter",
        metavar="<# iterations>",
        default=100,
        type="int",
        help="number of iterations (default: 100)",
    )
    grp_search.add_option(
        "--nquickiter",
        dest="nquickiter",
        metavar="<# quick iterations>",
        default=50,
        type="int",
        help="number of subproposals (default: 50)",
    )
    grp_search.add_option(
        "--freconroot",
        dest="freconroot",
        metavar="<fraction reconroot>",
        default=0.05,
        type="float",
        help="fraction of search proposals to reconroot (default: 0.05)",
    )
    grp_search.add_option(
        "--maxtime",
        dest="maxtime",
        metavar="<maximum runtime>",
        type="int",
        help="maximum runtime (per tree) in seconds",
    )
    parser.add_option_group(grp_search)

    grp_info = optparse.OptionGroup(parser, "Information")
    common.move_option(parser, "--version", grp_info)
    common.move_option(parser, "--help", grp_info)
    grp_info.add_option(
        "-V",
        "--verbose",
        dest="verbose",
        metavar="<verbosity level>",
        default="0",
        choices=["0", "1", "2", "3"],
        help="verbosity level (0=quiet, 1=low, 2=medium, 3=high)",
    )
    grp_info.add_option(
        "-l",
        "--log",
        dest="log",
        metavar="<log file>",
        default="-",
        help="log filename.  Use '-' to display on stdout.",
    )
    parser.add_option_group(grp_info)

    grp_debug = optparse.OptionGroup(parser, "Debug")
    grp_debug.add_option(
        "--debug",
        dest="debug",
        metavar="<debug mode>",
        default=0,
        type="int",
        help="debug mode (octal: 0=normal, "
        + "1=skips likelihood test, "
        + "2=skips cost filtering on pool, "
        + "4=computes likelihood for all trees in pool)",
    )
    parser.add_option_group(grp_debug)

    options, args = parser.parse_args()

    # =============================
    # check arguments

    # input gene tree files
    treefiles = common.get_input_files(parser, options, args)

    # required options
    common.check_req_options(parser, options, clade=False, species=False)
    options.verbose = int(options.verbose)

    # debug options
    if options.debug < 0 or options.debug > 7:
        parser.error("--debug must be in {0,...,7}: %d" % options.debug)
    global debug, DEBUG_SKIP_LIK, DEBUG_SKIP_COST, DEBUG_COMPUTE_ALL_LIK
    debug = bin(options.debug)[2:].zfill(3)
    DEBUG_SKIP_LIK = debug[-1] == "1"
    DEBUG_SKIP_COST = debug[-2] == "1"
    DEBUG_COMPUTE_ALL_LIK = debug[-3] == "1"
    if DEBUG_SKIP_LIK and DEBUG_COMPUTE_ALL_LIK:
        parser.error("cannot set debug flag 4 and 1: %d" % options.debug)

    # other options
    if options.alpha < 0 or options.alpha > 1:
        parser.error("--alpha/-p/--pval must be in [0,1]: %.5g" % options.alpha)

    if options.nboot < 1:
        parser.error("-b/--boot must be >= 1: %d" % options.nboot)

    if options.niter < 1:
        parser.error("--niter must be >= 1: %d" % options.niter)

    if options.nquickiter < 1:
        parser.error("--nquickiter must be >= 1: %d" % options.nquickiter)

    if options.freconroot < 0 or options.freconroot > 1:
        parser.error("--freconroot must be in [0,1]: %d" % options.freconroot)

    if options.reroot:
        print("-r/--reroot is deprecated (gene trees are automatically rerooted)", file=sys.stederr)

    return options, treefiles


# ==========================================================
# utilities

# string to binary representation
try:
    if bin(0):
        pass
except NameError as ne:

    def bin(x):
        """
        bin(number) -> string
        Stringifies an int or long in base 2.
        """
        if x < 0:
            return "-" + bin(-x)
        out = []
        if x == 0:
            out.append("0")
        while x > 0:
            out.append("01"[x & 1])
            x >>= 1
        try:
            return "0b" + "".join(reversed(out))
        except NameError as ne2:
            out.reverse()
        return "0b" + "".join(out)


def log_tree(gtree, log, oneline=True, writeDists=False):
    """print tree to log"""
    treeout = StringIO.StringIO()
    if oneline:
        gtree.write(treeout, oneline=oneline)
        log.log("tree: %s\n" % treeout.getvalue())
    else:
        if writeDists:
            treelib.draw_tree(gtree, out=treeout)
        else:
            treelib.draw_tree(gtree, out=treeout, minlen=5, maxlen=5)
        log.log("tree:\n %s\n" % treeout.getvalue())
    treeout.close()


# ==========================================================
# phylogeny functions


def unroot(gtree, newCopy=True):
    """returns unrooted gtree (with internal root always at the same place)"""

    if newCopy:
        gtree = gtree.copy()
    treelib.unroot(gtree, newCopy=False)
    treelib.reroot(
        gtree, gtree.nodes[sorted(gtree.leaf_names())[0]].parent.name, onBranch=False, newCopy=False
    )
    return gtree


# ==========================================================
# special cases


def check_small_tree(gtree, rooted=False):
    """check for small gene tree"""

    # global variables
    global options
    if options.verbose >= 1:
        global log

    # log
    if options.verbose >= 1:
        log.start("Checking for small input gene tree")
        log.log("")

    # small tree?
    genes = gtree.leaf_names()
    if len(genes) <= 2 or (len(genes) == 3 and not rooted):
        # log
        if options.verbose >= 1:
            log.log("tree size <= 2 or == 3 and unrooted -- search skipped\n")
            log.stop()
            log.log("")

        # is small tree
        return True

    # log
    if options.verbose >= 1:
        log.stop()
        log.log("")

    return False


def check_input_tree(gtree, smodule):
    """check if input gene tree achieved minimum cost"""

    # global variables
    global options
    global DEBUG_SKIP_LIK, DEBUG_SKIP_COST, DEBUG_COMPUTE_ALL_LIK
    global gtimer, runtime_prop, runtime_reconroot, runtime_cost, runtime_stat
    if options.verbose >= 1:
        global log

    # log
    if options.verbose >= 1:
        log.start("Checking input gene tree cost")
        log.log("")

    # compute cost
    gtimer.start()
    cost = smodule.compute_cost(gtree)
    runtime_cost += gtimer.stop()

    # input gene tree achieved minimum cost?
    if cost == smodule.mincost:
        # log
        if options.verbose >= 1:
            log.log("input achieved minimum cost: cost\t= %.6g\n" % cost)
            log.stop()
            log.log("")

        # input tree achieved minimum cost
        return True

    # log
    if options.verbose >= 1:
        log.stop()
        log.log("")

    return False


# search routines


def search_landscape(gtree, aln, module, smodule, rooted, seednum=1):
    """search gene tree landscape"""

    # global variables
    global options, seed
    if options.verbose >= 1:
        global log
    global DEBUG_SKIP_LIK, DEBUG_SKIP_COST, DEBUG_COMPUTE_ALL_LIK
    global gtimer, runtime_start, runtime_prop, runtime_reconroot, runtime_cost, runtime_stat
    global treehash0, usertree
    if options.usertreeext:
        global usertreehash, usertreehash_rooted, usertreehash_unrooted, searched_user_rooted0, searched_user_unrooted0

    # work with new copy
    gtree = gtree.copy()
    tree0 = gtree.copy()

    # optimize likelihood module
    if not DEBUG_SKIP_LIK:
        if options.verbose >= 1:
            log.start("Optimizing likelihood model")
        module.optimize_model(gtree, aln)
        if options.verbose >= 1:
            log.stop()
            log.log("")

    # optimize cost module
    if options.verbose >= 1:
        log.start("Optimizing cost model")
    smodule.optimize_model(gtree)
    if options.verbose >= 1:
        log.stop()
        log.log("")

    # log
    gtimer.start()
    cost0 = smodule.compute_cost(gtree)
    runtime_cost += gtimer.stop()

    # log
    if options.verbose >= 1:
        log.log("search: initial")
        log.log("search: cost\t= %.6g" % cost0)
        log_tree(gtree, log)
        log_tree(gtree, log, oneline=False)

    # user tree statistics
    if usertree:
        searched_user_rooted = searched_user_rooted0
        searched_user_unrooted = searched_user_unrooted0

        gtimer.start()
        usercost = smodule.compute_cost(usertree)
        runtime_cost += gtimer.stop()
        if options.verbose >= 1:
            log.log("user: cost\t= %.6g" % usercost)
            gtimer.start()
            userpval, userDlnl = module.compute_lik_test(usertree, options.test)
            runtime_stat += gtimer.stop()
            log.log("user: pval\t= %.6g" % userpval)
            log.log("user: Dlnl\t= %.6g" % userDlnl)
            log.log("\n")

    # initialize min values
    mintree, mincost, minpval, minDlnl = tree0, cost0, 1, 0

    # search functions
    search = phylo.TreeSearchMix(gtree)
    search.add_proposer(phylo.TreeSearchNni(gtree), 0.5)
    search.add_proposer(phylo.TreeSearchSpr(gtree), 0.5)
    uniques = set([treehash0])

    # do search
    nproposals = 0
    nuniques = 0
    npools = 0
    nemptypools = 0
    ndiffrecon = 0
    nrecon = 0
    for i in xrange(options.niter):  # outer search
        # end early if maxtime reached
        if (options.maxtime is not None) and (timer.time.time() - runtime_start > options.maxtime):
            if options.verbose >= 1:
                log.log("search: break (maxtime reached)\n")
            break

        # store initial search tree
        if rooted:
            treehash1 = phylo.hash_tree(search.tree)
        else:
            treehash1 = phylo.hash_tree(unroot(search.tree, newCopy=True))

        # random values -- have to reseed in case module has a random number generator
        random.seed(seed + i * 1024 * seednum)
        if NUMPY:
            nprnd.seed(seed + i * 1024 * seednum)
            randvec = nprnd.random(2 * options.nquickiter)
        else:
            randvec = [random.random() for _ in xrange(2 * options.nquickiter)]

        # search
        ntrees = 0
        pool = {}
        mincost_pool = mincost

        # note that reconroot is NOT propagated through the subproposals - doing so messes up the unique filter
        for j in xrange(options.nquickiter):  # inner search
            # propose tree
            gtimer.start()
            tree = search.propose()
            runtime_prop += gtimer.stop()

            # only allow unique proposals but if I have been rejecting too much allow some non-uniques through
            treehash = phylo.hash_tree(tree)
            if treehash in uniques and ntrees >= 0.1 * j:
                if options.verbose >= 3:
                    log.log("prescreen: iter %d" % j)
                    log.log("prescreen: revert")
                    log.log("")
                search.revert()
                continue

            # save tree
            nproposals += 1
            gtree = tree.copy()
            if treehash not in uniques:
                uniques.add(treehash)
                nuniques += 1
            ntrees += 1

            if options.verbose >= 3:
                log.log("prescreen: iter %d" % j)

            # reconroot (some percentage of the time depending on options.freconroot)
            if randvec[j + options.nquickiter] < options.freconroot:
                gtimer.start()
                gtree, cost = smodule.recon_root(gtree, newCopy=False, returnCost=True)
                runtime_reconroot += gtimer.stop()

                # did reconroot change the tree?
                nrecon += 1
                if phylo.hash_tree(gtree) == treehash:
                    if options.verbose > 3:
                        log.log("prescreen: recon\t= unchanged")
                else:
                    ndiffrecon += 1
                    if options.verbose > 3:
                        log.log("prescreen: recon\t= changed")
            else:
                gtimer.start()
                cost = smodule.compute_cost(gtree)
                runtime_cost += gtimer.stop()

            # log
            if options.verbose >= 3:
                log.log("prescreen: cost\t= %.6g" % cost)
                log_tree(gtree, log)
                log.log("")

            # store to pool if unique
            treehash_rooted = phylo.hash_tree(gtree)
            treehash_unrooted = phylo.hash_tree(unroot(gtree, newCopy=True))
            treehash = treehash_rooted if rooted else treehash_unrooted
            if treehash != treehash1 and treehash not in pool:
                pool[treehash] = (gtree, cost, j)
            if (
                (usertree)
                and (not searched_user_rooted)
                and (treehash_rooted == usertreehash_rooted)
            ):
                searched_user_rooted = True
            if (
                (usertree)
                and (not searched_user_unrooted)
                and (treehash_unrooted == usertreehash_unrooted)
            ):
                searched_user_unrooted = True

            # update mincost of pool and decide how to continue proposals from here
            if cost < mincost_pool:
                # make more proposals off this one
                mincost_pool = cost
            elif cost == mincost_pool:
                # flip a coin to decide whether to start from original or new proposal
                if randvec[j] < 0.5:
                    search.revert()
            else:
                # start from new proposal 10% of the time
                if randvec[j] < 0.9:
                    search.revert()

        # remove trees with higher costs
        if options.verbose >= 2:
            log.log("pool: size\t= %d" % len(pool))
        pool = pool.values()
        if DEBUG_SKIP_COST:
            fpool = pool
            nfpool = len(fpool)
        else:
            fpool = filter(lambda gtree, cost, ndx: cost <= mincost and cost < cost0, pool)
            nfpool = len(fpool)
            if options.verbose >= 2:
                log.log("pool: filtered size\t= %d" % nfpool)
        if options.verbose >= 2:
            log.log("")
        fpool.sort(key=lambda x: x[1])

        # propose a tree from the pool with minimum cost that passes threshold
        if DEBUG_SKIP_LIK:
            reject = False
            if nfpool > 0:
                mintree, mincost, minpval, minDlnl = fpool[0][0], fpool[0][1], 1, 0
        else:
            reject = True
            for j, (gtree, cost, ndx) in enumerate(fpool):
                gtimer.start()
                pval, Dlnl = module.compute_lik_test(gtree, options.test)
                runtime_stat += gtimer.stop()

                if (pval < options.alpha) or (cost == mincost and Dlnl > minDlnl):
                    # reject topology because
                    # (1) significantly worse topology or (2) worse Dlnl (smaller Dlnl is better)
                    if options.verbose >= 2:
                        log.log("pool: iter %d (%d)" % (j, ndx))
                        log.log("pool: reject")
                        log.log("pool: cost\t= %.6g" % cost)
                        log.log("pool: pval\t= %.6g" % pval)
                        log.log("pool: Dlnl\t= %.6g" % Dlnl)
                        log_tree(gtree, log)
                        log.log("")
                else:
                    # accept topology
                    reject = False
                    mintree, mincost, minpval, minDlnl = gtree, cost, pval, Dlnl
                    break

        # debug
        if DEBUG_COMPUTE_ALL_LIK:
            dpool = (fpool[j + 1 :] if j + 1 < nfpool else []) + [x for x in pool if x not in fpool]
            for (gtree, cost, ndx) in dpool:
                gtimer.start()
                pval, Dlnl = module.compute_lik_test(gtree, options.test)
                runtime_stat += gtimer.stop()

                if options.verbose >= 2:
                    log.log("pool: iter (%d)" % ndx)
                    log.log("pool: cost\t= %.6g" % cost)
                    log.log("pool: pval\t= %.6g" % pval)
                    log.log("pool: Dlnl\t= %.6g" % Dlnl)
                    log_tree(gtree, log)
                    log.log("")

        # reset search and log
        search.reset()
        search.set_tree(mintree.copy())
        npools += 1
        if nfpool == 0:
            nemptypools += 1
            if options.verbose >= 1:
                log.log("search: iter %d" % i)
                log.log("search: empty pool")
                log.log("")
            continue
        if reject:
            if options.verbose >= 1:
                log.log("search: iter %d" % i)
                log.log("search: reject")
                log.log("")
            continue
        if options.verbose >= 1:
            log.log("search: iter %d" % i)
            log.log("search: accept")
            log.log("search: cost\t= %.6g" % mincost)
            log.log("search: pval\t= %.6g" % minpval)
            log.log("search: Dlnl\t= %.6g" % minDlnl)
            log_tree(mintree, log)
            log_tree(mintree, log, oneline=False)
            log.log("")

        # end early if tree with minimum cost found
        if (smodule.mincost != -util.INF) and (mincost == smodule.mincost):
            if options.verbose >= 1:
                log.log("search: break (mincost reached)\n")
            break

    # has the tree changed?
    if rooted:
        treehash = phylo.hash_tree(mintree)
    else:
        # do a final reconroot
        gtimer.start()
        mintree, mincost = smodule.recon_root(mintree, newCopy=False, returnCost=True)
        runtime_reconroot += gtimer.stop()

        treehash = phylo.hash_tree(unroot(mintree, newCopy=True))

    # output search statistics
    if options.verbose >= 1:
        log.log("search: final")
        if usertree:
            log.log("search: visited user tree\t= %s" % ("yes" if searched_user_rooted else "no"))
            log.log(
                "search: visited (unrooted) user tree\t= %s"
                % ("yes" if searched_user_unrooted else "no")
            )
            log.log("search: equal user tree\t= %s" % ("yes" if treehash == usertreehash else "no"))
        if treehash == treehash0:
            log.log("search: changed\t= no")
        else:
            log.log("search: changed\t= yes")
            log.log("search: init cost\t= %.6g" % cost0)
            log.log("search: final cost\t= %.6g" % mincost)
            if not DEBUG_SKIP_LIK:
                log.log("search: pval\t= %.6g" % minpval)
                log.log("search: Dlnl\t= %.6g" % minDlnl)
            log_tree(mintree, log)
            log_tree(mintree, log, oneline=False)

        log.log("")
        log.log("num proposals:\t%d" % nproposals)
        log.log("num unique proposals:\t%d" % nuniques)
        if options.verbose >= 2:
            log.log("")
            log.log("num pools:\t%d" % npools)
            log.log("num empty pools:\t%d" % nemptypools)
            if npools == 0:
                log.log("empty pool rate:\tINF")
            else:
                log.log("empty pool rate:\t%f" % (float(nemptypools) / npools))
        log.log("")
        log.log("num reconroot:\t%d" % nrecon)
        log.log("num diff reconroot:\t%d" % ndiffrecon)
        if nrecon == 0:
            log.log("diff reconroot rate:\tINF")
        else:
            log.log("diff reconroot rate:\t%f" % (float(ndiffrecon) / nrecon))
        log.log("")

    # return optimum
    return mintree


# ==========================================================
# main


def main(options, treefiles):
    """main"""

    # global variables
    # global options, seed
    global seed
    global DEBUG_SKIP_LIK, DEBUG_SKIP_COST, DEBUG_COMPUTE_ALL_LIK
    global gtimer, runtime_start, runtime_prop, runtime_reconroot, runtime_cost, runtime_stat
    global treehash0, usertree

    # debug options
    if options.debug < 0 or options.debug > 7:
        parser.error("--debug must be in {0,...,7}: %d" % options.debug)
    global debug, DEBUG_SKIP_LIK, DEBUG_SKIP_COST, DEBUG_COMPUTE_ALL_LIK
    debug = bin(options.debug)[2:].zfill(3)
    DEBUG_SKIP_LIK = debug[-1] == "1"
    DEBUG_SKIP_COST = debug[-2] == "1"
    DEBUG_COMPUTE_ALL_LIK = debug[-3] == "1"
    if DEBUG_SKIP_LIK and DEBUG_COMPUTE_ALL_LIK:
        parser.error("cannot set debug flag 4 and 1: %d" % options.debug)

    # other options
    if options.alpha < 0 or options.alpha > 1:
        parser.error("--alpha/-p/--pval must be in [0,1]: %.5g" % options.alpha)

    if options.nboot < 1:
        parser.error("-b/--boot must be >= 1: %d" % options.nboot)

    if options.niter < 1:
        parser.error("--niter must be >= 1: %d" % options.niter)

    if options.nquickiter < 1:
        parser.error("--nquickiter must be >= 1: %d" % options.nquickiter)

    if options.freconroot < 0 or options.freconroot > 1:
        parser.error("--freconroot must be in [0,1]: %d" % options.freconroot)

    if options.reroot:
        print("-r/--reroot is deprecated (gene trees are automatically rerooted)", file=sys.stederr)

    # parse arguments
    # options, treefiles = parse_args()
    options.verbose = int(options.verbose)
    if options.verbose >= 1:
        global log

    # global timer
    gtimer = timer.Timer()

    # import likelihood module
    if DEBUG_SKIP_LIK:
        rooted = False
    else:
        folder = ".".join(options.module.split(".")[:-1])
        if folder != "":
            exec("import %s" % folder, globals(), locals())
        module = eval("%s(options.extra)" % options.module)
        rooted = module.rooted

    # import cost module
    folder = ".".join(options.smodule.split(".")[:-1])
    if folder != "":
        exec("import %s" % folder, globals(), locals())
    smodule = eval("{}()".format(options.smodule))

    # log file
    if options.verbose >= 1:
        if options.log == "-":
            log = timer.globalTimer()
        else:
            outlog = util.open_stream(options.log, "w")
            log = timer.Timer(outlog)

        log.log("TreeFix version: %s\n" % VERSION)
        log.log("TreeFix executed with the following arguments:")
        log.log(
            "%s %s\n"
            % (
                os.path.basename(sys.argv[0]),
                " ".join(map(lambda x: x if x.find(" ") == -1 else '"%s"' % x, sys.argv[1:])),
            )
        )

        log.log("Likelihood Model:\t%s (%s)" % (options.module, module.VERSION))
        log.log("Reconciliation Model:\t%s (%s)" % (options.smodule, smodule.VERSION))
        log.log("\n\n")

        if DEBUG_SKIP_LIK:
            log.log("debug: skip likelihood test")
        if DEBUG_SKIP_COST:
            log.log("debug: skip cost requirement")
        if DEBUG_COMPUTE_ALL_LIK:
            log.log("debug: compute likelihoods for all trees in pool")
        if any(map(lambda x: x == "1", debug)):
            log.log("\n")

    # process genes trees
    for treefile in treefiles:
        # seed random generator
        if options.seed:
            seed = options.seed
        else:
            seed = int(timer.time.time())
        random.seed(seed)
        if NUMPY:
            nprnd.seed(seed)

        # runtime statistics
        runtime_start = timer.time.time()
        runtime_prop = 0
        runtime_reconroot = 0
        runtime_cost = 0
        runtime_stat = 0

        # start log
        if options.verbose >= 1:
            log.start("Working on file '%s'" % treefile)
            log.log("random seed: %s\n" % seed)

        # setup files
        alnfile = util.replace_ext(treefile, options.oldext, options.alignext)
        outfile = util.replace_ext(treefile, options.oldext, options.newext)

        # read input tree
        try:
            gtrees = treelib.read_trees(treefile)
        except:
            print('ERROR: problem reading gene tree: "%s"' % treefile,file=sys.stderr)
            return 1
        if len(gtrees) != 1:
            print("ERROR: treefile contains multiple trees: %s" % treefile, file=sys.stederr)
            return 1
        gtree = gtrees[0]

        # remove bootstraps and dists if present
        for node in gtree:
            node.dist = 0
            if "boot" in node.data:
                del node.data["boot"]
        if "boot" in gtree.default_data:
            del gtree.default_data["boot"]

        # optimize cost module
        if options.verbose >= 1:
            log.start("Optimizing cost model")
        smodule.optimize_model(gtree)
        if options.verbose >= 1:
            log.stop()
            log.log("")

        # reroot input tree
        gtree = smodule.recon_root(gtree, newCopy=False)

        # special cases -- no need to search
        #   small gene tree
        #   input gene tree achieved minimum cost
        if check_small_tree(gtree, rooted=rooted) or check_input_tree(gtree, smodule):
            # output tree
            gtree.write(outfile, oneline=False)

            # add bootstraps
            if options.nboot > 1:
                phylo.add_bootstraps(gtree, [gtree], rooted=True)

            # log
            if options.verbose >= 1:
                log.stop()
                log.log("")

            # skip rest of algorithm
            continue

        # read alignment
        try:
            aln = alignlib.fasta.read_fasta(alnfile)
        except:
            print('ERROR: problem reading alignment: "%s"' % alnfile, file=sys.stderr)
            return 1
        if set(aln) != set(gtree.leaf_names()):
            print("ERROR: gene tree and alignment contain different genes", file=sys.stderr)
            return 1
        alnlen = len(aln.values()[0])

        usertree = None

        # store input tree and whether input tree matches user tree
        treehash0_rooted = phylo.hash_tree(gtree)
        treehash0_unrooted = phylo.hash_tree(unroot(gtree, newCopy=True))
        treehash0 = treehash0_rooted if rooted else treehash0_unrooted
        if usertree:
            if set(usertree.leaf_names()) != set(gtree.leaf_names()):
                print("ERROR: gene tree and user tree contain different genes", file=sys.stderr)
                return 1
            searched_user_rooted0 = treehash0_rooted == usertreehash_rooted
            searched_user_unrooted0 = treehash0_unrooted == usertreehash_unrooted

        # setup output options.nboot == 1:
            boot = False
        else:
            boot = True
            boottreefile = util.replace_ext(treefile, options.oldext, options.boottreeext)
            out = util.open_stream(boottreefile, "w")

        # main algorithm
        for bootnum in xrange(options.nboot):  # bootstrap search
            # end early if maxtime reached
            if (options.maxtime is not None) and (
                timer.time.time() - runtime_start > options.maxtime
            ):
                if options.verbose >= 1:
                    log.log("boot: break (maxtime reached)\n")
                break

            # get bootstrapped alignment (equal number of columns, sample with replacement)
            if boot:
                if options.verbose >= 1:
                    log.start("boot: %d" % bootnum)
                    log.log("")

                random.seed(seed + bootnum * 4096)
                if NUMPY:
                    nprnd.seed(seed + bootnum * 4096)
                    cols = nprnd.randint(alnlen, size=alnlen)
                else:
                    cols = [random.randint(0, alnlen - 1) for _ in xrange(alnlen)]
                baln = alignlib.subalign(aln, cols)
            else:
                baln = aln

            # search
            mintree = search_landscape(
                gtree, baln, module, smodule, rooted, seednum=bootnum + boot + 1
            )

            # output bootstrap tree
            if boot:
                mintree.write(out, oneline=True)
                out.write("\n")
                if options.verbose >= 1:
                    log.stop()
                    log.log("")

        # write optimal tree or get optimal tree with bootstrap support
        if not boot:
            mintree.write(outfile)
        else:
            out.close()
            if options.verbose >= 1:
                log.start("boot: final")
                log.log("")

            # search
            mintree = search_landscape(gtree, aln, module, smodule, rooted)

            # add bootstraps
            phylo.add_bootstraps(mintree, treelib.iter_trees(boottreefile), rooted=True)

            # log final tree with bootstraps
            log_tree(mintree, log)
            log_tree(mintree, log, oneline=False)

            # output tree
            mintree.write(outfile)
            if options.verbose >= 1:
                log.stop()
                log.log("")

        # output runtime statistics
        if options.verbose >= 1:
            log.log("proposal runtime:\t%f" % runtime_prop)
            log.log("reconroot runtime:\t%f" % runtime_reconroot)
            log.log("cost runtime:\t%f" % runtime_cost)
            log.log("statistic runtime:\t%f" % runtime_stat)

        # stop log
        if options.verbose >= 1:
            log.stop()
            log.log("\n\n")

    # close log
    if options.verbose >= 1 and options.log != "-":
        outlog.close()


# main function
if __name__ == "__main__":
    sys.exit(main())
