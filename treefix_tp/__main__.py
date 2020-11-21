#!/usr/bin/env python

# wrapper script around TreeFix for use with fitch reconciliation model

# python libraries
import os, sys, optparse, argparse
import subprocess

# treefix libraries
from treefix_tp import common
from treefix_tp import treefix
from treefix_tp import PROGRAM_VERSION_TEXT

#==========================================================
# parser

VERSION = PROGRAM_VERSION_TEXT
#commands.getoutput("treefix --version").rsplit()[-1]

def build_parser():
    """parse input arguments"""

    parser = optparse.OptionParser(
        usage = "usage: %prog [options] <gene tree> ...",

        version = "%prog " + VERSION,

        description =
             "TreeFix-TP is a phylogenetic program for improving viral phylogenetic tree reconstructions using " +
             "a test statistic for likelihood equivalence and a transmission aware cost function. " +
             "See http://github.com/samsledje/TreeFix-TP for details.",

        epilog =
             "Written by Samuel Sledzieski (samuel.sledzieski@uconn.edu) and Mukul Bansal (mukul.bansal@uconn.edu), University of Connecticut. "+
             "(c) 2018. Released under the terms of the GNU General Public License."+
             "                              "+
             "Treefix written by Yi-Chieh Wu (yjw@mit.edu), Massachusetts Institute of Technology. " +
             "(c) 2011. Released under the terms of the GNU General Public License. ")

    grp_io = optparse.OptionGroup(parser, "Input/Output")
    grp_io.add_option("-A","--alignext", dest="alignext",
                      metavar="<alignment file extension>",
                      default=".align",
                      help="alignment file extension (default: \".align\")")
    grp_io.add_option("-o", "--oldext", dest="oldext",
                      metavar="<old tree file extension>",
                      default=".tree",
                      help="old tree file extension (default: \".tree\")")
    grp_io.add_option("-n", "--newext", dest="newext",
                      metavar="<new tree file extension>",
                      default=".treefix.tree",
                      help="new tree file extension (default: \".treefix.tree\")")
    grp_io.add_option("-r", "--reroot", dest="reroot",
                      action="store_true", default=False,
                      metavar="<reroot tree>",
                      help="set to reroot the input tree")
    parser.add_option_group(grp_io)

    default_module = "treefix_tp.models.raxmlmodel.RAxMLModel"
    grp_model = optparse.OptionGroup(parser, "Likelihood Model")
    grp_model.add_option("-m", "--module", dest="module",
                         metavar="<module for likelihood calculations>",
                         default=default_module,
                         help="module for likelihood calculations " +\
                              "(default: \"%s\")" % default_module)
    grp_model.add_option("-e", "--extra", dest="extra",
                         metavar="<extra arguments to module>",
                         help="extra arguments to pass to program")
    parser.add_option_group(grp_model)

    grp_test = optparse.OptionGroup(parser, "Likelihood Test")
    grp_test.add_option("-t", "--test", dest="test",
                        metavar="<test statistic>",
                        choices=["AU", "NP", "BP", "KH", "SH", "WKH", "WSH"],
                        default="SH",
                        help="test statistic for likelihood equivalence (default: \"SH\")")
    grp_test.add_option("--alpha", dest="alpha",
                        metavar="<alpha>",
                        default=0.05, type="float",
                        help="alpha threshold (default: 0.05)")
    grp_test.add_option("-p", "--pval", dest="alpha",
                        metavar="<alpha>",
                        type="float",
                        help="same as --alpha")
    parser.add_option_group(grp_test)

    default_smodule = "treefix_tp.models.fitchmodel.FitchModel"
    grp_smodel = optparse.OptionGroup(parser, "Transmission Cost Model")
    grp_smodel.add_option("-M", "--smodule", dest="smodule",
                          metavar="<module for transmission cost calculation>",
                          default=default_smodule,
                          help="module for transmission cost calculation " +\
                               "(default: \"%s\")" % default_smodule)
    parser.add_option_group(grp_smodel)

    grp_search = optparse.OptionGroup(parser, "Search Options")
    grp_search.add_option("-b", "--boot", dest="nboot",
                          metavar="<# bootstraps>",
                          default=1, type="int",
                          help="number of bootstraps to perform (default: 1)")
    grp_search.add_option("-x", "--seed", dest="seed",
                          metavar="<seed>",
                          type="int",
                          help="seed value for random generator")
    grp_search.add_option("--niter", dest="niter",
                          metavar="<# iterations>",
                          default=100, type="int",
                          help="number of iterations (default: 100)")
    grp_search.add_option("--nquickiter", dest="nquickiter",
                          metavar="<# quick iterations>",
                          default=50, type="int",
                          help="number of subproposals (default: 50)")
    grp_search.add_option("--freconroot", dest="freconroot",
                          metavar="<fraction reconroot>",
                          default=0.05, type="float",
                          help="fraction of search proposals to reconroot (default: 0.05)")
    grp_search.add_option("--maxtime", dest="maxtime",
                          metavar="<maximum runtime>",
                          type="int",
                          help="maximum runtime (per tree) in seconds")
    parser.add_option_group(grp_search)


    grp_info = optparse.OptionGroup(parser, "Information")
    common.move_option(parser, "--version", grp_info)
    common.move_option(parser, "--help", grp_info)
    grp_info.add_option("-V", "--verbose", dest="verbose",
                        metavar="<verbosity level>",
                        default="0", choices=["0","1","2","3"],
                        help="verbosity level (0=quiet, 1=low, 2=medium, 3=high)")
    grp_info.add_option("-l", "--log", dest="log",
                        metavar="<log file>",
                        default="-",
                        help="log filename.  Use '-' to display on stdout.")
    parser.add_option_group(grp_info)

    grp_debug = optparse.OptionGroup(parser, "Debug")
    grp_debug.add_option("--debug", dest="debug",
                         metavar="<debug mode>",
                         default=0, type="int",
                         help="debug mode (octal: 0=normal, " +\
                              "1=skips likelihood test, " +\
                              "2=skips cost filtering on pool, " +\
                              "4=computes likelihood for all trees in pool)")
    parser.add_option_group(grp_debug)

    return parser

def main():

    # parse arguments
    parser = build_parser()
    options, args = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    if len(args) == 0:
        parser.print_help()
        sys.exit(1)

    # required options
    common.check_req_options(parser, options, clade=False, species=False)
    options.verbose = int(options.verbose)

    # debug options
    if options.debug < 0 or options.debug > 7:
        parser.error("--debug must be in {0,...,7}: %d" % options.debug)
    global debug, DEBUG_SKIP_LIK, DEBUG_SKIP_COST, DEBUG_COMPUTE_ALL_LIK
    debug = bin(options.debug)[2:].zfill(3)
    DEBUG_SKIP_LIK = (debug[-1] == "1")
    DEBUG_SKIP_COST = (debug[-2] == "1")
    DEBUG_COMPUTE_ALL_LIK = (debug[-3] == "1")
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
        print("-r/--reroot is deprecated (gene trees are automatically rerooted)",file=sys.stederr)

    options.boottreeext = ".treefix.boot.trees"
    options.usertreeext = None

    #log
    print("TreeFix-TP version: %s\n" % VERSION)
    print("TreeFix-TP called as:")
    print(' '.join(sys.argv))

    treefix.treefix(options, args)

# main function
if __name__ == "__main__":
    main()
