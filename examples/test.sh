#!/usr/bin/bash
#
# This is an example of how to use TreeFix-VP to correct a tree topology.
#

# Make sure tools are compiled and installed before running the commands in
# this tutorial.  See INSTALL.txt for more information.

# Or you can run from the source directory:

cd ..
python setup.py build_ext --inplace

# if using Linux
cp bin/ranger-dtl-U.linux bin/ranger-dtl-U
# if using Mac
cp bin/ranger-dtl-U.mac bin/ranger-dtl-U
# if using Windows (Cygwin)
cp bin/ranger-dtl-U.exe bin/ranger-dtl-U

cd examples
export PATH=$PATH:../bin
export PYTHONPATH=$PYTHONPATH:../python

#=============================================================================
# Compute the corrected gene tree using RAxML SH statistics and ranger-dtl-U cost model

# show help information
treefixVP -h

# Options:
#   Input/Output:
#     -i <input file>, --input=<input file>
#                         file containing a viral phylogeny
#     -A <alignment file extension>, --alignext=<alignment file extension>
#                         alignment file extension (default: ".align")
#     -o <old tree file extension>, --oldext=<old tree file extension>
#                         old tree file extension (default: ".tree")
#     -n <new tree file extension>, --newext=<new tree file extension>
#                         new tree file extension (default: ".treefix.tree")

#   Likelihood Model:
#     -m <module for likelihood calculations>, --module=<module for likelihood calculations>
#                         module for likelihood calculations (default:
#                         "treefix.models.raxmlmodel.RAxMLModel")
#     -e <extra arguments to module>, --extra=<extra arguments to module>
#                         extra arguments to pass to program

#   Likelihood Test:
#     -t <test statistic>, --test=<test statistic>
#                         test statistic for likelihood equivalence (default:
#                         "SH")
#     --alpha=<alpha>     alpha threshold (default: 0.05)
#     -p <alpha>, --pval=<alpha>
#                         same as --alpha

#   Transmission Cost Model:
#     -M <module for transmission cost calculation>, --smodule=<module for transmission cost calculation>
#                         module for transmission cost calculation (default:
#                         "treefix.models.fitchmodel.FitchModel")

#   Search Options:
#     -x <seed>, --seed=<seed>
#                         seed value for random generator
#     --niter=<# iterations>
#                         number of iterations (default: 1000)
#     --nquickiter=<# quick iterations>
#                         number of subproposals (default: 100)

#   Information:
#     --version           show program's version number and exit
#     -h, --help          show this help message and exit
#     -V <verbosity level>, --verbose=<verbosity level>
#                         verbosity level (0=quiet, 1=low, 2=medium, 3=high)
#     -l <log file>, --log=<log file>
#                         log filename.  Use '-' to display on stdout.

treefixVP \
    -A .align \
    -o .tree.old \
    -n .tree.new\
    -V1 -l VP_test.log \
    VP_test.tree.old

#=============================================================================
# Clean up

rm VP_test{.new.tree,.log}
