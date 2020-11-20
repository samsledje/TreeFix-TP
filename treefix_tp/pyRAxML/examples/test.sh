#!/usr/bin/bash
#
# This is an example of how to use TreeFix-RAxML to compute the SH statistic.
#

# Make sure tools are compiled and installed before running the commands in
# this tutorial.  See INSTALL.txt for more information.

# Or you can run from the source directory:

cd ..
python setup.py build_ext --inplace

cd examples
export PYTHONPATH=$PYTHONPATH:../python

#=============================================================================
# Compute the SH statistics for variations on an original gene trees.

# Usage: ./test_mod.py [OPTION] <gene tree>
#
# Options:
#  -h, --help            show this help message and exit
#  -T <tree file extension>, --treeext=<tree file extension>
#                        tree file extension (default: ".tree")
#  -A <alignment file extension>, --alignext=<alignment file extension>
#                        alignment file extension (default: ".align")
#  --niter=<# iterations>
#                        number of iterations (default: 5)
#  -e <extra arguments to initialize RAxML>, --extra=<extra arguments to initialize RAxML>
#                        extra arguments to pass to program
#  -p <p-value>, --pval=<p-value>
#                        p-value threshold (default: 0.05)
#  -o <output file>, --output=<output file>

./test_mod.py \
    -T .nt.raxml.tree \
    -A .nt.align.phylip \
    sim-fungi/0/0.nt.raxml.tree
