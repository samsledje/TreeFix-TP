#!/bin/bash
#
# This is an example of how to use TreeFix-VP to correct a tree topology.
#

# Make sure tools are compiled and installed before running the commands in
# this tutorial.  See INSTALL.txt for more information.

#=============================================================================
# Compute the corrected phylogeny using RAxML SH statistics and fitch.linux cost model

# show help information
echo TreeFix-VP help:
sleep 2
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

printf '\n\n'
echo The following dataset was simulated using FAVITES \(https://github.com/niemasd/FAVITES\)
sleep 2
echo RAxML was used to reconstruct a maximum likelihood tree from the multiple sequence alignment
sleep 2
echo The MSA and ML tree are used as the inputs of TreeFix-VP to generate the corrected viral phylogeny
sleep 2

treefixVP \
    -A .fasta \
    -o .raxml \
    -n .treefix\
    -V1 -l test_TP.log \
    test_TP.raxml

prefix=Score:
echo The transmission cost for the RAxML tree is
foo=$(../bin/fitch.linux test_VP.raxml)
echo "${foo#"$prefix "}"
sleep 1
echo The transmission cost for the TreeFix-VP tree is
foo=$(../bin/fitch.linux test_VP.treefix)
echo "${foo#"$prefix "}"
sleep 1
echo The actual number of tranmissions is
foo=$(../bin/fitch.linux test_VP.true_tree)
echo "${foo#"$prefix "}"
sleep 1
echo More information can be found in VP_test.log
