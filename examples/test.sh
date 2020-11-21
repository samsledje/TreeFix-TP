#!/bin/bash
#
# This is an example of how to use TreeFix-TP to correct a tree topology.
#

# Make sure tools are compiled and installed before running the commands in
# this tutorial.  See INSTALL.txt for more information.

#=============================================================================
# Compute the corrected phylogeny using RAxML SH statistics and fitch.linux cost model

# show help information
echo TreeFix-TP help:
treefix-tp --help

sleep 2
printf '\n\n'
echo The following dataset was simulated using FAVITES \(https://github.com/niemasd/FAVITES\)
sleep 2
echo RAxML was used to reconstruct a maximum likelihood tree from the multiple sequence alignment
sleep 2
echo The MSA and ML tree are used as the inputs of TreeFix-VP to generate the corrected viral phylogeny
sleep 2

treefix-tp \
    -A .fasta \
    -o .raxml \
    -n .treefix\
    -V1 -l test_TP.log \
    test_TP.raxml

prefix=Score:
echo The transmission cost for the RAxML tree is
foo=$(../bin/fitch.linux test_TP.raxml)
echo "${foo#"$prefix "}"
sleep 1
echo The transmission cost for the TreeFix-TP tree is
foo=$(../bin/fitch.linux test_TP.treefix)
echo "${foo#"$prefix "}"
sleep 1
echo The actual number of tranmissions is
foo=$(../bin/fitch.linux test_TP.true_tree)
echo "${foo#"$prefix "}"
sleep 1
echo More information can be found in test_TP.log
