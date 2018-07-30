TreeFix-VP
https://github.com/samsledje/TreeFix-VP
Mukul Bansal, Samuel Sledzieski, Yi-Chieh Wu,
with libraries contributed by Matthew Rasmussen

=============================================================================
ABOUT

TreeFixVP is a program for reconstructing highly accurate viral phylogenetic trees.
TreeFixVP incorporates host information with sequence data to reconstruct the tree using 
a mixed maximum likelihood / maximum parsimony framework. Given a maximum likelihood phylogeny
and a multiple sequence alignment, TreeFixVP searches among multiple topologies which are supported
equally by sequence data (using the Shimodaiara-Hasegawa test statistic). The selected topology
is the one that minimizes then umber of necessary transmissions.

TreeFix-DTL citation:
Bansal*, Wu*, Alm, Kellis. Improving the Accuracy of Gene Tree Reconstruction
in Prokaryotes: Strategies and Impact. Submitted.

By default, TreeFix-VP uses p-values based on the SH test statistic,
as computed by RAxML.  If you use this default, please also cite

Stamatakis. RAxML-VI-HPC: Maximum Likelihood-based Phylogenetic Analyses
with Thousands of Taxa and Mixed Models. Bioinformatics 22(21):2688-2690, 2006

The original RAxML source code (v7.0.4) is written by Alexandros Stamatakis
and available at http://sco.h-its.org/exelixis/software.html.

This package includes the Python source code of the TreeFix-VP program
(including the TreeFix program and TreeFix-VP wrapper),
modified RAxML source code, the fitch.linux executable,
as well as several library interfaces for Python.


=============================================================================
DETAILS

TreeFixVP uses the Shimodaira-Hasegawa (SH) test statistic with RAxML site-wise likelihoods to compute
p-values for each candidate tree. TreeFixVP scores each statistically equivalent candidate tree using Fitch's algorithm.
Given host labels at the leaf nodes, the Fitch module computes a score equivalent to the
minimum necessary number of transmissions needed to label the internal nodes.

=============================================================================
REQUIREMENTS

Python (2.5.4 or greater)
C compiler (gcc)
SWIG (1.3.29 or greater)
Numpy (1.5.1 or greater)
Scipy (0.7.1 or greater)
Additionally, Python modules are required for computing the p-value for likelihood equivalence.

=============================================================================
Usage: treefixVP [options] <gene tree> ...

TreeFix-VP is a phylogenetic program for improving viral phylogenetic tree
reconstructions using a test statistic for likelihood equivalence and a
transmission aware cost function.

Options:
  Input/Output:
    -A <alignment file extension>, --alignext=<alignment file extension>
                        alignment file extension (default: ".align")
    -o <old tree file extension>, --oldext=<old tree file extension>
                        old tree file extension (default: ".tree")
    -n <new tree file extension>, --newext=<new tree file extension>
                        new tree file extension (default: ".treefix.tree")

  Likelihood Model:
    -m <module for likelihood calculations>, --module=<module for likelihood calculations>
                        module for likelihood calculations (default:
                        "treefix.models.raxmlmodel.RAxMLModel")
    -e <extra arguments to module>, --extra=<extra arguments to module>
                        extra arguments to pass to program

  Likelihood Test:
    -t <test statistic>, --test=<test statistic>
                        test statistic for likelihood equivalence (default:
                        "SH")
    --alpha=<alpha>     alpha threshold (default: 0.05)
    -p <alpha>, --pval=<alpha>
                        same as --alpha

  Transmission Cost Model:
    -M <module for transmission cost calculation>, --smodule=<module for transmission cost calculation>
                        module for transmission cost calculation (default:
                        "treefix.models.fitchmodel.FitchModel")

  Search Options:
    -x <seed>, --seed=<seed>
                        seed value for random generator
    --niter=<# iterations>
                        number of iterations (default: 1000)
    --nquickiter=<# quick iterations>
                        number of subproposals (default: 100)

  Information:
    --version           show program's version number and exit
    -h, --help          show this help message and exit
    -V <verbosity level>, --verbose=<verbosity level>
                        verbosity level (0=quiet, 1=low, 2=medium, 3=high)
    -l <log file>, --log=<log file>
                        log filename.  Use '-' to display on stdout.

TreeFix Written by Yi-Chieh Wu (yjw@mit.edu), Massachusetts Institute of Technology.
(c) 2018. Released under the terms of the GNU General Public License.

Written by Mukul Bansal (mukul.bansal@uconn.edu) and Samuel Sledzieski (samuel.sledzieski@uconn.edu), University of Connecticut.
(c) 2018. Released under the terms of the GNU General Public License.

#=============================================================================
# Examples

See examples/test.sh for an example of how to use TreeFix-VP.
