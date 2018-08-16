# TreeFix-VP


TreeFixVP is a program for reconstructing highly accurate viral phylogenetic trees.
TreeFixVP incorporates host information with sequence data to reconstruct the tree using 
a mixed maximum likelihood / maximum parsimony framework. Given a maximum likelihood phylogeny
and a multiple sequence alignment, TreeFixVP searches among multiple topologies which are supported
equally by sequence data (using the Shimodaiara-Hasegawa test statistic). The selected topology
is the one that minimizes the number of necessary transmissions.

The input for TreeFix-VP is a multiple sequence alignment in .fasta format, and a maximum likelihood tree in .newick format.

* The latest version of TreeFix-VP can be found [here](https://github.com/samsledje/TreeFix-VP/releases).

* Installation instructions for TreeFix-VP can be found [here](https://github.com/samsledje/TreeFix-VP/blob/master/docs/INSTALL.txt).

* The manual can be found [here](https://github.com/samsledje/TreeFix-VP/blob/master/docs/TreeFix-VP-Manual.pdf), which includes information on requirements which must be installed.

* A test script which runs TreeFix-VP on example data can be found [here](https://github.com/samsledje/TreeFix-VP/tree/master/examples).

* We recommend using a maximum likelihood tree generated using [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html).

* An explanation of the theory behind TreeFix-VP can be found [here](https://github.com/samsledje/TreeFix-VP/blob/master/docs/Bansal_CAME_2017.pdf).

* If you have any problems installing or using TreeFix-VP, please submit an issue [here](https://github.com/samsledje/TreeFix-VP/issues) or email samuel.sledzieski@uconn.edu.
