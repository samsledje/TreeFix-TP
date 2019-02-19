# TreeFix-VP


TreeFix-TP is a program for reconstructing highly accurate infectious disease phylogenetic trees.
TreeFix-TP incorporates host information with sequence data to reconstruct the tree using 
a mixed maximum likelihood / maximum parsimony framework. Given a maximum likelihood phylogeny
and a multiple sequence alignment, TreeFix-TP searches among multiple topologies which are supported
equally by sequence data (using the Shimodaiara-Hasegawa test statistic). The selected topology
is the one that minimizes the number of necessary transmissions.

The input for TreeFix-TP is a multiple sequence alignment in .fasta format, and a maximum likelihood tree in .newick format.

* The latest version of TreeFix-TP can be found [here](https://github.com/samsledje/TreeFix-TP/releases).

* Installation instructions for TreeFix-TP can be found [here](https://github.com/samsledje/TreeFix-TP/blob/master/docs/INSTALL.txt).

* The manual can be found [here](https://github.com/samsledje/TreeFix-TP/blob/master/docs/TreeFix-VP-Manual.pdf), which includes information on requirements which must be installed.

* A test script which runs TreeFix-TP on example data can be found [here](https://github.com/samsledje/TreeFix-TP/tree/master/examples).

* We recommend using a maximum likelihood tree generated using [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html).

* An explanation of the theory behind TreeFix-TP can be found [here](https://github.com/samsledje/TreeFix-TP/blob/master/docs/Bansal_CAME_2017.pdf).

* If you have any problems installing or using TreeFix-TP, please submit an issue [here](https://github.com/samsledje/TreeFix-TP/issues) or email samuel.sledzieski@uconn.edu.
