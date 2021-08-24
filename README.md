# TreeFix-TP
[![TreeFix-TP](https://img.shields.io/github/v/release/samsledje/TreeFix-TP)](https://github.com/samsledje/TreeFix-TP/releases)
[![DOI](https://zenodo.org/badge/140891112.svg)](https://zenodo.org/badge/latestdoi/140891112)
[![License](https://img.shields.io/github/license/samsledje/TreeFix-TP)](https://github.com/samsledje/TreeFix-TP/blob/main/LICENSE)



TreeFix-TP is a program for reconstructing highly accurate infectious disease phylogenetic trees.
TreeFix-TP incorporates host information with sequence data to reconstruct the tree using 
a mixed maximum likelihood / maximum parsimony framework. Given a maximum likelihood phylogeny
and a multiple sequence alignment, TreeFix-TP searches among multiple topologies which are supported
equally by sequence data (using the Shimodaiara-Hasegawa test statistic). The selected topology
is the one that minimizes the number of necessary transmissions.

The input for TreeFix-TP is a multiple sequence alignment in .fasta format, and a maximum likelihood tree in .newick format.

* The latest version of TreeFix-TP can be found [here](https://github.com/samsledje/TreeFix-TP/releases).

* Installation instructions for TreeFix-TP can be found [here](https://github.com/samsledje/TreeFix-TP/blob/master/docs/INSTALL.txt).

* The manual can be found [here](https://github.com/samsledje/TreeFix-TP/blob/master/docs/TreeFix-TP-Manual.pdf), which includes information on requirements which must be installed.

* TreeFix-TP is described in *TreeFix-TP: Phylogenetic Error-Correction for Infectious Disease Transmission Network Inference*. Samuel Sledzieski, Chengchen Zhang, Ion Mandoiu, and Mukul S. Bansal. Pacific Symposium on Biocomputing (PSB) 2021: [Proceedings, pages 119-130](https://psb.stanford.edu/psb-online/proceedings/psb21/sledzieski.pdf). A previous version of the manuscript appeared on [bioRxiv](https://www.biorxiv.org/content/10.1101/813931v1).

* A test script which runs TreeFix-TP on example data can be found [here](https://github.com/samsledje/TreeFix-TP/tree/master/examples).

* We recommend using a maximum likelihood tree generated using [RAxML](https://sco.h-its.org/exelixis/web/software/raxml/index.html).

* If you have any problems installing or using TreeFix-TP, please submit an issue [here](https://github.com/samsledje/TreeFix-TP/issues) or email mukul.bansal@uconn.edu.
