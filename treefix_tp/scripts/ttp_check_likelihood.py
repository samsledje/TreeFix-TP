import sys
from treefix_tp.models.raxmlmodel import RAxMLModel
from treefix_tp.rasmus.treelib import read_trees
from treefix_tp.compbio.fasta import read_fasta

def main():
    try:
        alnFile, oldTreFile, newTreFile = sys.argv[1:]
    except ValueError:
        print('usage: python ttp-check-likelihood.py [alignment] [old tree file] [new tree file]')
        sys.exit(0)

    aln = read_fasta(alnFile)
    oldT = read_trees(oldTreFile)[0]
    newT = read_trees(newTreFile)[0]

    module = RAxMLModel(None)
    module.optimize_model(oldT, aln)
    pval, Dlnl = module.compute_lik_test(newT, stat="SH")
    print("pval: {}".format(pval))
    print("Dlnl: {}".format(Dlnl))

if __name__ == "__main__":
    main()
