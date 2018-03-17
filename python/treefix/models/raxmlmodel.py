#
# Python module for SH test using RAxML likelihoods
# Requires raxml libraries (http://compbio.mit.edu/treefix/#raxml)
#

# treefix libraries
from treefix.models import StatModel

# treefix_raxml library
import treefix_raxml as raxml

# python libraries
import os, sys
import optparse
import tempfile

# rasmus libraries
from rasmus import util

# compbio libraries
from compbio import phylip

class RAxMLModel(StatModel):
    """Computes test statistics using RAxML site-wise likelihoods"""

    def __init__(self, extra):
        """Initializes the RAxML model"""
        StatModel.__init__(self, extra)

        self.VERSION = "0.2.4"
        self._raxml = raxml.RAxML()
        self.rooted = self._raxml.rooted

        parser = optparse.OptionParser(prog="RAxMLModel")
        parser.add_option("-m", "--model", dest="model",
                          metavar="<model>",
                          default="GTRGAMMA",
                          help="model of nucleotide or amino acid substitution (default: GTRGAMMA)")
        parser.add_option("-e", "--eps", dest="eps",
                          metavar="<eps>",
                          default=2.0, type="float",
                          help="model optimization precision in log likelihood units (default 2.0)")
        self.parser = parser

        StatModel._parse_args(self, extra)

    def __del__(self):
        """Cleans up the RAxML model"""
        del self._raxml

    def optimize_model(self, gtree, aln):
        """Optimizes the RAxML model"""
        StatModel.optimize_model(self, gtree, aln)

        fd, treefile = tempfile.mkstemp('.tree')
        os.close(fd)
        gtree.write(treefile)

        fd, seqfile = tempfile.mkstemp('.align')
        os.close(fd)
        out = util.open_stream(seqfile, "w")
        phylip.write_phylip_align(out, aln, strip_names=False)
        out.close()

        self._raxml.optimize_model(treefile, seqfile,
                                   "-m %s -e %s -n test" % (self.model, self.eps))

        os.remove(treefile)
        os.remove(seqfile)

    def compute_lik_test(self, gtree, stat="SH", alternative=None):
        """Computes the test statistic 'stat' using RAxML likelihoods"""
        return self._raxml.compute_lik_test(gtree, stat, alternative)
