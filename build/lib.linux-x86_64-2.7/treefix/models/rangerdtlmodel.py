#
# Python module for ranger-dtl cost
#

# treefix libraries
from treefix.models import CostModel

# python libraries
import optparse
import os, sys, subprocess
import re
import tempfile

# rasmus libraries
from rasmus import treelib, util

# compbio libraries
from compbio import phylo

#=============================================================================
# command

# uses ranger-dtl-U v1.0
##cmd = os.path.join(os.path.realpath(os.path.dirname(__file__)),
##                   "ranger-dtl-U.linux")
##cmd = "ranger-dtl-U.linux"

patt = "The minimum reconciliation cost is: (?P<cost>\d+) \(Duplications: (?P<D>\d+), Transfers: (?P<T>\d+), Losses: (?P<L>\d+)\)"

#============================================================================

class DTLModel(CostModel):
    """Computes DTL costs"""

    def __init__(self, extra):
        """Initializes the model"""
        CostModel.__init__(self, extra)

        self.VERSION = "0.1.2"
        self.mincost = 0

        parser = optparse.OptionParser(prog="DTLModel")
        parser.add_option("--cmd", dest="cmd",
                          metavar="<ranger-dtl-U command>",
                          default="ranger-dtl-U",
                          help="ranger-dtl-U command (default: \"ranger-dtl-U\")")
        parser.add_option("-D", "--dupcost", dest="dupcost",
                          metavar="<duplication cost>",
                          default=2, type="int",
                          help="duplication cost, integer only (default: 2)")
        parser.add_option("-T", "--transcost", "--transfercost", dest="transcost",
                          metavar="<transfer cost>",
                          default=3, type="int",
                          help="transfer cost, integer only (default: 3)")
        parser.add_option("-L", "--losscost", dest="losscost",
                          metavar="<loss cost>",
                          default=1, type="int",
                          help="loss cost, integer only (default: 1)")
        parser.add_option("--seed", dest="seed",
                          metavar="<seed>",
                          type="int",
                          help="user defined random number generator seed")
        parser.add_option("--tmp", dest="tmp",
                          metavar="<tmp directory>",
                          help="directory for temporary files (must exist)")

        self.parser = parser

        CostModel._parse_args(self, extra)

        # check temporary directory
        if self.tmp:
            print self.temp
            print "\n\n\n"
            if not os.path.exists(os.path.realpath(self.tmp)):
                raise Exception("--tmp directory does not exist")

        # make temporary file
        fd, self.treefile = tempfile.mkstemp(dir=self.tmp)
        os.close(fd)

        # hack for cygwin (ranger-dtl-U cannot handle system files)
        if sys.platform == 'cygwin':
            cwd = os.getcwd()
            if self.treefile.startswith(cwd):
                # remove working path (and backslash)
                self.treefile = self.treefile[len(cwd)+1:]
            else:
                raise Exception("--tmp must be a relative path when using cygwin")

    def __del__(self):
        """Cleans up the model"""
        # delete temporary file
        os.remove(self.treefile)

    def optimize_model(self, gtree, stree, gene2species):
        """Optimizes the model"""
        CostModel.optimize_model(self, gtree, stree, gene2species)

        if self.dupcost < 0:
            self.parser.error("-D/--dupcost must be >= 0")
        if self.transcost < 0:
            self.parser.error("-T/--transcost must be >= 0")
        if self.losscost < 0:
            self.parser.error("-L/--losscost must be >= 0")

        # ensure gtree and stree are both rooted and binary
        if not (treelib.is_rooted(gtree) and treelib.is_binary(gtree)):
            raise Exception("gene tree must be rooted and binary")
        if not (treelib.is_rooted(stree) and treelib.is_binary(stree)):
            raise Exception("species tree must be rooted and binary")
        try:
            junk = phylo.reconcile(gtree, stree, gene2species)
        except:
            raise Exception("problem mapping gene tree to species tree")

    def recon_root(self, gtree, newCopy=True, returnCost=False):
        """
        Returns input gene tree and min DTL cost.

        Note: The optimally rooted gene tree from ranger-dtl-U
              contains species names, so this function does NOT reroot.
              It simply delegates to compute_cost.
        """
        if newCopy:
            gtree = gtree.copy()

        if returnCost:
            mincost = self.compute_cost(gtree)
            return gtree, mincost
        else:
            return gtree

    def compute_cost(self, gtree):
        """Returns the DTL cost"""

        # write species tree and gene tree using species map
        treeout = util.open_stream(self.treefile, 'w')
        self.stree.write(treeout, oneline=True, writeData=lambda x: "")
        treeout.write("\n[&U]")
        gtree.write(treeout, namefunc=lambda name: self.gene2species(name),
                    oneline=True, writeData=lambda x: "")
        treeout.write("\n")
        treeout.close()

        # create command
        args = [self.cmd,
                '-i', self.treefile,
                '-D', str(self.dupcost),
                '-T', str(self.transcost),
                '-L', str(self.losscost)]
        if self.seed:
            args.extend(['--seed', str(self.seed)])

        # execute command
        try:
            proc = subprocess.Popen(args,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT,
                                    universal_newlines=True)
            ret = proc.wait()
        except:
            raise Exception("ranger-dtl-U failed")
        if ret != 0:
            raise Exception("ranger-dtl-U failed with returncode %d" % ret)

        # parse output
        cost = None
        for line in proc.stdout:
            m = re.match(patt, line)
            if m:
                cost = int(m.group("cost"))
                break
        assert cost is not None

        return cost
