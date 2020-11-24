#
# Python module for fitch cost
#

# python libraries
import optparse
import os, sys, subprocess
import re
import tempfile

# treefix libraries
from treefix_tp.models import CostModel

# rasmus libraries
from treefix_tp.rasmus import treelib, util

# compbio libraries
from treefix_tp.compbio import phylo

# =============================================================================
# command

# uses fitch.linux
CMD = os.path.join(os.path.realpath(os.path.dirname(__file__)),
                   "../../bin/fitch.linux")
#CMD = "fitch.linux"

patt = r"Score: (?P<cost>\d+)"

# ============================================================================


class FitchModel(CostModel):
    """Computes cost of running Fitch's algorithm for a given tree"""

    def __init__(self, extra=None):
        """Initializes the model"""
        CostModel.__init__(self, extra)

        self.VERSION = "1.0"
        self.mincost = 0

        parser = optparse.OptionParser(prog="FitchModel")
        parser.add_option(
            "--cmd",
            dest="cmd",
            metavar="<fitch command>",
            default=CMD,
            help='fitch command (default: "fitch")',
        )
        parser.add_option(
            "--seed",
            dest="seed",
            metavar="<seed>",
            type="int",
            help="user defined random number generator seed",
        )
        parser.add_option(
            "--tmp",
            dest="tmp",
            metavar="<tmp directory>",
            help="directory for temporary files (must exist)",
        )

        self.parser = parser

        CostModel._parse_args(self, extra)

        # check temporary directory
        if self.tmp:
            if not os.path.exists(os.path.realpath(os.path.abspath(self.tmp))):
                raise Exception("--tmp directory does not exist")

        # make temporary file
        fd, self.treefile = tempfile.mkstemp(dir=self.tmp)
        os.close(fd)

        # hack for cygwin (ranger-dtl-U cannot handle system files)
        if sys.platform == "cygwin":
            cwd = os.getcwd()
            if self.treefile.startswith(cwd):
                # remove working path (and backslash)
                self.treefile = self.treefile[len(cwd) + 1 :]
            else:
                raise Exception("--tmp must be a relative path when using cygwin")

    def __del__(self):
        """Cleans up the model"""
        # delete temporary file
        os.remove(self.treefile)

    def optimize_model(self, gtree):
        """Optimizes the model"""

        # ensure gtree and stree are both rooted and binary
        if not (treelib.is_rooted(gtree) and treelib.is_binary(gtree)):
            raise Exception("gene tree must be rooted and binary")

    def recon_root(self, gtree, newCopy=True, returnCost=False):
        """
        Returns input gene tree and min transmission cost.

        """
        if newCopy:
            gtree = gtree.copy()

        if returnCost:
            mincost = self.compute_cost(gtree)
            return gtree, mincost
        else:
            return gtree

    def compute_cost(self, gtree):
        """Returns the cost (minimal number of transmissions)"""

        # write  gene tree using species map
        treeout = util.open_stream(self.treefile, "w")
        gtree.write(treeout, oneline=True, writeData=lambda x: "")
        # gtree.write(treeout, namefunc=lambda name: self.gene2species(name),
        # oneline=True, writeData=lambda x: "")
        treeout.write("\n")
        treeout.close()

        # create command
        args = [self.cmd, self.treefile]
        # if self.seed:
        #     args.extend(['--seed', str(self.seed)])

        # execute command
        try:
            proc = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True
            )
            ret = proc.wait()
        except:
            raise Exception("fitch.linux failed")
        if ret != 0:
            raise Exception("fitch.linux failed with returncode %d" % ret)

        # parse output
        cost = None
        for line in proc.stdout:
            m = re.match(patt, line)
            if m:
                cost = int(m.group("cost"))
                break
        assert cost is not None

        return cost
