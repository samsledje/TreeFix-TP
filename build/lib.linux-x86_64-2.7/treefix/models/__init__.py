#
# Python module for statistical models and reconciliation costs
#

# python libraries
import optparse, sys, os

# rasmus libraries
from rasmus import treelib, util

# compbio libraries
from compbio import phylo

class Model(object):
    def __init__(self, extra):
        self.VERSION = "-"

    def __del__(self):
        """Cleans up the model"""

    def print_help(self, file=sys.stdout):
        """Print help"""
        if self.parser:
            self.parser.print_help(file)

    def _parse_args(self, extra):
        """Parses the extra parameter arguments"""
        if self.parser is None:
            if extra is not None:
                raise Exception("--extra parameters not allowed if no parser defined")
        else:
            options, args = self.parser.parse_args(extra.split() if extra is not None else [])
            for k, v in vars(options).iteritems():
                setattr(self, k, v)


class StatModel(Model):
    def __init__(self, extra):
        """Initializes the model"""
        Model.__init__(self, extra)

        self.rooted = True
        self.parser = None

    def optimize_model(self, gtree, aln):
        """Optimizes the underlying model in the module given the tree and seq (alignment)"""
        pass

    def compute_lik_test(self, gtree, stat, alternative=None):
        """
        Computes the test statistic for tree likelihood equivalence

        Returns the p-value and Dlnl (delta lnl = best lnl - current lnl).
        Some test statistics require the user to specify the alternative hypothesis.
        """
        raise


class CostModel(Model):
    def __init__(self, extra):
        """Initializes the model"""
        Model.__init__(self, extra)

        self.mincost = -util.INF
        self.parser = None

    def optimize_model(self, gtree, stree, gene2species):
        """Optimizes the underlying model in the module given the tree"""
        self.gene2species = gene2species
        self.stree = stree

    def _reroot_helper(self, gtree, newCopy=True, returnEdge=False):
        """
        Yields rerooted trees.
        Adapted from phylo.recon_root.
        """

        # make a consistent unrooted copy of gene tree
        if newCopy:
            gtree = gtree.copy()

        if len(gtree.leaves()) == 2:
            raise StopIteration

        oldroot = gtree.root.name
        treelib.unroot(gtree, newCopy=False)
        treelib.reroot(gtree,
                       gtree.nodes[sorted(gtree.leaf_names())[0]].parent.name,
                       onBranch=False, newCopy=False)

        # make rerooting order consistent using hash ordering
        phylo.hash_order_tree(gtree, self.gene2species)

        # get list of edges to root on
        edges = []
        def walk(node):
            edges.append((node, node.parent))
            if not node.is_leaf():
                node.recurse(walk)
                edges.append((node, node.parent))
        for child in gtree.root.children:
            walk(child)

        # try initial root
        treelib.reroot(gtree, edges[0][0].name, newCopy=False)
        gtree.rename(gtree.root.name, oldroot)
        if returnEdge:
            yield gtree, edges[0]
        else:
            yield gtree
        rootedge = sorted(edges[0])

        # try rerooting on everything
        for edge in edges[1:]:
            if sorted(edge) == rootedge:
                continue
            rootedge = sorted(edge)

            node1, node2 = edge
            if node1.parent != node2:
                node1, node2 = node2, node1
            assert node1.parent == node2, "%s %s" % (node1.name, node2.name)

            # new root and cost
            treelib.reroot(gtree, node1.name, newCopy=False, keepName=True)
            if returnEdge:
                yield gtree, edge
            else:
                yield gtree

    def recon_root(self, gtree, newCopy=True, returnCost=False):
        """
        Reroots the tree by minimizing the cost function.
        Rerooted tree must keep same node names as the original tree.
        Adapted from phylo.recon_root.
        """

        # try all rerootings
        mincost = util.INF
        for gtree, edge in self._reroot_helper(gtree, newCopy=newCopy, returnEdge=True):
            cost = self.compute_cost(gtree)
            if cost < mincost:
                mincost = cost
                minroot = edge

        # root tree by minroot
        if edge != minroot:
            node1, node2 = minroot
            if node1.parent != node2:
                node1, node2 = node2, node1
            assert node1.parent == node2
            treelib.reroot(gtree, node1.name, newCopy=False, keepName=True)

        if returnCost:
            return gtree, mincost
        else:
            return gtree

    def compute_cost(self, gtree):
        """Returns the species tree aware cost."""
        raise
