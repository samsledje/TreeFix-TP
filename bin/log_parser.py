import sys
import re
import numpy as np
import pandas as pd
import dendropy as dp
from dendropy.calculate import treecompare
import optparse

NAME = 'TreeFix Log Parser'
AUTHOR = 'Samuel Sledzieski'
VERSION = '1.0'

parser = optparse.OptionParser(usage = 'python treefix_log_parser.py [options] <log file>')
parser.add_option('--out', dest='out_path', default=None, help='Path for output')
parser.add_option('--near', dest='near_percent', default=20, help='Trees within <--near>% of the optimal cost will be captured')
parser.add_option('--true', dest='true_tree_path', default=None, help='Can provide a true tree to compare multiple optimal trees with')
parser.add_option('--include_near', action='store_true', dest='include_near', default=False, help='Include nearby optimal trees in summary statistics')
parser.add_option('--separate_trees', action='store_true', dest='separate_trees', default=False, help='Create two output files separating trees and statistics')

options, args = parser.parse_args()
if len(args) == 0 or len(args) > 1:
    parser.print_help()
    sys.exit(1)

print('### {} Version {} ###'.format(NAME, VERSION))

file_path = args[0]
out_path = options.out_path
near_percent = options.near_percent
true_tree_path = options.true_tree_path
include_near = options.include_near
separate_trees = options.separate_trees

print('Logfile: {}'.format(file_path))
if true_tree_path is not None:
    try:
        true_tree = dp.Tree.get(path=true_tree_path, schema='newick')
    except:
        print('True tree path is not a valid tree file')
        sys.exit(1)
else:
    true_tree = False


if out_path:
    f = open(out_path, 'w+')
else:
    f = sys.stdout

if separate_trees:
    if not out_path:
        print('Cannot separate trees if using stdout')
        sys.exit(1)
    g = open('{}.trees'.format(out_path), 'w+')
else:
    g = f

final = False
current_iter = 0
best_cost = float("inf")
accepted_iters = []
best_iters = []
best_trees = []
near_iters = []
near_trees = []

print('Parsing log file...')
with open(file_path, 'r') as h:
    for line in h:
        if line.strip().startswith('search: cost'):
            current_cost = int(re.search(r'\d+$', line.strip()).group())
            if current_cost < best_cost:
                best_cost = current_cost
near_cost = float(best_cost * (1 + (float(near_percent)/100)))

with open(file_path, 'r') as h:
    for line in h:
        line = line.strip()
        if final and line.startswith('search: changed') and line.endswith('no'):
            break
        if line.startswith('search: iter'):
            current_iter = int(re.search(r'\d+$', line).group())
        elif line.startswith('search: final') and not line.startswith('search: final cost'):
            current_iter = 'Final'
            final = True
        elif line.startswith('search: cost') or line.startswith('search: final cost'):
            current_cost = int(re.search(r'\d+$', line).group())
            if current_cost == best_cost:
                best_iters.append(current_iter)
            elif current_cost <= near_cost:
                near_iters.append(current_iter)
        elif current_iter in best_iters and line.startswith('tree:') and line.endswith(';'):
                tree_string = re.search('tree: (.*)', line).group(1)
                new_tree = tree_string
                best_trees.append(new_tree)
                if final:
                    break
        elif current_iter in near_iters and line.startswith('tree:') and line.endswith(';'):
                tree_string = re.search('tree: (.*)', line).group(1)
                new_tree = tree_string
                near_trees.append(new_tree)

assert len(best_iters) == len(best_trees)

# Best and Near Trees
g.write('#### Best Trees ####\n')
for i in range(len(best_iters)):
    g.write('>{}\n'.format(str(best_iters[i])))
    g.write('{}\n'.format(best_trees[i]))

g.write('\n#### Near Best Trees ####\n')
for i in range(len(near_iters)):
    g.write('>{}\n'.format(str(near_iters[i])))
    g.write('{}\n'.format(near_trees[i]))

if include_near:
    combined_iters = best_iters + near_iters
    combined_trees = best_trees + near_trees
else:
    combined_iters = best_iters
    combined_trees = best_trees

N = len(combined_trees)

f.write('\n#### Summary ####\n')
f.write('There are {} trees with cost {}\n'.format(len(best_iters), best_cost))
f.write('There are {} more trees within {}% of the best cost ({} < cost <= {})\n'.format(len(near_iters), near_percent, best_cost, near_cost))
test_t = dp.Tree.get_from_string(combined_trees[0], 'newick')
f.write('Number of taxa is {}\n'.format(len(test_t.leaf_nodes())))

# Pairwise Robinson Foulds
print('Calculating pairwise Robinson-Foulds distances...')
f.write('\n#### Pairwise Robinson-Foulds Distances ####\n')
dp_trees = [dp.Tree.get_from_string(i, 'newick') for i in combined_trees]
tlist = dp.TreeList(dp_trees)
rf_pair_matrix = np.zeros((N,N))
for i in range(N):
    for j in range(N):
        if i != j:
            rf_pair_matrix[i,j] = treecompare.symmetric_difference(dp_trees[i], dp_trees[j])
pd.set_option('display.max_columns', None)
pd_pair_matrix = pd.DataFrame(data=rf_pair_matrix,index=combined_iters, columns=combined_iters)
f.write('{}\n'.format(str(pd_pair_matrix)))
if len(combined_iters) > 1:
    f.write('\nMinimum (non-zero) RF Distance: {}\n'.format(rf_pair_matrix[np.nonzero(rf_pair_matrix)].min()))
    f.write('Maximum RF Distance: {}\n'.format(rf_pair_matrix.max()))
    f.write('Mean RF Distance: {}\n'.format(rf_pair_matrix[np.nonzero(rf_pair_matrix)].mean()))
    f.write('Std. Dev. of RF Distance: {}\n'.format(np.sqrt(rf_pair_matrix[np.nonzero(rf_pair_matrix)].var())))

# Consensus Tree Methods
print('Calculating consensus trees...')
f.write('\n#### Strict Consensus Tree ####\n')
strict_con_tree = tlist.consensus(min_freq=1.0)
f.write('{}\n'.format(strict_con_tree.as_string('newick')))
strict_stats = np.zeros((N,3))
for i in range(N):
    fp, fn = treecompare.false_positives_and_negatives(strict_con_tree,dp_trees[i])
    strict_stats[i] = [fp, fn, fp+fn]
pd_strict_stats = pd.DataFrame(data=strict_stats,index=combined_iters,columns=['False Positive', 'False Negative', 'RF Distance'])
f.write('{}\n'.format(str(pd_strict_stats)))
f.write('\nMinimum RF Distance: {}\n'.format(strict_stats[:,2].min()))
f.write('Maximum RF Distance: {}\n'.format(strict_stats[:,2].max()))
f.write('Mean RF Distance: {}\n'.format(strict_stats[:,2].mean()))
f.write('Std. Dev. of RF Distance: {}\n'.format(np.sqrt(strict_stats[:,2].var())))

f.write('\n#### Majority Rule Consensus Tree ####\n')
maj_con_tree = tlist.consensus(min_freq=0.5)
f.write('{}\n'.format(maj_con_tree.as_string('newick')))
maj_stats = np.zeros((N,3))
for i in range(N):
    fp, fn = treecompare.false_positives_and_negatives(maj_con_tree, dp_trees[i])
    maj_stats[i] = [fp, fn, fp+fn]
pd_maj_stats = pd.DataFrame(data=maj_stats,index=combined_iters,columns=['False Positive', 'False Negative', 'RF Distance'])
f.write('{}\n'.format(str(pd_maj_stats)))
f.write('\nMinimum RF Distance: {}\n'.format(maj_stats[:,2].min()))
f.write('Maximum RF Distance: {}\n'.format(maj_stats[:,2].max()))
f.write('Mean RF Distance: {}\n'.format(maj_stats[:,2].mean()))
f.write('Std. Dev. of RF Distance: {}\n'.format(np.sqrt(maj_stats[:,2].var())))

# Comparison with True Tree
if true_tree:
    print('Calculating Robinson-Foulds distances to the true tree...')
    f.write('\n#### Comparison to True Tree ####\n')
    true_matrix = np.zeros((N+2,3))
    true_plus_all = dp_trees + [strict_con_tree,maj_con_tree,true_tree]
    truelist = dp.TreeList(true_plus_all)
    for i in range(N+2):
        fp, fn = treecompare.false_positives_and_negatives(true_tree, true_plus_all[i])
        true_matrix[i] = [fp,fn,fp+fn]
    combined_consensus_iters = combined_iters+['Strict Consensus', 'Majority Consensus']
    combined_consensus_trees = combined_trees+[strict_con_tree, maj_con_tree]
    pd_true_matrix = pd.DataFrame(data=true_matrix,index=combined_consensus_iters,columns=['False Positive', 'False Negative', 'RF Distance'])
    f.write('{}\n'.format(str(pd_true_matrix)))
    f.write('\nMinimum RF Distance: {}\n'.format(true_matrix[:,2].min()))
    f.write('Maximum RF Distance: {}\n'.format(true_matrix[:,2].max()))
    f.write('Mean RF Distance: {}\n'.format(true_matrix[:,2].mean()))
    f.write('Std. Dev. of RF Distance: {}\n'.format(np.sqrt(true_matrix[:,2].var())))
    f.write('Tree {} is closest to the true tree\n'.format(combined_consensus_iters[np.argmin(true_matrix[:,2])]))
    f.write('{}\n'.format(combined_consensus_trees[np.argmin(true_matrix[:,2])]))

f.close()
if out_path and not separate_trees:
    print('Optimal and Near-Optimal Trees and Summary Statistics written to {}'.format(out_path))
elif out_path and separate_trees:
    print('Optimal and Near-Optimal Trees written to {}.trees'.format(out_path))
    print('Summary Statistics written to {}'.format(out_path))
    g.close()
