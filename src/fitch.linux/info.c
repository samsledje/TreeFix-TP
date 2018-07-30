/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "newick-tools.h"

static void show_tree_info(int tip_count,
                           int inner_count,
                           int min_inner_degree,
                           int max_inner_degree)
{
  printf("Leaves (tip nodes): %d\n"
         "Inner nodes: %d\n"
         "Total nodes: %d\n"
         "Edges: %d\n"
         "Minimum inner node degree: %d\n"
         "Maximum inner node degree: %d\n",
         tip_count,
         inner_count, 
         tip_count + inner_count,
         tip_count + inner_count - 1,
         min_inner_degree,
         max_inner_degree);
}

static void rtree_info(rtree_t * root)
{
  int tip_count = root->leaves;
  int inner_count = tip_count - 1;
  int max_inner_degree = 3;
  int min_inner_degree = 2;

  show_tree_info(tip_count,
                 inner_count,
                 min_inner_degree,
                 max_inner_degree);

  double * outbuffer = (double *)xmalloc((tip_count+inner_count-1) * 
                                          sizeof(double));

  int count = rtree_query_branch_lengths(root, outbuffer);

  double min,max,mean,median,var,stdev;
  stats(outbuffer,count,&min,&max,&mean,&median,&var,&stdev);
  printf("Min. branch length: %f\n", min);
  printf("Max. branch length: %f\n", max);
  printf("Mean branch length: %f\n", mean);
  printf("Median branch length: %f\n", median);
  printf("Variance branch length: %f\n", var);
  printf("Standard deviation branch length: %f\n", stdev);
  printf("Longest lineage: %f\n", rtree_longest_path(root));

  free(outbuffer);
}

void cmd_info(void)
{
  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  rtree_t * rtree = rtree_parse_newick(opt_treefile);

  if (rtree)
  {
    if (!opt_quiet)
      printf("Loaded binary rooted tree...\n");

    /* show info */
    rtree_info(rtree);

    /* deallocate tree structure */
    rtree_destroy(rtree);
  }
}

