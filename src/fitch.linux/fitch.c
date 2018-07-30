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
#include "bitset.h"
#include "hashmap.h"

static char * progname;
static char progheader[80];
static char * cmdline;

/* global error message buffer */
char errmsg[200] = {0};

/* number of mandatory options for the user to input */
static const char mandatory_options_count = 1;
static const char * mandatory_options_list = " --tree_file";

/* options */
long opt_help;
long opt_version;
int opt_quiet;
char * opt_treefile;
long opt_treeshow;
long opt_info;
long opt_fitch;


static struct option long_options[] =
{
  {"help",                 no_argument,       0, 0 },  /*  0 */
  {"version",              no_argument,       0, 0 },  /*  1 */
  {"quiet",                no_argument,       0, 0 },  /*  2 */
  {"tree_file",            required_argument, 0, 0 },  /*  3 */
  {"tree_show",            no_argument,       0, 0 },  /*  4 */
  {"info",                 no_argument,       0, 0 },  /*  5 */
  {"fitch",                no_argument,       0, 0 },  /*  6 */
  { 0, 0, 0, 0 }
};

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  int mand_options = 0;

  /* set defaults */

  progname = argv[0];

  opt_help = 0;
  opt_version = 0;
  opt_quiet = 1;
  opt_treefile = NULL;
  opt_treeshow = 0;
  opt_info = 0;
  opt_fitch = 0;
  

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        opt_quiet = 1;
        break;

      case 3:
        free(opt_treefile);
        opt_treefile = optarg;
        break;

      case 4:
        opt_treeshow = 1;
        break;

      case 5:
        opt_info = 1;
        break;

      case 6:
          opt_fitch = 1;
          break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for mandatory options */
  if (opt_treefile)
  {
    mand_options++;
  }

  /* check for number of independent commands selected */
  if (opt_version)
  {
    mand_options = mandatory_options_count;
    commands++;
  }
  if (opt_help)
  {
    mand_options = mandatory_options_count;
    commands++;
  if (opt_treeshow)
    commands++;
  if (opt_info)
    commands++;
  if (opt_fitch)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  /* if no command specified, turn on --help */
  if (!commands)
  {
    opt_help = 1;
    return;
  }
  
  /* check for mandatory options */
  if (mand_options != mandatory_options_count)
    fatal("Mandatory options are:\n\n%s", mandatory_options_list);

}

void cmd_help()
{
  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help                           Display help information.\n"
          "  --version                        Display version information.\n"
          "  --quiet                          Only output warnings and fatal errors to stderr.\n"
          "  --fitch                          Calculate number of transmissions using Fitch's algorithm.\n"
          "  --tree_file FILENAME             Tree file in newick format. This is a mandatory argument.\n"
         );
}

int args_getdouble2(char * arg, double * a, double * b)
{
  int len;

  int ret = sscanf(arg, "%lf,%lf%n", a, b, &len);
  
  if ((ret != 2) || ((unsigned int)(len)) < strlen(arg))
    return 0;

  return 1;
}

void cmd_tree_show()
{
  FILE * out;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");

  rtree_t * rtree = rtree_parse_newick(opt_treefile);

    rtree_show_ascii(out,rtree);
    rtree_destroy(rtree);
  }

}

void fitch_prepare_leaf (rtree_t * node) {
    char * label = node->label;
    node->host = strtok(label, "_");
    node->strain = strtok(NULL, "");
}

int fitch_recurse_up (rtree_t * rtree, int score, int set_size) {

  // Base Case
  if (rtree->leaves == 1) return 0;

  // Recursive Step
  int s = fitch_recurse_up(rtree->left, score, set_size) 
            + fitch_recurse_up(rtree->right, score, set_size);
  rtree->host_set = bitset_create(set_size);
  bitset_intersection(rtree->host_set,
                      rtree->left->host_set,
                      rtree->right->host_set);
  rtree->host = "intersection";
  if (bitset_is_empty(rtree->host_set)) {
    bitset_union( rtree->host_set,
                  rtree->left->host_set,
                  rtree->right->host_set);
    rtree->host = "union";
    s++;
  };
  return s;
}

void fitch_recurse_down (rtree_t * rtree) {
  if (rtree->leaves == 1) return;
  if (bitset_size(rtree->host_set) == 1) {
    //rtree->host = "one element";
  }
  else {
    //rtree->host = "multple";
  }
  fitch_recurse_down(rtree->left);
  fitch_recurse_down(rtree->right);
}

void cmd_fitch() {

  rtree_t * rtree = rtree_parse_newick(opt_treefile);
  
  //Set host and strain for each leaf
  //create host string --> int mappings
  int n_tips = rtree->leaves;
  rtree_t ** node_list = (rtree_t **)calloc(rtree->leaves,sizeof(rtree_t *)); 
  rtree_query_tipnodes(rtree, node_list);

  map_t * host_map = hashmap_new();
  
  for (int i = 0; i < n_tips; i++) {
    rtree_t * node = node_list[i];
    fitch_prepare_leaf(node);
    hashmap_put(host_map, node->host);
  }
  
  // Number of hosts, size of bitset will be the first multiple of (8 * sizeof(int)) greater than the number of hosts.
  //int set_size = bst_size(host_map);
  int set_size = hashmap_length(host_map);

  //Initialize leaf host sets

  for (int i = 0; i < n_tips; i++) {
    rtree_t * node = node_list[i];
    node->host_set = bitset_create(set_size);
    bitset_insert(node->host_set, hashmap_get(host_map, node->host));
  }

  //Upward step of Fitch's algorithm
  int score = fitch_recurse_up(rtree, 0, set_size);
  printf("Score: %d\n", score);
  
  //Downward step of Fitch's algorithm
  //fitch_recurse_down(rtree);

  // Show final tree
  // FILE * out;
  // out = opt_outfile ? xopen(opt_outfile,"w") : stdout;
  // rtree_show_ascii(out,rtree);

  //Deallocate memory
  rtree_destroy(rtree);
  hashmap_free(host_map);
  free(node_list);
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fill_header()
{
  snprintf(progheader, 80,
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           sysconf(_SC_NPROCESSORS_ONLN));
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "Modified for UConn Computational Biology Lab\n");
  fprintf(stdout, "https://github.com/xflouris/newick-tools\n");
  fprintf(stdout,"\n");

}

int main (int argc, char * argv[])
{
  if (argc < 2) {
    printf("Please pass the path to a tree file on the command line.\n");
    exit(1);
  }
  opt_treefile = argv[1];
  cmd_fitch();
}
