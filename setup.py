#!/usr/bin/env python
#
# setup for TreeFix-TP library packages
#
# use the following to install:
#   python setup.py build
#   python setup.py install
#

import os, sys, shutil
from setuptools import setup, Extension, find_packages

# version control
import treefix_tp
VERSION = treefix_tp.PROGRAM_VERSION_TEXT

# raxml sources
srcs = [os.path.join('treefix_tp/pyRAxML/src',fn) for fn in os.listdir('treefix_tp/pyRAxML/src')
        if (not os.path.isdir(fn)) and fn.endswith('.c')]
raxml_module = Extension('treefix_tp.pyRAxML._raxml',
                         sources=['treefix_tp/pyRAxML/raxml.i',
                                  'treefix_tp/pyRAxML/tmaps.i'] + srcs,
                         )

setup(
    name='TreeFix-TP',
    version=VERSION,
    description='TreeFix-TP',
    long_description = """TreeFix-TP is a program for reconstructing highly accurate transmission phylogenies, i.e., phylogenies depicting the evolutionary relationships between infectious disease strains (viral or bacterial) transmitted between different hosts. TreeFix-TP is designed for scenarios where multiple strain sequences have been sampled from each infected host, and it uses the host assignment of each sequence sample to error-correct a given maximum likelihood phylogeny of the strain sequences. Specifically, given a maximum likelihood phylogeny, the multiple sequence alignment on which the phylogeny was built, and the host assignment for each sequence, TreeFix-TP searches around the maximum likelihood phylogeny to find an alternate error-corrected phylogeny which is equally well-supported by the sequence data and minimizes the number of necessary inter-host transmissions.""",
    author='Samuel Sledzieski, Mukul Bansal',
    author_email='mukul.bansal@uconn.edu',
    url='https://www.github.com/samsledje/TreeFix-TP',

    classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Operating System :: POSIX',
          'Programming Language :: Python',
          'Topic :: Education',
          ],
    packages = find_packages(),
    py_modules=[],
    data_files=[('bin', ['bin/fitch.linux', 'bin/log_parser.py'])],
    entry_points={
        "console_scripts": [
            "treefix-tp = treefix_tp.__main__:main",
            "ttp-parse-log = treefix_tp.scripts.ttp_parse_log:main",
            "ttp-check-cost = treefix_tp.scripts.ttp_check_cost:main",
            "ttp-check-likelihood = treefix_tp.scripts.ttp_check_likelihood:main"
            ],
        },
    ext_modules=[raxml_module]
    )


