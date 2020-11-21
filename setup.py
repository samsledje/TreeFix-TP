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
    long_description = """
            """,
    author='Mukul Bansal, Samuel Sledzieski, Yi-Chieh Wu',
    author_email='samuel.sledzieski@uconn.edu',
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
            ],
        },
    ext_modules=[raxml_module]
    )


