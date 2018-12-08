#!/usr/bin/env python
#
# setup for TreeFix-VP library packages
#
# use the following to install:
#   python setup.py build
#   python setup.py install
#

import os, sys, shutil
from setuptools import setup, Extension

# version control
sys.path.insert(0, os.path.realpath(
            os.path.join(os.path.dirname(__file__), "python")))
from treefix import treefixVP
VERSION = treefixVP.PROGRAM_VERSION_TEXT

# find correct ranger-dtl-U executable
if not os.path.exists('bin/fitch.linux'):
    if sys.platform == 'darwin':
        ranger_dtl_script = 'bin/fitch.linux'
    elif sys.platform == 'cygwin':  #windows
        ranger_dtl_script = 'bin/fitch.linux'
    else:
        ranger_dtl_script = 'bin/fitch.linux'
    shutil.copy(ranger_dtl_script, 'bin/fitch.linux')

# platform dependency
extra_link_args = []
if sys.platform != 'darwin':
    extra_link_args.append('-s')

# raxml sources
srcs = [os.path.join('src/raxml',fn) for fn in os.listdir('src/raxml')
        if (not os.path.isdir(fn)) and fn.endswith('.c')]
raxml_module = Extension('treefix_raxml._raxml',
                         sources=['python/treefix_raxml/raxml.i'] + srcs,
                         extra_link_args=extra_link_args
                         )
                         
setup(
    name='treefixVP',
    version=VERSION,
    description='TreeFix-VP',
    long_description = """
            """,
    author='Mukul Bansal, Samuel Sledzieski, Yi-Chieh Wu',
    author_email='samuel.sledzieski@uconn.edu',
    url='https://www.github.com/samsledje/treefixVP',

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

    package_dir = {'': 'python'},
    packages=['treefix',
              'treefix.models',
              'treefix.deps',
              'treefix.deps.rasmus',
              'treefix.deps.compbio',
	      'treefix_raxml',
	      'treefix_raxml.deps.rasmus',
	      'treefix_raxml.deps.compbio'],
    py_modules=[],
    scripts=['bin/treefix_for_VP', 'bin/treefixVP', 'bin/fitch.linux', 'bin/log_parser'],
    ext_modules=[raxml_module]
    )


