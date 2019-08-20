# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 07:49:16 2017

@author: epracht
"""

import os
import subprocess
import re
from setuptools import setup, find_packages

VERSION_PY = """
# This file is originally generated from Git information by running 'setup.py
# install'. You won't find it in your repository, just in the installed package.
__version__ = '%s'
"""


# def update_version_py():
#     if not os.path.isdir(".git"):
#         print("This does not appear to be a Git repository.")
#         return
#     try:
#         p = subprocess.Popen(["git", "describe", "--tags", "--dirty", "--always"], stdout=subprocess.PIPE)
#
#     except EnvironmentError:
#         print("unable to run git, leaving twixreader/version.py alone")
#         return
#
#     stdout = p.communicate()[0].rstrip()
#
#     if p.returncode != 0:
#         print("unable to run git, twixreader/version.py alone")
#         return
#
#     ver = stdout
#     f = open("twixreader/version.py", "w")
#     f.write(VERSION_PY % ver)
#     f.close()
#     print("setting version to '%s'", ver)
#
#
# def get_version():
#     update_version_py()
#     try:
#         f = open("twixreader/version.py")
#     except EnvironmentError:
#         return None
#     for line in f.readlines():
#         mo = re.match("__version__ = '([^']+)'", line)
#         if mo:
#             ver = mo.group(1)
#             return ver
#     return None


setup(name='twixtools',
      version="1.0",  # get_version(),
      description="file reader for Siemens twix(.dat)-files",
      long_description=open('README.md').read(),
      keywords='twix,siemens,mri',
      author='Philipp Ehses',
      author_email='philipp.ehses@dzne.de',
      url='',
      license='',
      packages=['twixtools'],
      zip_safe=False,
      )
