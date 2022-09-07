# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 07:49:16 2017

@author: Philipp Ehses
"""


from setuptools import setup


setup(name='twixtools',
      version="1.0",
      description="file reader for Siemens twix(.dat)-files",
      long_description=open('README.md').read(),
      keywords='twix,siemens,mri',
      author='Philipp Ehses',
      author_email='philipp.ehses@dzne.de',
      url='https://github.com/pehses/twixtools/',
      license='MIT License',
      packages=['twixtools','twixtools/contrib'],
      scripts=['twixtools/twixzip.py','utils/convert_to_cfl.py', 'utils/rotate_3Dcfl.py'],
      zip_safe=False,
      test_suite="test",
      )
