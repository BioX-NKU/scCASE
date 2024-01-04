#!/usr/bin/env python
#-*- coding:utf-8 -*-

from setuptools import setup, find_packages

setup(name='scCASE',
    version='0.2.1',
    keywords=("pip", "scCASE", "single-cell"),
    url="https://hub.njuu.cf/BioX-NKU/",
    author="BioX-NKU",
    packages=find_packages(),

    python_requires='>=3.8.13',

    scripts=['scCASE.py'],
    license='MIT Licence',
    install_requires=[
        'pandas',
        'scanpy',
        'numpy',
        'anndata',
        'scikit-learn',
        'louvain',
        'leidenalg',
        'igraph'
        'episcanpy',
        'kneed>=0.7.0',],
    classifiers=['Intended Audience :: Science/Research',
      'License :: OSI Approved :: MIT License',
      'Programming Language :: Python :: 3.8',
      'Operating System :: MacOS :: MacOS X',
      'Operating System :: Microsoft :: Windows',
      'Operating System :: POSIX :: Linux',
      'Topic :: Scientific/Engineering :: Bio-Informatics']
     )




