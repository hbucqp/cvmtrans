#!/usr/bin/env python3
#-*- coding:utf-8 -*-

#############################################
# File Name: setup.py
# Author: Jake Cui
# Mail: cqp@cau.edu.cn
# Created Time:  2024-04-28 22:04:41
#############################################

import os
import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup


requirements = [
    'Bio', 'pandas', 'setuptools', 'numpy', 'dnaio', 'pysam'
]


about = {}
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'cvmtrans', '__init__.py'), 'r', encoding='utf-8') as f:
    exec(f.read(), about)


# Get the long description from the relevant file
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="cvmtrans",
    version=about['__version__'],
    keywords=["pip", "transposon", "transposon sequecing"],
    description="Transposon data process",
    long_description=long_description,
    long_description_content_type='text/markdown',
    license="MIT Licence",
    url=about['__url__'],
    author=about['__author__'],
    author_email=about['__author_email__'],
    packages=find_packages(),
    include_package_data=True,
    # package_data={'': ['*']},
    platforms="any",
    install_requires=requirements,
    classifiers=[
        # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Development Status :: 3 - Alpha',
        # Define that your audience are developers
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',   # Again, pick a license
        # Specify which pyhton versions that you want to support
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11'
    ],
    entry_points={
        'console_scripts': [
            'cvmtrans=cvmtrans.cvmtrans:main',
            'cvmtrans_cut_tags=cvmtrans.cvmtrans_cut_tags:main',
            'cvmtrans_extract_reads=cvmtrans.cvmtrans_extract_reads:main',
        ],
    },
)
