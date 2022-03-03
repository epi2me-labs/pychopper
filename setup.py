#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup
from glob import glob

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = [
    'edlib',
    'parasail',
    'matplotlib',
    'seaborn',
    'tqdm',
    'six',
    'pandas',
    'pytest',
    'sphinx',
    'sphinx_rtd_theme',
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pychopper',
    version='2.7.0',
    description="A tool to identify full length cDNA reads.",
    long_description=readme,
    author="ONT Applications Group",
    author_email='Apps@nanoporetech.com',
    url='',
    packages=[
        'pychopper',
        'pychopper.phmm_data',
        'pychopper.primer_data'
    ],
    package_dir={'pychopper':
                 'pychopper'},
    package_data={'pychopper': ['primer_data/*.fas', 'phmm_data/*.*']},
    include_package_data=True,
    install_requires=requirements,
    zip_safe=False,
    keywords='pychopper',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
    ],
    tests_require=test_requirements,
    scripts=[x for x in glob('scripts/*.py') if x != 'scripts/__init__.py']
)
