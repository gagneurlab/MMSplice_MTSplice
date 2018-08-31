#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ 
    'keras',
    'kipoi',
    'numpy',
    'pandas',
    'concise',
    'cyvcf2',
    'gffutils',
    'pyfaidx',
    'tqdm'
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Jun Cheng",
    author_email='chengju@in.tum.de',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        # "Programming Language :: Python :: 2",
        # 'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Predict splicing variant effect from VCF",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    # package_data={
    #     'eis.models':
    #     ['Acceptor.h5',
    #     'Donor.h5',
    #     'Exon.h5',
    #     'Intron3.h5',
    #     'Intron5.h5']
    # }, # Done with MANIFEST.in
    include_package_data=True,
    keywords='eis',
    name='eis',
    packages=find_packages(include=['eis']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/s6juncheng/eis',
    version='0.2.0',
    zip_safe=False,
)
