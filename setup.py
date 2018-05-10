#!/usr/bin/env python
"""
MultiQC_bcbio is a plugin for MultiQC, providing additional tools which are
specific to bcbio_nextgen.py pipeline.

For more information about bcbio_nextgen, see http://github.com/chapmanb/bcbio_nextgen
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '0.2.6'

setup(
    name = 'multiqc_bcbio',
    version = version,
    author = 'Lorena Pantano, Brad Chapman, Vlad Saveliev, Phil Ewels',
    author_email = 'lpantano@iscb.org',
    description = "MultiQC plugin for bcbio_nextgen pipeline",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/MultiQC/MultiQC_bcbio',
    download_url = 'https://github.com/MultiQC/MultiQC_bcbio/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    zip_safe=False,
    install_requires = [
        'multiqc>=1.0',
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'bcbio = multiqc_bcbio.bcbio:MultiqcModule',
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_bcbio:multiqc_bcbio_config',
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)

