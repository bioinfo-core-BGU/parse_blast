#from distutils.core import setup
from setuptools import find_packages, setup

setup(
    name                = 'parse-blast-general',
    version             = '0.0.1',
    author              = 'Menachem Sklarz',
    author_email        = 'sklarz@bgu.ac.il',
    maintainer          = 'Menachem Sklarz',
    maintainer_email    = 'sklarz@bgu.ac.il',
    url                 = 'https://github.com/bioinfo-core-BGU/parse_blast',
    description         = 'A script for parsing tabular BLAST output',
    license             = 'Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description    =  open('README').read(),
    download_url        = 'https://github.com/bioinfo-core-BGU/parse_blast.git',
    platforms           = ["POSIX","Windows"],
    scripts             = ['bin/parse_blast.R',
                            'bin/compare_blast_parsed_reports.R'],
    )
    

    
