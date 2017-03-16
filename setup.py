import os
from setuptools import setup

setup(
    name = "kvasirHGT",
    version = "v0.6.4",
    author = "Kevin Bonham, PhD",
    author_email = "kevbonham@gmail.com",
    description = "A package to identify HGT in bacterial genomes",
    license = "MIT",
    keywords = ["HGT", "biology", "bacteria", "genomics"],
    url = "http://github.com/kescobo/kvasir",
    download_url = 'https://github.com/kescobo/kvasir/archive/v0.6.4-beta.tar.gz',
    packages=['kvasir', 'tests'],
    scripts=[os.path.join('bin', 'kv_blast.py'),
             os.path.join('bin', 'kv_analysis.py'),
             os.path.join('bin', 'kv_import.py'),
             os.path.join('bin', 'kv_distance.py')],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
)
