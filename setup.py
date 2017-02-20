import os
from setuptools import setup

# Utility function to read the README file.

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "kvasir",
    version = "0.6",
    author = "Kevin Bonham, PhD",
    author_email = "kevbonham@gmail.computation",
    description = ("A package to identify HGT in bacterial genomes"),
    license = "MIT",
    keywords = "HGT biology bacteria genomics",
    url = "http://github.com/kescobo/kvasir",
    packages=['kvasir', 'tests'],
    scripts=[os.path.join('bin', 'blast.py'),
             os.path.join('bin', 'import_genomes.py')],
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
)
