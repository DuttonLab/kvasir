#Kvasir

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53206.svg)](http://dx.doi.org/10.5281/zenodo.53206) [![Join the chat at https://gitter.im/kescobo/kvasir](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/kescobo/kvasir?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Dependencies:
* Python 2.7
* MongoDB
  * pymongo
* BioPython
* BLAST+ CLI

**Identification of horizontal gene transfer between sequenced microbial genomes**

Kvasir takes as the input a folder containing genomes in genbank format. The protein coding genes from these genomes are loaded into a database, and blasted against each other.

### Usage

Change the values in `settings.py` to point at your input folder, output folder, and the name you want for your database.

Launch a local `mongod` instance:
```
$ mongod --dbpath path/to/db
```

Run functions in `run.py`. Eventually, this will get more streamlined, for now...

```
Python 2.7.11 (default, Dec 14 2015, 10:44:13)
[GCC 4.2.1 Compatible Apple LLVM 7.0.2 (clang-700.1.81)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import run
>>> run.import_data()
>>> run.blast_db()
>>> run.analyze(0.99)
```
