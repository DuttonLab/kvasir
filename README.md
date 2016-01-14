#Kvasir

[![Join the chat at https://gitter.im/kescobo/kvasir](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/kescobo/kvasir?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Dependencies:
* Python 2.7
* MongoDB
  * pymongo
* BioPython
* ~~Scikit Bio~~
* ~~Pandas~~

**Identification of horizontal gene transfer between sequenced microbial genomes**

##Running Kvasir
With dependencies installed, fire up a Mongod instance. In the terminal:

`mongod --dbpath path/to/db`

Run Kvasir by invoking run_kvasir.py in your working directory:

`python run_kvasir.py /path/to/gb_files name_of_mongoDB`

##The following is out of date... Will get back to this soon

####~~DataImport~~:
* ~~Imports genbank-formated annotated genomes into Mongo database.
* ~~.gb files require "locus_tag" feature. If your genomes don't have it, FixGbk.py shoul take care of it for you~~
* ~~Mongo database has "collections" and "documents" - a different collection is generated for each species (each separate genbank file) and documents representing each CDS. ~~
    * ~~CDS documents are like python dictionaries, and contain entries for species, DNA and amino acid protein sequences, contig and location info, and annotation information.~~
    * ~~each document is assigned a uniqe `_id` attribute within the species, so every gene is uniquely identified by a `(species, _id)` tuple~~

####~~MakeBlastDb~~
~~Generates a multi-fasta file containing every gene in the mongo database, generates a BLASTable database, then deletes the temprorary file~~

####~~KsavirBlast~~
* ~~For each species, generates a temporary fasta file and BLASTs against every other gene in the database~~
* ~~BLAST generates and xml document, which is parsed for unique hits~~
* ~~new "hits" entry is added to each gene document in MongoDB, which contains a list of `(species, _id)` tuples for each hit (these are used in the next script to gather info about hits)~~

####~~Outputs~~
~~Still a work in progress. So far, have a bunch of output formates working... will detail later.~~
