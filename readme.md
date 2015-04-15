#Demeter

Dependencies:
* Python 2.7
* MongoDB
  * pymongo
* BioPython

### Identification of horizontal gene transfer between sequenced microbial genomes

####DataImport:
Import genbank-formated annotated genomes into database

####DataTypes
Framework for loading genomic information, broken down into species, contig/chromasome, gene, and groups of related genes (HGT_island)

####KsavirBlast
* Compares each gene in a genome to those of other genomes in database
* Of hits, combines those within some number of bases (to identify hits in multi-gene islands)
* Output of genes, species and islands identified as HGT
