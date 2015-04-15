##!/usr/bin/env python

from pymongo import MongoClient
from subprocess import Popen
from subprocess import _active
import os
import sys
import atexit


from subprocess import Popen
import atexit

started = []

def auto_popen(*args, **kw):
    p = Popen(*args, **kw)
    started.append(p)
    return p

def testing():
    mongod = auto_popen(["mongod", "--dbpath", '/Users/KBLaptop/computation/db/'], shell=True)

    print mongod.pid
    client = MongoClient()
    db = client['location_import_test']
    collection = db['test_collection']

    all_species = db.collection_names(False)
        
    for species in all_species:
        current_species_collection = db[species]
        for gene in current_species_collection.find():
            print gene['locus_taag']
            pass

    mongod.terminate()

def cleanup():
    for proc in started:
        if proc.poll() is not None:
            try: proc.kill()
            except: pass

atexit.register(cleanup)
testing()
