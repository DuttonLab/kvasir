##!/usr/bin/env python

from pymongo import MongoClient
from subprocess import Popen
from subprocess import _active
import os
import sys
import atexit

started = []

def auto_popen(*args, **kw):
    p = Popen(*args, **kw)
    started.append(p)
    return p

def testing():
    mongod = auto_popen(
            ["mongod", "--dbpath", '/Users/KBLaptop/computation/db/'],
        )
    print mongod.pid
    
    client = MongoClient()
    db = client['location_import_test']
    collection = db['test_collection']

    all_species = db.collection_names(False)
        
    for species in all_species:
        current_species_collection = db[species]
        for gene in current_species_collection.find():
            print gene['locus_taag']

    print mongod.pid
    mongod.terminate()

def cleanup():
    print 'Hi there'
    for proc in started:
        print proc
        print proc.poll()
        if proc.poll() is not None:
            print 'hello'
            try: proc.kill()
            except: pass

atexit.register(cleanup)
testing()