#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

import sys
import KvDataStructures as kv

kv.mongo_init(sys.argv[1])
kv.reset_database(kv.db.name)