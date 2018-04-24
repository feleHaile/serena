#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:35:50 2018

@author: je
"""

from multiprocessing import Pool
def process_file (filename) :
    #... traiter ici le fichier filename
    return filename

pool = Pool()
files = [file1, file2, file3, file4, file5, file6, file7]
results = pool.imap( process_file , files )

for result in results :
    print (result)