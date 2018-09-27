#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:35:50 2018

@author: je

"""

######################### WITH MULTIPROCESSING ##################################

from multiprocessing import Pool

def process_file (filename) :
    #... traiter ici le fichier filename
    return filename

pool = Pool(processes=4)
files = []
results = pool.imap( process_file , files )

for result in results :
    print (result)

######################### WITH CONCURRENT FUTURES ##################################

import concurrent.futures
import time 

NUM_WORKERS = 4

lstFolders= []

def work(folderName):
    """
    Thread Function
    """
    #... traiter ici le dossier foldername

start_time = time.time()

# Submit the jobs to the thread pool executor.
with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
# Map the futures returned from executor.submit() to their destination windows.
# The _example.compute function modifies no Python objects and releases the GIL. It can execute concurrently.
    futures = {executor.submit(work, folderName) for folderName in lstFolders}
    concurrent.futures.wait(futures)
    
end_time = time.time()        
     
print("Time for proccess : %ssecs" % (end_time - start_time))
    
    
 

     
        
     
    