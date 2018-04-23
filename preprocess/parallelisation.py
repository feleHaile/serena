#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 11 11:35:50 2018

@author: je
"""

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
"""

import concurrent.futures
import multiprocessing

# Submit the jobs to the thread pool executor.
with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
# Map the futures returned from executor.submit()
# to their destination windows.
#
# The _example.compute function modifies no Python
# objects and releases the GIL. It can execute
# concurrently.
    future_to_window = {executor.submit(compute, data, res): (res, window) for data, res, window in jobs()}
# As the processing jobs are completed, get the
# results and write the data to the appropriate
# destination window.
    for future in concurrent.futures.as_completed(future_to_window):
        result, window = future_to_window[future]
        dst.write(arr, window=window)
        
#    def work(folderName):
#        """
#        Thread Function
#        """
#        dnToTOA(os.path.join(inPath,folderName), outPath, sensor)
#    
#    NUM_WORKERS = 7
# 
#    start_time = time.time()
#     
#    with concurrent.futures.ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
#        
#        futures = {executor.submit(work, folderName) for folderName in lstFolders}
#        
#        concurrent.futures.wait(futures)
#     
#    end_time = time.time()        
#     
#    print("Time for proccess : %ssecs" % (end_time - start_time))