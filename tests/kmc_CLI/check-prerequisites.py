#!/usr/bin/env python3

import hashlib
import os
import sys
import subprocess

if len(sys.argv) < 2:
    print("Usage: {} <data_dir>".format(sys.argv[0]))
    sys.exit(1)

script_path = os.path.dirname(os.path.abspath(__file__))

def compute_md5(path):
    with open(path, "rb") as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)        
        return file_hash.hexdigest()

def verify_md5(baseDir, filePath, pattern):    
    cached_path_md5_file = os.path.join(baseDir, "cached", "md5", filePath + ".md5")    
    
    if not os.path.exists(os.path.dirname(cached_path_md5_file)):
        os.makedirs(os.path.dirname(cached_path_md5_file))
    if not os.path.exists(cached_path_md5_file):        
        md5 = compute_md5(os.path.join(baseDir, filePath))
        with open (cached_path_md5_file, "w") as f:
            f.write(md5)
    else:
        with open (cached_path_md5_file) as f:
            md5 = f.read()

    return md5 == pattern

with open(os.path.join(script_path, "prerequisite-files.txt")) as f:
    for line in f:
        line = line.strip()        
        path = line[0:-33]
        md5 = line[-32:]
        
        full_path = os.path.join(sys.argv[1], path)
        if not os.path.exists(full_path):
            print("{} does not exists".format(full_path))
            sys.exit(1)
        print("checking {}".format(full_path))
    
        if not verify_md5(sys.argv[1], path, md5):        
            print("Invalid md5 for {}".format(full_path))
            sys.exit(1)