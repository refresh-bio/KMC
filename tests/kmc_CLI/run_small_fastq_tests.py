#!/usr/bin/env python3

import subprocess
import os
import sys
import filecmp

if len(sys.argv) < 5:
    print("Usage: {} <data_dir> <kmc_path> <kmc_tools_path> <kmc_dump_path>".format(sys.argv[0]))
    sys.exit(1)

script_path = os.path.dirname(os.path.abspath(__file__))

data_dir = sys.argv[1]
kmc = sys.argv[2]
kmc_tools = sys.argv[3]
kmc_dump = sys.argv[4]

inputDir = os.path.join(data_dir, "cached", "small.fastq")

def read_params_from_dir_name(dir_name):
    parts = dir_name.split(".")
    res = ""
    for part in parts:
        if part.startswith("k"):
            res += " -k{}".format(part[1:])
        elif part.startswith("ci"):
            res += " -ci{}".format(part[2:])
        elif part.startswith("cx"):
            res += " -cx{}".format(part[2:])            
        elif part.startswith("cs"):
            res += " -cs{}".format(part[2:])            
        elif part.startswith("b"):
            res += " -b"
    return res

def is_sorted(kmc_db):
    command = "{} info o".format(kmc_tools)
    print(command)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode("utf-8")
    for line in stdout.split('\n'):
        if "database format" in line:
            return "KMC1" in line
    print("Error: cannot read kmc db format")
    sys.exit(1)

def compare_stats(kmc_stats, pattern_file_path):
    kmc = dict()
    pattern = dict()
    for line in kmc_stats.split('\n'):
        if "No. of k-mers below min. threshold" in line:
            kmc["n_ci"] = line.split(":")[-1].strip()
        elif "No. of k-mers above max. threshold" in line:
            kmc["n_cx"] = line.split(":")[-1].strip()
        elif "No. of unique k-mers" in line:
            kmc["n_unique"] = line.split(":")[-1].strip()
        elif "No. of unique counted k-mers" in line:
            kmc["n_unique_cnt"] = line.split(":")[-1].strip()
        elif "Total no. of k-mers" in line:
            kmc["n_kmers"] = line.split(":")[-1].strip()
        elif "Total no. of reads" in line:
            kmc["n_seqs"] = line.split(":")[-1].strip()

    with open(pattern_file_path) as f:
        for line in f:
            line = line.strip()
            if "Total unique k-mers below min cutoff:" in line:
                pattern["n_ci"] = line.split(":")[-1].strip()
            elif "Total unique k-mers above max cutoff:" in line:
                pattern["n_cx"] = line.split(":")[-1].strip()
            elif "Total unique k-mers:" in line:
                pattern["n_unique"] = line.split(":")[-1].strip()

            elif "Total unique counted k-mers:" in line:
                pattern["n_unique_cnt"] = line.split(":")[-1].strip()
            elif "Total k-mers:" in line:
                pattern["n_kmers"] = line.split(":")[-1].strip()
            elif "Total sequences:" in line:
                pattern["n_seqs"] = line.split(":")[-1].strip()                                
    return kmc == pattern

def run_for_file(input):
    pattern_dir = input + ".pattern-res"
    dirs_to_run = []
    for x in os.listdir(pattern_dir):
        if os.path.isfile(os.path.join(pattern_dir, x)):
            print("Unexpected file " + os.path.join(pattern_dir, x))
            sys.exit(1)
        dirs_to_run.append(os.path.basename(x))
    
    dirs_to_run.sort()
    for x in dirs_to_run:
        kmc_params = read_params_from_dir_name(x)
        command = "/usr/bin/time -v {} {} -v -hp {} o .".format(kmc, kmc_params, input)
        print(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True)
        stdout, stderr = proc.communicate()
        stdout = stdout.decode("utf-8")
        print(stderr.decode("utf-8"))
        print(stdout)
        

        if not compare_stats(stdout, os.path.join(pattern_dir, x, "kmers.stats")):
            print("Error: summary statistics does not match")
            sys.exit(1)

        db = "o"

        pattern_dump_path = os.path.join(pattern_dir, x, "kmers")
        if not is_sorted("o"):
            command = "/usr/bin/time -v {} -hp transform o sort o.sorted".format(kmc_tools)
            print(command)
            subprocess.call([command], shell = True)
            db = "o.sorted"

            #verify kmc_tools transform db dump -s db.sorted.dump
            command = "/usr/bin/time -v {} -hp transform o dump -s o.sorted.dump".format(kmc_tools)
            print(command)
            subprocess.call([command], shell = True)

            print("compare o.sorted.dump (from kmc_tools dump -s) and {}".format(pattern_dump_path))
            if not filecmp.cmp("o.sorted.dump", pattern_dump_path):
                print("Error: o.sorted.dump (from kmc_tools dump -s) and {} differs".format(db, pattern_dump_path))
                sys.exit(1)

        command = "/usr/bin/time -v {} -hp transform {} dump {}.dump".format(kmc_tools, db, db)
        print(command)
        subprocess.call([command], shell = True)
                
        print("compare {}.dump and {}".format(db, pattern_dump_path))
        if not filecmp.cmp("{}.dump".format(db), pattern_dump_path):
            print("{}.dump and {} differs".format(db, pattern_dump_path))
            sys.exit(1)

        command = "/usr/bin/time -v {} {} {}.dump".format(kmc_dump, db, db)
        print(command)
        subprocess.call([command], shell = True)

        print("compare {}.dump (from kmc_dump) and {}".format(db, pattern_dump_path))
        if not filecmp.cmp("{}.dump".format(db), pattern_dump_path):
            print("{}.dump (from kmc_dump) and {} differs".format(db, pattern_dump_path))
            sys.exit(1)
        
for x in os.listdir(inputDir):
    fullPath = os.path.join(inputDir, x)
    if not os.path.isfile(fullPath):
        continue
    run_for_file(fullPath)    

