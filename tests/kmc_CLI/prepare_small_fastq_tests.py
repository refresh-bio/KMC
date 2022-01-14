#!/usr/bin/env python3

import subprocess
import os
import sys
import queue
import threading

if len(sys.argv) < 5:
    print("Usage: {} <data_dir> <trivial k-mer counter path> <fileToGetPartFrom> <no. lines>".format(sys.argv[0]))
    sys.exit(1)

script_path = os.path.dirname(os.path.abspath(__file__))

data_dir = sys.argv[1]
count_progam = sys.argv[2]
fileToGetPartFrom = os.path.join(data_dir, sys.argv[3])
no_lines = sys.argv[4]

partFilePath = os.path.join(data_dir, "cached", "small.fastq", os.path.basename(fileToGetPartFrom) + ".part.{}.fq".format(no_lines))
if not os.path.exists(os.path.dirname(partFilePath)):
    os.makedirs(os.path.dirname(partFilePath))
if not os.path.exists(partFilePath):
    subprocess.call(["zcat {} | head -n {} > {}".format(fileToGetPartFrom, no_lines, partFilePath)], shell=True)

q=queue.Queue()
print_q = queue.Queue()


def count(k, ci, cx, cs):
    o = os.path.join(partFilePath + ".pattern-res", "k{}.ci{}.cx{}.cs{}".format(k, ci, cx, cs))    
    if not os.path.exists(o):
        os.makedirs(o)    
    o = os.path.join(o, "kmers")
    if not os.path.exists(o) or not os.path.exists(o+".stats"):
        cmd = "{} -k {} -ci {} -cx {} -cs {} {} {}".format(count_progam, k, ci, cx, cs, partFilePath, o)
        print_q.put((cmd))        
        out = open(o + ".stdout", 'w')
        err = open(o + ".stderr", 'w')
        subprocess.call([cmd], shell=True,stdout=out, stderr=err)


def print_worker():
	while True:
		(cmd) = print_q.get()
		print(cmd)
		print_q.task_done()

def worker():
	while True:
		(k, ci, cx, cs) = q.get()
		count(k, ci, cx, cs)
		q.task_done()

num_thread = 32
for i in range(0,num_thread):
	t=threading.Thread(target=worker)
	t.daemon=True
	t.start()

t=threading.Thread(target=print_worker)
t.daemon=True
t.start()

for k in range(1, 257):
    q.put((k, 2, 1000000000, 255))    
	
q.join()
print_q.join()