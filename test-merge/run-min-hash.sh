#!/bin/bash

data=data.txt

(cd data; ./get-data.sh)

kmc=../bin/kmc
kmc_tools=../bin/kmc_tools

#In the case of min hash there is no need to generate common mapping because it is always mapped the same way
#lets also force reopening kmc temp file at each operation

$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024163_1.part.fastq.gz o1 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024163_2.part.fastq.gz o2 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024164_1.part.fastq.gz o3 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024164_2.part.fastq.gz o4 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024165_1.part.fastq.gz o5 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024165_2.part.fastq.gz o6 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024166_1.part.fastq.gz o7 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024166_2.part.fastq.gz o8 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024167_1.part.fastq.gz o9 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024167_2.part.fastq.gz o10 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024168_1.part.fastq.gz o11 .
$kmc -k28 -ci1 -sssmin_hash -v --reopen-tmp data/ERR024168_2.part.fastq.gz o12 .

#example how to merge, the code is in merge.cpp
#./merge merged.dump map.stbm o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12
t=16
./merge-v2-min-hash merged.dump o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12

#verification
echo "sorting"
sort merged.dump > merged.dump.sorted

$kmc -k28 -cs1000000 -ci1 @$data pat .
$kmc_tools transform pat dump -s pat.dump

diff merged.dump.sorted pat.dump

# we may also generate verification using min-hash, should give the same results
$kmc -k28 -sssmin_hash --reopen-tmp -cs1000000 -ci1 @$data pat_min_hash .
$kmc_tools transform pat dump -s pat_min_hash.dump

diff pat_min_hash.dump pat.dump
