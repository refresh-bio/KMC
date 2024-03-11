#!/bin/bash

data=data.txt

(cd data; ./get-data.sh)

kmc=../bin/kmc
kmc_tools=../bin/kmc_tools
#first generate mapping that will be used for each k-mer counting below
#use 5% of each input file
#the output mapping will be stored in map.stbm file
$kmc -k28 -ci1 --sig-to-bin-map-stats5 --only-generate-sig-to-bin-mappingmap.stbm @$data ignored .

#now use this mapping for each sample
# tell kmc to use map.stbm instead of computing maping by itself
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024163_1.part.fastq.gz o1 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024163_2.part.fastq.gz o2 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024164_1.part.fastq.gz o3 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024164_2.part.fastq.gz o4 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024165_1.part.fastq.gz o5 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024165_2.part.fastq.gz o6 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024166_1.part.fastq.gz o7 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024166_2.part.fastq.gz o8 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024167_1.part.fastq.gz o9 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024167_2.part.fastq.gz o10 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024168_1.part.fastq.gz o11 .
$kmc -k28 -ci1 --sig-to-bin-mappingmap.stbm data/ERR024168_2.part.fastq.gz o12 .

#example how to merge, the code is in merge.cpp
#./merge merged.dump map.stbm o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12
t=16
#./merge merged.dump map.stbm $t o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12

./merge-v2 merged.dump map.stbm o1 o2 o3 o4 o5 o6 o7 o8 o9 o10 o11 o12

#verification
echo "sorting"
sort merged.dump > merged.dump.sorted

$kmc -k28 -cs1000000 -ci1 @$data pat .
$kmc_tools transform pat dump -s pat.dump

diff merged.dump.sorted pat.dump
