#!/bin/bash

data=data.txt

(cd data; ./get-data.sh)

kmc=../bin/kmc

for k in {13..100}; do
	$kmc -k$k -ci1 --reopen-tmp -sssmin_hash data/ERR024163_1.part.fastq.gz test .
	./test-ckmer-in-api test
done
