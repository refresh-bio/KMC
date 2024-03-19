#!/bin/bash

if [ -e data.tar.gz ]; then
	echo "already downloaded"
else
    wget wget http://mkokot.ddns.net:8088/kmc-bin-ordering-exmaple-data/data.tar.gz
	tar -xvf data.tar.gz
fi
