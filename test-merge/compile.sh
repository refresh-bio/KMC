#!/bin/bash

git submodule update --init --recursive
(cd ..; make kmc kmc_tools -j)

g++ -O3 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive merge.cpp -I ../kmc_api/ ../kmc_api/*.cpp -o merge
g++ -O3 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive merge-v1.cpp -I ../kmc_api/ ../kmc_api/*.cpp -o merge-v1
g++ -O3 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive merge-v2.cpp -I ../kmc_api/ ../kmc_api/*.cpp -o merge-v2
g++ -O3 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive merge-v2-min-hash.cpp -I ../kmc_api/ ../kmc_api/*.cpp -o merge-v2-min-hash


g++ -O3 -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive test-ckmer-in-api.cpp -I ../kmc_api/ ../kmc_api/*.cpp -o test-ckmer-in-api

