[![Build Status][]][]

KMC
===

KMC is a disk-based program for counting k-mers from (possibly gzipped) FASTQ/FASTA files. The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Requirements
============

-   gcc-4.7 or higher
-   cmake
-   make

Note that on Mac OSX, gcc defaults to clang which does not currently support OpenMP. OpenMP is an essential part of KMC so you will need to install gcc-4.7 (check out [homebrew][]). Once gcc is installed, make sure you export it before you attempt to build KMC.

    export CC=/usr/local/bin/gcc-4.9 && export CXX=/usr/local/bin/g++-4.9

Installation
============

To compile KMC:

    cd KMC
    mkdir build
    cd build
    cmake ..
    make

After compilation you will obtain two binaries:

  - bin/kmc - the main program for counting k-mer occurrences 
  - bin/kmc\_dump - the program listing k-mers in a database produced by KMC.

External Libraries
==================

The following libraries come with KMC in a binary (64-bit compiled for x86 platform) form. If your system needs other binary formats, you should put the following libraries in kmer\_counter/libs:

  - asmlib - for fast memcpy operation (http://www.agner.org/optimize/asmlib-instructions.pdf)
  - libbzip2 - for support for bzip2-compressed input FASTQ/FASTA files (http://www.bzip.org/)
  - zlib - for support for gzip-compressed input FASTQ/FASTA files (http://www.zlib.net/)

NOTE: asmlib is free only for non commercial purposes. If needed, you can contact the author of asmlib or compile KMC without asmlib. If you want to compile KMC without asmlib use `cmake -DADD_ASM_LIB=off ..` instead of `cmake ..`.

There are a few optimizations to also take note of:

  - You redefine maximal length of k-mer, which is 256 in the current version
  - KMC is highly optimized and spends only as many bytes for k-mer (rounded up to 8) as necessary, so using large values of MAX\_K does not affect the KMC performance for short k-mers
  - KMC creates a hundreds of temporary files, while default limit for opened files is small for under MAC OS platform. To increase this number use `ulimit -n 2048` before running KMC

Directory structure
===================

-   kmer\_counter - source code of KMC program
-   kmer\_counter/libs - compiled binary versions of libraries used by KMC
-   kmc\_api - C++ source codes implementing API; must be used by any program that wants to process databases produced by KMC
-   kmc\_dump - source codes of kmc\_dump program listing k-mers in databases produced by KMC

License
=======

-   KMC software distributed under GNU GPL 2 licence.
-   libbzip2 is open-source (BSD-style license)
-   gzip is free, open-source
-   asmlib is under the licence GNU GPL 3 or higher

Note: for commercial usage of asmlib follow the instructions in 'License conditions' (http://www.agner.org/optimize/asmlib-instructions.pdf) or compile KMC without asmlib. In case of doubt, please consult the original documentations.

Warranty
========

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

  [Build Status]: https://travis-ci.org/marekkokot/KMC.svg?branch=develop
  [![Build Status][]]: https://travis-ci.org/marekkokot/KMC
  [homebrew]: http://brew.sh
