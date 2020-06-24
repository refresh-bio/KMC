KMC
=
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/kmc/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/KMC/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/famsa.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/kmc)

KMC is a disk-based programm for counting k-mers from (possibly gzipped) FASTQ/FASTA files.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

For accessing k-mers stored in database produced by KMC there is an API (kmc_api directory). Note that for KMC versions 0.x and 1.x dababase format differs from produced by KMC version 2.x. From version 2.2.0 API is unified for both formats and all new features/bug fixes are present only for 2.x branch (standalone API for older KMC version is not longer under development, so new version of API  should be used even for databases produced by older KMC version).

Installation
=
The following libraries come with KMC in a binary (64-bit compiled for x86 platform) form.
If your system needs other binary formats, you should put the following libraries in kmer_counter/libs:
* libbzip2 - for support for bzip2-compressed input FASTQ/FASTA files (http://www.bzip.org/)
* zlib - for support for gzip-compressed input FASTQ/FASTA files (http://www.zlib.net/)

The following libraries come with KMC in a source coude form.
 * pybind11 - used to create python wrapper of KMC API (https://github.com/pybind/pybind11)

If needed, you can also redefine maximal length of k-mer, which is 256 in the current version.

Note: KMC is highly optimized and spends only as many bytes for k-mer (rounded up to 8) as
necessary, so using large values of MAX_K does not affect the KMC performance for short k-mers.

Some parts of KMC use C++14 features, so you need a compatible C++ compiler, e.g., gcc 4.9+ or clang 3.4+

After that, you can run make to compile kmc and kmc_dump applications.

##### Additional infromation for MAC OS installation

For compilation under MAC OS there is makefile_mac.
Usage:

    make -f makefile_mac

There might be a need to change g++ path in makefile_mac.
If needed we recommend install g++ with brew (http://brew.sh/).

Note that KMC creates a hundreds of temporary files, while default limit for opened files is small for under MAC OS platform.
To increase this number use following command before running KMC:

    ulimit -n 2048

Directory structure
=
 * bin           - main directory of KMC (programs after compilation will be stored here)
 * kmer_counter  - source code of kmc program
 * kmer_counter/libs - compiled binary versions of libraries used by KMC
 * kmc_api       - C++ source codes implementing API; must be used by any program that wants to process databases produced by kmc
 * kmc_dump      - source codes of kmc_dump program listing k-mers in databases produced by kmc
 * py_kmc_api    - python wrapper for kmc API

Python wrapper for KMC API
=
Python wrapper for KMC API was created using pybind11.
**Warning:** python binding is experimental. The library used to create binding as well as public interface may change in the future.
**Warning 2:** python wrapper for C++ KMC API is **much** slower (much, much more than I have been expecting) than native C++ API. In fact the first attempt to create python wrapper was to use [ctypes](https://docs.python.org/2/library/ctypes.html), but it turned out it was even slower than in case when pybind11 is used.
The wrapper is designed and was tested only for python3. The main goal was to make it as similar to C++ API as possible. For this reason the API may be not [pythonic]
(https://blog.startifact.com/posts/older/what-is-pythonic.html) enough for regular python programmer. Suggestions or pull requests to make it more robust are welcome.

Python module wrapping KMC API must be compiled.
 * for windows there is a visual studio project (note that there will be probably the need to change include directories and library directories to point python include and libs location)
 * for linux one should run ```make py_kmc_api```
 * for mac on should run ```make -f makefile_mac py_kmc_api```

As a result of pybind11 *.so file (for linux and mac os) or *.pyd (for windows) is created and may be used as a python module. *.pyd file is in fact DLL file, the only difference is its extension.
  * for windows following file is created: x64/Release/py_kmc_api.pyd
  * for linux/mac os the following file is created: bin/py_kmc_api`python3-config --extension-suffix`

To be able to use this file one should make it visible for python. One way to do this is to extend PYTHONPATH environment variable.
For linux/mac os one may just
```
source py_kmc_api/set_path.sh
```
while, for windows:
```
py_kmc_api\set_path.bat
```
it will export apropriate file. The example of Python wrapper for KMC API is presented in file:
py_kmc_api/py_kmc_dump.py

Detailed API describtion is avaiable at [wiki](https://github.com/refresh-bio/KMC/wiki/Python-wrapper-for-KMC-API)


Binaries
=
After compilation you will obtain two binaries:
* bin/kmc - the main program for counting k-mer occurrences
* bin/kmc_dump - the program listing k-mers in a database produced by kmc
* bin/kmc_tools - the program allowing to manipulate kmc databases (set operations, transformations, etc.)


License
=
* KMC software distributed under GNU GPL 3 licence.

* libbzip2 is open-source (BSD-style license)

* gzip is free, open-source

* pybind11 (https://github.com/pybind/pybind11) is open-source (BDS-style license)

In case of doubt, please consult the original documentations.

Warranty
=
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING
THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

