KMC
=
[![GitHub downloads](https://img.shields.io/github/downloads/refresh-bio/kmc/total.svg?style=flag&label=GitHub%20downloads)](https://github.com/refresh-bio/KMC/releases)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/kmc.svg?style=flag&label=Bioconda%20downloads)](https://anaconda.org/bioconda/kmc)
[![Biocontainer downloads](https://img.shields.io/endpoint?url=https%3A%2F%2Fmmseqs.com%2Fbiocontainer.php%3Fcontainer%3Dkmc)](https://biocontainers.pro/tools/kmc)
[![GitHub Actions CI](../../actions/workflows/main.yml/badge.svg)](../../actions/workflows/main.yml) [![Join the chat at https://gitter.im/refresh-bio/KMC](https://badges.gitter.im/refresh-bio/KMC.svg)](https://gitter.im/refresh-bio/KMC?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

KMC is a disk-based program for counting k-mers from (possibly gzipped) FASTQ/FASTA files.
KMC is one of many projects developed by [REFRESH Bioinformatics Group](http://sun.aei.polsl.pl/REFRESH/).

For accessing k-mers stored in database produced by KMC there is an API (kmc_api directory). Note that for KMC versions 0.x and 1.x dababase format differs from produced by KMC version 2.x. From version 2.2.0 API is unified for both formats and all new features/bug fixes are present only for 2.x branch (standalone API for older KMC version is not longer under development, so new version of API  should be used even for databases produced by older KMC version).

Quick start
=
#### Getting the executable
The simplest way to get the KMC is to download newest release for appropriate operating system from [KMC releases](https://github.com/refresh-bio/KMC/releases).
#### Counting the k-mers from a single fastq file
```
./kmc -k27 input.fastq 27mers .
```
The command above will count all the 27-mers occurring in ```input.fastq``` at least twice (configurable with ```-ci``` switch). 
The result will be stored in a KMC database, which is split into two files: ```27mers.kmc_pre``` and ```27mers.kmc_suf```.
KMC will create hundreds of intermediate files. In the case of the above command, those will be created in the current working directory(the ```.``` at the end of the command). 
It may be more convinient to use dedicated directory for KMC temporary files, for example:
```
mkdir kmc_tmp # create directory for kmc temporary files
./kmc -k27 input.fastq 27mers kmc_tmp
```
#### Create text dump from KMC database binary format
Having the k-mers counted it is possible to dump KMC binary database to textual form with ```kmc_tools```.
```
./kmc_tools transform 27mers dump 27mers.txt
```

Installation details
##### Compile from sources
```
git clone --recurse-submodules https://github.com/refresh-bio/KMC.git
make -j32
```
=
The following libraries come with KMC in a binary (64-bit compiled for x86 platform) form.
If your system needs other binary formats, you should put the following libraries in kmc_core/libs:
* zlib - for support for gzip-compressed input FASTQ/FASTA files

The following libraries come with KMC in a source coude form.
 * pybind11 - used to create python wrapper of KMC API (https://github.com/pybind/pybind11)

If needed, you can also redefine maximal length of k-mer, which is 256 in the current version.

Note: KMC is highly optimized and spends only as many bytes for k-mer (rounded up to 8) as
necessary, so using large values of MAX_K does not affect the KMC performance for short k-mers.

Some parts of KMC use C++17 features, so you need a compatible C++ compiler

After that, you can run make to compile kmc and kmc_dump applications.

##### Additional infromation for MAC OS installation

There might be a need to change g++ path in makefile_mac.
If needed we recommend install g++ with brew (http://brew.sh/).

Note that KMC creates a hundreds of temporary files, while default limit for opened files is small for under MAC OS platform.
To increase this number use following command before running KMC:

    ulimit -n 2048

Directory structure
=
 * bin           - after compilation executables and libraries after compilation will be stored here
 * include       - after compilation header file to use kmc core through the C++ API will be stored here
 * kmc_core      - source code of kmc core library
 * kmc_CLI       - source code of kmc command line interface
 * kmc_tools     - source codes of kmc_tools program
 * kmc_core/libs - libraries used by KMC
 * kmc_api       - C++ source codes implementing API to access KMC databases; must be used by any program that wants to process databases produced by kmc
 * kmc_dump      - source codes of kmc_dump program listing k-mers in databases produced by kmc (deprecated, use kmc_tools instead)
 * py_kmc_api    - python wrapper for kmc API
 * tests         - tests files

Use the KMC directly from code through the API
=
It is possible to use the KMC directly from C++ code through.
Detailed API description is available at [wiki](https://github.com/refresh-bio/KMC/wiki/Use-the-KMC-directly-from-code-through-the-API)

Python wrapper for KMC API
=
Python wrapper for KMC API was created using pybind11.
**Warning:** python binding is experimental. The library used to create binding as well as public interface may change in the future.
**Warning 2:** python wrapper for C++ KMC API is **much** slower (much, much more than I have been expecting) than native C++ API. In fact the first attempt to create python wrapper was to use [ctypes](https://docs.python.org/2/library/ctypes.html), but it turned out it was even slower than in case when pybind11 is used.
The wrapper is designed and was tested only for python3. The main goal was to make it as similar to C++ API as possible. For this reason the API may be not [pythonic]
(https://blog.startifact.com/posts/older/what-is-pythonic.html) enough for regular python programmer. Suggestions or pull requests to make it more robust are welcome.

Python module wrapping KMC API must be compiled.
 * for windows there is a visual studio project (note that there will be probably the need to change include directories and library directories to point python include and libs location)
 * for linux or mac one should run ```make py_kmc_api```

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

Detailed API description is available at [wiki](https://github.com/refresh-bio/KMC/wiki/Python-wrapper-for-KMC-API)

Binaries
=
After compilation you will obtain two binaries:
* bin/kmc - the main program for counting k-mer occurrences
* bin/kmc_dump - the program listing k-mers in a database produced by kmc
* bin/kmc_tools - the program allowing to manipulate kmc databases (set operations, transformations, etc.)
* bin/libkmc_core.a - compiled KMC code sources
* py_kmc_api.cpython-39-x86_64-linux-gnu.so - compiled python wrapper for KMC API


License
=
* KMC software distributed under GNU GPL 3 licence.

* gzip is free, open-source

* pybind11 (https://github.com/pybind/pybind11) is open-source (BDS-style license)

In case of doubt, please consult the original documentations.

Archival source codes, binaries and documentation
= 
Archival source codes, binaries and documentation are available at [wiki](https://github.com/refresh-bio/KMC/wiki/Archive-versions-of-sources,-binaries-and-documentations).

Warranty
=
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE,
TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING
THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Citing
[Marek Kokot, Maciej Długosz, Sebastian Deorowicz, KMC 3: counting and manipulating k-mer statistics, Bioinformatics, Volume 33, Issue 17, 01 September 2017, Pages 2759–2761, https://doi.org/10.1093/bioinformatics/btx304](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399)

[Sebastian Deorowicz, Marek Kokot, Szymon Grabowski, Agnieszka Debudaj-Grabysz, KMC 2: fast and resource-frugal k-mer counting, Bioinformatics, Volume 31, Issue 10, 15 May 2015, Pages 1569–1576, https://doi.org/10.1093/bioinformatics/btv022](https://academic.oup.com/bioinformatics/article/31/10/1569/177467)

[Deorowicz, S., Debudaj-Grabysz, A. & Grabowski, S. Disk-based k-mer counting on a PC. BMC Bioinformatics 14, 160 (2013). https://doi.org/10.1186/1471-2105-14-160](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-160)
