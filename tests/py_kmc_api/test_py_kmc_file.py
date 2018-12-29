#!/usr/bin/env python3
'''
A series of test for PyKMCFile class.
'''

import sys
import os
import subprocess
import kmer_utils
import init_sys_path
import py_kmc_api as pka
import pytest




@pytest.fixture(scope="module", autouse=True)
def create_kmc_db():
    '''
    Set up tests and clean up after.
    '''
    kmer_len = 17
    memory = 1 #GB
    cutoff_min = 1
    sig_len = 9
    reads_src = 'input.fastq'
    reads = (
        'GGCATTGCATGCAGTNNCAGTCATGCAGTCAGGCAGTCATGGCATGCAACGACGATCAGTCATGGTCGAG',
        'GGCATTGCATGCAGTNNCAGTCATGCAGTCAGGCAGTCATGGCATGCAACGACGATCAGTCATGGTCGAG',
        'GTCGATGCATCGATGCTGATGCTGCTGTGCTAGTAGCGTCTGAGGGCTA'
    )
    _save_reads_as_fastq(reads, reads_src)

    kmers = _cout_kmers(reads, kmer_len)
    absent_kmers = _generate_not_existing_kmers(kmers, kmer_len)
    _run_kmc(cutoff_min, kmer_len, memory, sig_len, reads_src)

    result = {
        'kmers': kmers,
        'kmer_len': kmer_len,
        'sig_len': sig_len,
        'absent_kmers': absent_kmers
    }
    yield result

    os.remove(reads_src)
    os.remove('kmc_db.kmc_pre')
    os.remove('kmc_db.kmc_suf')

def _cout_kmers(reads, kmer_len):
    ''' Simple k-mer counting routine. '''
    kmers = {}
    for read in reads:
        for start in range(0, len(read) - kmer_len + 1):
            kmer = read[start:start+kmer_len]
            if 'N' in kmer:
                continue
            rev = kmer_utils.rev_comp(kmer)
            if rev < kmer:
                kmer = rev
            if kmer in kmers.keys():
                kmers[kmer] += 1
            else:
                kmers[kmer] = 1
    return kmers

def _save_reads_as_fastq(reads, file_name):
    ''' Save reads from input to file named file_name. '''
    file = open(file_name, 'w')
    for read in reads:
        file.write("@TEST\n")
        file.write(read + "\n")
        file.write("+TEST\n")
        file.write("I"*len(read) + "\n")
    file.close()


def _generate_not_existing_kmers(kmers, kmer_len):
    ''' Generate k-mers that are not present in the database.

        :kmers: existing k-mers
        :kmer_len: length of k-mers

    '''
    def increment_kmer(kmer, start):
        ''' Increments k-mer to next lexographical.

        Start from pos :start: (from end, i.e. start = 0 means last k-mer symbol). '''
        def replace_char(string, pos, new_char):
            ''' Create new string with character at :pos: changed to :new_char:. '''
            if pos < 0:
                pos = len(string) + pos
            return string[:pos] + new_char + string[pos+1:]
        for i in range(start, len(kmer)):
            if kmer[-1-i] == 'A':
                return replace_char(kmer, -1 - i, 'C')
            if kmer[-1-i] == 'C':
                return replace_char(kmer, -1 - i, 'G')
            if kmer[-1-i] == 'T':
                return replace_char(kmer, -1 - i, 'T')
            kmer = replace_char(kmer, -1 - i, 'T')
        return kmer

    absent_kmers = []

    for i in range(0, kmer_len):
        for kmer_str in kmers.keys():
            inc_kmer = increment_kmer(kmer_str, i)
            if not inc_kmer in kmers.keys():
                absent_kmers.append(inc_kmer)
    return absent_kmers

def _run_kmc(cutoff_min, kmer_len, memory, sig_len, reads_src):
    ''' Runs kmc. '''
    if init_sys_path.is_linux() or init_sys_path.is_mac():
        kmc_path = os.path.join(os.path.dirname(__file__), '../../bin/kmc')
    elif init_sys_path.is_windows():
        kmc_path = os.path.join(os.path.dirname(__file__), '../../x64/Release/kmer_counter.exe')

    subprocess.call([kmc_path,
                     '-ci{}'.format(cutoff_min),
                     '-k{}'.format(kmer_len),
                     '-m{}'.format(memory),
                     '-p{}'.format(sig_len),
                     reads_src,
                     'kmc_db',
                     '.'
                    ])



def _open_for_listing():
    ''' Open kmc database for listing and check if opened sucessfully. '''
    kmc_file = pka.KMCFile()
    assert kmc_file.OpenForListing('kmc_db')
    return kmc_file

def _open_for_ra():
    ''' Open kmc database for random access and check if opened sucessfully. '''
    kmc_file = pka.KMCFile()
    assert kmc_file.OpenForRA('kmc_db')
    return kmc_file

def test_info(create_kmc_db):
    '''
    Test if some fields in object returned from Info are set properly.

    '''
    pattern = create_kmc_db
    kmc_file = _open_for_listing()
    info = kmc_file.Info()
    assert info.kmer_length == pattern['kmer_len']
    assert info.mode == 0 # no Quake mode (quake not supported anymore)
    assert info.counter_size == 1
    assert info.signature_len == pattern['sig_len']
    assert info.min_count == 1
    assert info.both_strands
    assert info.total_kmers == len(pattern['kmers'])

def test_kmc_file_next_kmer(create_kmc_db):
    ''' Test if all counted k-mers are returned by KMC API using NextKmer method. '''
    pattern = create_kmc_db['kmers']
    kmc_file = _open_for_listing()
    counter = pka.Count()
    kmer = pka.KmerAPI(create_kmc_db['kmer_len'])
    res = {}
    while kmc_file.ReadNextKmer(kmer, counter):
        res[str(kmer)] = counter.value
    assert res == pattern

def test_get_counters_for_read(create_kmc_db):
    ''' Test case for GetCountersForRead method of KMCFile. '''
    kmers = create_kmc_db['kmers']
    read = "GGCATTGCATGCAGTNNCAGTCATGCAGTCAGGCAGTCATGGCATGCGTAAACGACGATCAGTCATGGTCGAG"
    pattern = []
    kmer_len = create_kmc_db['kmer_len']
    for i in range(0, len(read) - kmer_len + 1):
        kmer = read[i:i+kmer_len]
        if 'N' in kmer:
            pattern.append(0)
            continue
        rev = kmer_utils.rev_comp(kmer)
        if rev < kmer:
            kmer = rev
        if not kmer in kmers.keys():
            pattern.append(0)
        else:
            pattern.append(kmers[kmer])

    kmc_file = _open_for_ra()
    res = pka.CountVec()
    kmc_file.GetCountersForRead(read, res)
    assert res.value == pattern


def test_check_kmer(create_kmc_db):
    '''
    Test case for CheckKmer method.

    Check if are k-mers from input are present in the database and
    if some not present in the input are absent in the output.
    '''
    kmers = create_kmc_db['kmers']
    kmer_len = create_kmc_db['kmer_len']
    kmer = pka.KmerAPI(kmer_len)
    counter = pka.Count()
    kmc_file = _open_for_ra()
    for kmer_str, count in kmers.items():
        kmer.from_string(kmer_str)
        assert kmc_file.CheckKmer(kmer, counter)
        assert counter.value == count
    absent_kmers = create_kmc_db['absent_kmers']
    for kmer_str in absent_kmers:
        kmer.from_string(kmer_str)
        assert not kmc_file.CheckKmer(kmer, counter)
