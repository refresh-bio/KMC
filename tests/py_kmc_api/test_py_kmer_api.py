#!/usr/bin/env python3
'''
A series of test for KmerAPI class.
'''
import sys
import os
import kmer_utils
import init_sys_path
import py_kmc_api as pka

def _py_kmer_api_from_string(kmer_str):
    kmer = pka.KmerAPI(len(kmer_str))
    kmer.from_string(kmer_str)
    return kmer

def test_kmer_from_string_to_string():
    ''' Tests if KmerAPI object is created correctly from
     string and then converted back to string. '''
    kmers = ('A', 'C', 'GT', 'TGC', 'ACGGTTAGG',
             'GCATCATGCAGTCTGAGCAACGTATGCTGAGCTGATGCTGACACTGATGCAAC')
    for kmer_str in kmers:
        kmer = _py_kmer_api_from_string(kmer_str)
        assert kmer_str == str(kmer)
        assert kmer_str == kmer.to_string()


def test_kmer_lt():
    ''' Tests if operator< (__lt__) works correctly. '''
    kmer1_str = 'ACGACGACG'
    kmer2_str = 'ACGACGACT'
    kmer1 = _py_kmer_api_from_string(kmer1_str)
    kmer2 = _py_kmer_api_from_string(kmer2_str)
    assert kmer1 < kmer2
    assert not kmer1 == kmer2
    assert not kmer2 < kmer1

def test_kmer_eq():
    ''' Tests if operator== (__eq__) works correctly. '''
    kmer_str = 'ACGACGACG'
    kmer1 = _py_kmer_api_from_string(kmer_str)
    kmer2 = _py_kmer_api_from_string(kmer_str)
    assert kmer1 == kmer2
    assert not kmer1 < kmer2
    assert not kmer2 < kmer1


def test_ascii_symbol():
    ''' Tests if ascii symbol of k-mer is correctly returned by ascii_symbol method. '''
    kmers = ('A', 'C', 'GT', 'TGC', 'ACGGTTAGG',
             'GCATCATGCAGTCTGAGCAACGTATGCTGAGCTGATGCTGACACTGATGCAAC')
    for kmer_str in kmers:
        kmer = _py_kmer_api_from_string(kmer_str)
        for i in range(len(kmer_str)):
            assert kmer.get_asci_symbol(i) == kmer_str[i]

def test_get_num_symbol():
    ''' Tests if get_num_bybol ('A' -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3) works correctly. '''
    kmers = ('A', 'C', 'GT', 'TGC', 'ACGGTTAGG',
             'GCATCATGCAGTCTGAGCAACGTATGCTGAGCTGATGCTGACACTGATGCAAC')
    mapping = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3,
    }
    for kmer_str in kmers:
        kmer = _py_kmer_api_from_string(kmer_str)
        for i in range(len(kmer_str)):
            assert kmer.get_num_symbol(i) == mapping[kmer_str[i]]


def test_reverse():
    ''' Tests if conversion of k-mer to its reverse complement works correctly. '''
    kmer_str = 'AAAAACGAAATTTA'
    kmer = _py_kmer_api_from_string(kmer_str)
    kmer.reverse()
    kmer_str_rev = str(kmer)
    assert kmer_str_rev == 'TAAATTTCGTTTTT'

def test_get_signature():
    ''' Test if valid signature is returned for k-mer. '''
    kmer_str = 'ACGGATGCAGTGCTAGCGGTGGCC'
    kmer = _py_kmer_api_from_string(kmer_str)
    sig = kmer.get_signature(7)
    assert sig == 1436
    kmer_str = 'AAAAAAAAAAAAAAAAGC'
    kmer = _py_kmer_api_from_string(kmer_str)
    for sig_len in (5, 11):
        assert kmer.get_signature(sig_len) == (1 << (2*sig_len))
