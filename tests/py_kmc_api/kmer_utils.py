#!/usr/bin/env python3
''' Some helpful utility functions.'''


def rev_comp(kmer):
    ''' Gets rev. comp of a k-mer.'''
    res = ""
    mapping = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    for i in range(0, len(kmer)):
        res += mapping[kmer[-i - 1]]
    return res

def get_minimizer_no_canonical(kmer, minimizer_len):
    ''' Get no canonical minimizer of a k-mer.'''
    minimizer = kmer[0:minimizer_len]
    for i in range(1, len(kmer) - minimizer_len + 1):
        candidate = kmer[i:minimizer_len+i]
        if candidate < minimizer:
            minimizer = candidate
    return minimizer

def get_minimizer(kmer, minimizer_len):
    ''' Get canonical minimizer of a k-mer.'''
    rev = rev_comp(kmer)
    candidate1 = get_minimizer_no_canonical(kmer, minimizer_len)
    candidate2 = get_minimizer_no_canonical(rev, minimizer_len)
    return min(candidate1, candidate2)

def kmer_to_uint(kmer):
    ''' Convers k-mer string to int '''
    mapping = {
        'A': 0,
        'C': 1,
        'G': 2,
        'T': 3
    }
    res = 0
    for symb in kmer:
        res <<= 2
        res += mapping[symb]
    return res


