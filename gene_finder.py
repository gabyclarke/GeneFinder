# -*- coding: utf-8 -*-
"""
CODE FOR SOFTDES MINI PROJECT 1: GENE FINDER
SPRING 2015

@author: Gaby Clarke

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

stop_codons = ['TAA', 'TAG', 'TGA']


def shuffle_string(s):
    """Shuffles the characters in the input string"""

    return ''.join(random.sample(s, len(s)))


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    complements = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return complements[nucleotide]


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    reverse = dna[::-1]
    reverse_complement = ''

    for i in reverse:
        reverse_complement += get_complement(i)
    return reverse_complement


def get_codons(dna):
    """ Takes a DNA sequence and breaks it into codons.

        dna: a DNA sequence
        returns: a list of codons
    >>> get_codons("ATGTGAA")
    ['ATG', 'TGA', 'A']
    >>> get_codons("ATGAAATGA")
    ['ATG', 'AAA', 'TGA']
    """

    return [dna[i:i + 3] for i in range(0, len(dna), 3)]



def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """

    codons = get_codons(dna)
    rest = ''

    for i in codons:
        if i in stop_codons:
            break
        else:
            rest += i
    return rest


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    codons = get_codons(dna)
    ORFs = []
    i = 0
    while i < len(codons):
    # for i in range(len(codons)): # I would much rather use a for loop than a while loop...
	    if codons[i] == 'ATG':
	    	ORFs.append(rest_of_ORF(dna[i*3:]))
	    	i += len(rest_of_ORF(dna))/3
	    i += 1
    return ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    return find_all_ORFs_oneframe(dna) + find_all_ORFs_oneframe(dna[1:]) + find_all_ORFs_oneframe(dna[2:])


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    
    return find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    # return find_all_ORFs_both_strands(dna).sort() # HAHA NOPE
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    # doctest.run_docstring_examples(find_all_ORFs, globals())
