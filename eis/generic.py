import numpy as np
import sys
import warnings
import pandas as pd

clip = lambda x: np.clip(x, 0.00001, 0.99999)
def logit(x):
    x = clip(x)
    return np.log(x)-np.log(1-x)

def expit(x): return 1. / (1. + np.exp(-x))

def get_var_side(var):
    ''' Get exon variant side
    '''
    varstart, ref, alt, start, end, strand = var
    varend = varstart + len(ref) - 1
    # for insertion deletion, find the actual start of variant
    # e.g. A->AGG: start from G position, CA->CAGG:start from G position, 
    # ATT->A: start from T position, CAT->CA: start from T position
    # For SNP var.POS is the actual mutation position
    if len(ref) != len(alt):
    	# indels
    	varstart = varstart + min(len(ref), len(alt)) 
    if strand == "+":
        if varstart < start:
            return "left"
        elif varend > end:
            return "right"
        else:
            return None
    else:
        if varstart < start:
            return "right"
        elif varend > end:
            return "left"
        else:
            return None

bases = ['A', 'C', 'G', 'T']

def onehot(seq):
    X = np.zeros((len(seq), len(bases)))
    for i, char in enumerate(seq):
        if char == "N":
            pass
        else:
            X[i, bases.index(char.upper())] = 1
    return X

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])


class Variant(object):

    def __init__(self,
                 CHROM,
                 POS,
                 REF,
                 ALT,
                 ID=None,
                 strand="*",
                 side=None
                 ):
        """ side: variant side of reference, e.g. splice juction. 
        Sequence will be trimed or padded based on this for indels
        strand: specify if the variant is intepreted at any strand. e.g. in the case of splicing
        """
        self.CHROM = CHROM
        self.POS = POS
        self.ID = ID
        self.strand = strand
        if strand == "-":
            self.REF = self._reverse_complement(REF.upper())
            self.ALT = self._reverse_complement(ALT.upper())
        else:
            self.REF = REF.upper()
            self.ALT = ALT.upper()
        self._side = side
        self.len_diff = len(self.REF) - len(self.ALT)

    @property
    def is_snp(self):
        return len(self.REF) == 1 and len(self.ALT) == 1

    @property
    def side(self):
        return self._side

    @side.setter
    def side(self, value):
        self._side = value

    @property
    def to_dict(self):
        if self.strand == '-':
            # convert back
            REF = self._reverse_complement(self.REF)
            ALT = self._reverse_complement(self.ALT)
        else:
            REF = self.REF
            ALT = self.ALT
        return {'CHROM': self.CHROM,
               'POS': self.POS,
               'ID': self.ID,
               'REF': REF,
               'ALT': ALT,
               'STR':self.CHROM+ ":" + str(self.POS)+":" + REF + ":['" + ALT + "']"}

    def _reverse_complement(self, dna):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
        return ''.join([complement[base] for base in dna[::-1]])

    # def __str__(self):
    #     return "{0}:{1}, REF={2}, ALT={3}".format(self.CHROM, self.POS, self.REF, self.ALT)
    
    def __repr__(self):
        return "Variant(CHROM={0}, POS={1}, REF={2}, ALT={3}, ID={4})".format(self.CHROM, self.POS, self.REF, self.ALT, self.ID)



class SpliceSite(object):
    ''' A splice site with flanking intron and exon sequence
    Args:
    order: order of splice site (donor or acceptor) in a transcript counted from 5' to 3'. 
    variant: Variant class instance. Variant associated with this site.
    '''

    def __init__(self,
                 chrom,
                 start,
                 stop,
                 strand,
                 transcript_id,
                 gene_id,
                 biotype,
                 order=None,
                 variant=None, 
                 encode=True):
        self.chrom = chrom
        self.grange = (start, stop)
        self.strand = strand
        self.transcriptID = transcript_id
        self.geneID = gene_id
        self.biotype = biotype
        self.order = order
        self._seq = None
        self.variant = variant
        self.encode = encode

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, value):
        self._seq = value

    def get_seq(self, fasta):
        seq = fasta.get_seq(self.chrom,
                            self.grange,
                            self.strand)
        if self.encode:
            return onehot(seq.upper())
        else:
            return seq.upper()


