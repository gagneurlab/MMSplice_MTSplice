import numpy as np
import sys
import pdb
import warnings
import pandas as pd

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


class Exon(SpliceSite):

    def __init__(self,
                 *args,
                 **kwargs):
        ''' Exon with flanking intron
        Args:
            order: order of splice site (donor or acceptor) in a transcript counted from 5' to 3'. 
        '''
        super().__init__(*args, **kwargs)

    def get_seq(self, fasta):
        seq = fasta.get_seq(self.chrom,
                            self.grange,
                            self.strand)
        return seq.upper()

    def get_mut_seq(self, fasta):
        if self.variant is None:
            warnings.warn("Variant not known, return empty string", UserWarning)
            return ""
        else:
            assert self.variant.side in (None, "left", "right")

            # check reference
            ref_check = fasta.get_seq(self.chrom,
                                      (self.variant.POS, self.variant.POS + len(self.variant.REF) - 1),  # -1 because 1 based
                                      self.strand)
            if ref_check != self.variant.REF:
                warnings.warn("Reference not match, cannot mutate, return original sequence.", UserWarning)
                print(self.variant)
                return self.get_seq(fasta)

            seq = self.get_seq(fasta)
            p = self._var_pos()

            ###########
            # seq[p] is variant for snp
            ############
            # if self.variant.POS == 91731667:
            # pdb.set_trace()

            if p < 0 or p >= len(seq):  # p is 1 based, len(seq) 0 based
                return seq
            elif self.variant.len_diff == 0:
                # SNP
                mut_seq = seq[:max(0, p)] + self.variant.ALT + seq[p + len(self.variant.REF):]
                assert len(mut_seq) == len(seq)
                return mut_seq
            elif self.variant.len_diff < 0:
                # insertion
                mut_seq = seq[:max(0, p)] + self.variant.ALT + seq[p + len(self.variant.REF):]
                seq = self.get_seq(fasta)
                assert len(seq) - len(mut_seq) == self.variant.len_diff
                if self.variant.side is None:
                    # only substitute
                    return mut_seq
                elif self.variant.side == 'left':
                    return mut_seq[-self.variant.len_diff:]
                else:
                    return mut_seq[:self.variant.len_diff]
            elif self.variant.len_diff > 0:
                # deletion
                if self.variant.side is None:
                    mut_seq = seq[:max(0, p)] + self.variant.ALT + seq[p + self.variant.len_diff + 1:]
                    assert len(seq) - len(mut_seq) == self.variant.len_diff
                    return mut_seq
                elif self.variant.side == 'left':
                    if self.strand == "+":
                        mut_seq = fasta.get_seq(self.chrom,
                                                (self.grange[0] - self.variant.len_diff, self.grange[1]),
                                                self.strand)
                    else:
                        mut_seq = fasta.get_seq(self.chrom,
                                                (self.grange[0], self.grange[1] + self.variant.len_diff),
                                                self.strand)
                    mut_seq = mut_seq[:p + self.variant.len_diff] + self.variant.ALT + seq[p + self.variant.len_diff + 1:]
                    assert len(mut_seq) == len(seq)
                    return mut_seq
                else:
                    if self.strand == "+":
                        mut_seq = fasta.get_seq(self.chrom,
                                                (self.grange[0], self.grange[1] + self.variant.len_diff),
                                                self.strand)
                    else:
                        mut_seq = fasta.get_seq(self.chrom,
                                                (self.grange[0] - self.variant.len_diff, self.grange[1]),
                                                self.strand)
                    mut_seq = seq[:p] + self.variant.ALT + mut_seq[p + self.variant.len_diff + 1:]
                    assert len(mut_seq) == len(seq)
                    return mut_seq

    def _var_pos(self):
        """ Get variant relative position in the sequence
        """
        if self.strand == "+":
            return self.variant.POS - self.grange[0]
        else:
            return self.grange[1] - self.variant.POS - len(self.variant.REF) + 1


def mutate_seq(seq, valrel, base, is_rc):
    """ valrel: variant relative position on the sequence
    base: base to substitute
    is_rc: is reverse complement
    """
    pass


class Target(object):
    """ Read (miso) target file, counts or PSI
    """

    def __init__(self,
                 target_file,
                 label_col='event_name',
                 iscounts=True):
        self.label_col = label_col
        self._read_target(target_file)

    def _read_target(self, target_file):
        dt = pd.read_csv(target_file, index_col=0)
        event_names = dt[self.label_col].tolist()
        self._index = event_names
        dt = dt.drop(self.label_col, axis=1)
        tissues = dt.columns
        dt = dt.as_matrix()
        dt = np.stack((dt, 1 - dt), axis=2)  # might bug if only one tissue
        self.target = dt
        self.tissues = tissues

    def get_target(self, name):
        try:
            inx = self._index.index(name)
            return self.target[inx]
        except:
            dim = self.target.shape
            return nans((dim[1:]))


def nans(shape, dtype=float):
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a


class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = None

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout


def clip(x): return np.clip(x, 0.00001, 0.99999)
def logit(x):
    x = clip(x)
    return np.log(x) - np.log(1 - x)


def expit(x): return 1. / (1. + np.exp(-x))

def scale_factor(delta_score, ref_psi, delta_psi):
    # delta_score: predicted delta_score on logit scale
    # delta_psi: measured delta_psi
    delta_score_true = logit(ref_psi + delta_psi) - logit(ref_psi)  # measured delta logit
    # remove NAs
    keep = np.isfinite(delta_score_true)
    return np.dot(delta_score[keep], delta_score_true[keep]) / np.dot(delta_score[keep], delta_score[keep])




