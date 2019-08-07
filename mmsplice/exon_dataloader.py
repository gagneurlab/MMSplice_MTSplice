import logging
import pandas as pd
from pybedtools import Interval
from concise.preprocessing import encodeDNA
from kipoi.data import Dataset
from kipoiseq.extractors import VariantSeqExtractor
from mmsplice.utils import Variant

logger = logging.getLogger('mmsplice')
logger.addHandler(logging.NullHandler())


class ExonSeqVcfSeqExtrator:
    """
    Extracts variant applied sequence of exon with give overhang. If variant
    is in exon extractor extends exon but if variant is in intron it keeps
    intron (overhang) length fixed.
    """

    def __init__(self, fasta_file):
        self.variant_seq_extractor = VariantSeqExtractor(fasta_file)
        self.fasta = self.variant_seq_extractor.fasta

    def extract(self, interval, variants, sample_id=None, overhang=(100, 100)):
        """
        Args:
          interval (pybedtools.Interval): zero-based interval of exon
            without overhang.
        """
        down_interval = Interval(
            interval.chrom, interval.start - overhang[0],
            interval.start, strand=interval.strand)
        up_interval = Interval(
            interval.chrom, interval.end,
            interval.end + overhang[1], strand=interval.strand)

        down_seq = self.variant_seq_extractor.extract(
            down_interval, variants, anchor=interval.start)
        up_seq = self.variant_seq_extractor.extract(
            up_interval, variants, anchor=interval.start)

        exon_seq = self.variant_seq_extractor.extract(
            interval, variants, anchor=0, fixed_len=False)

        if interval.strand == '-':
            down_seq, up_seq = up_seq, down_seq

        return down_seq + exon_seq + up_seq


class SeqSpliter:
    """
    Splits given seq for each modules.

    Args:
      exon_cut_l: number of bp to cut out at the begining of an exon
      exon_cut_r: number of bp to cut out at the end of an exon
        (cut out the part that is considered as acceptor site or donor site)
      acceptor_intron_cut: number of bp to cut out at the end of
        acceptor intron that consider as acceptor site
      donor_intron_cut: number of bp to cut out at the end of donor intron
        that consider as donor site
      acceptor_intron_len: length in acceptor intron to consider
        for acceptor site model
      acceptor_exon_len: length in acceptor exon to consider
        for acceptor site model
      donor_intron_len: length in donor intron to consider for donor site model
      donor_exon_len: length in donor exon to consider for donor site model
    """

    def __init__(self, exon_cut_l=0, exon_cut_r=0,
                 acceptor_intron_cut=6, donor_intron_cut=6,
                 acceptor_intron_len=50, acceptor_exon_len=3,
                 donor_exon_len=5, donor_intron_len=13, pattern_warning=True):
        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.acceptor_intron_cut = acceptor_intron_cut
        self.donor_intron_cut = donor_intron_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_exon_len = donor_exon_len
        self.donor_intron_len = donor_intron_len
        self.pattern_warning = pattern_warning

    def split(self, x, overhang, exon_row='', pattern_warning=True):
        """
        Split seqeunce for each module.

        Args:
          seq: seqeunce to split.
          overhang: (acceptor, donor) overhang.
        """
        pattern_warning = self.pattern_warning and pattern_warning

        intronl_len, intronr_len = overhang
        # need to pad N if left seq not enough long
        lackl = self.acceptor_intron_len - intronl_len
        if lackl >= 0:
            x = "N" * (lackl + 1) + x
            intronl_len += lackl + 1
        lackr = self.donor_intron_len - intronr_len
        if lackr >= 0:
            x = x + "N" * (lackr + 1)
            intronr_len += lackr + 1

        acceptor_intron = x[:intronl_len - self.acceptor_intron_cut]

        acceptor_start = intronl_len - self.acceptor_intron_len
        acceptor_end = intronl_len + self.acceptor_exon_len
        acceptor = x[acceptor_start: acceptor_end]

        exon_start = intronl_len + self.exon_cut_l
        exon_end = -intronr_len - self.exon_cut_r
        exon = x[exon_start: exon_end]

        donor_start = -intronr_len - self.donor_exon_len
        donor_end = -intronr_len + self.donor_intron_len
        donor = x[donor_start: donor_end]

        donor_intron = x[-intronr_len + self.donor_intron_cut:]

        if not exon:
            exon = 'N'

        if pattern_warning:
            if donor[self.donor_exon_len:self.donor_exon_len + 2] != "GT" \
               and overhang[1]:
                logger.warning('None GT donor: %s' % str(exon_row))

            if acceptor[self.acceptor_intron_len - 2:self.acceptor_intron_len] != "AG" \
               and overhang[0]:
                logger.warning('None AG acceptor: %s' % str(exon_row))

        splits = {
            "acceptor_intron": acceptor_intron,
            "acceptor": acceptor,
            "exon": exon,
            "donor": donor,
            "donor_intron": donor_intron
        }

        return splits


class BaseExonSplicingDataloader:

    def __init__(self, fasta_file, split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None):
        self.fasta_file = fasta_file
        self.split_seq = split_seq
        self.encode = encode
        self.spliter = seq_spliter or SeqSpliter()
        self.vseq_extractor = ExonSeqVcfSeqExtrator(fasta_file)
        self.fasta = self.vseq_extractor.fasta

    def _variant_to_dict(self, variant):
        return {
            'CHROM': variant.CHROM,
            'POS': variant.POS,
            'REF': variant.REF,
            'ALT': variant.ALT,
            'STR': "%s:%s:%s:['%s']" % (variant.CHROM, str(variant.POS),
                                        variant.REF, variant.ALT[0])
        }

    def _exon_to_dict(self, row, exon, overhang):
        return {
            'chrom': exon.chrom,
            'start': exon.start,
            'end': exon.end,
            'strand': exon.strand,
            'exon_id': row['exon_id'],
            'gene_id': row['gene_id'],
            'gene_name': row['gene_name'],
            'transcript_id': row['transcript_id'],
            'left_overhang': overhang[0],
            'right_overhang': overhang[1],
            'annotation': '%s:%d-%d:%s' % (exon.chrom, exon.start,
                                           exon.end, exon.strand)
        }

    def _encode_seq(self, seq):
        return {k: encodeDNA([v]) for k, v in seq.items()}

    def _next(self, row, exon, variant, overhang):
        seq = self.fasta.extract(Interval(
            exon.chrom, exon.start - overhang[0],
            exon.end + overhang[1], strand=exon.strand)).upper()
        mut_seq = self.vseq_extractor.extract(
            exon, [variant], overhang=overhang).upper()

        if exon.strand == '-':
            overhang = (overhang[1], overhang[0])

        if self.split_seq:
            seq = self.spliter.split(seq, overhang, exon)
            mut_seq = self.spliter.split(mut_seq, overhang, exon,
                                         pattern_warning=False)
            if self.encode:
                seq = self._encode_seq(seq)
                mut_seq = self._encode_seq(mut_seq)

        return {
            'inputs': {
                'seq': seq,
                'mut_seq': mut_seq
            },
            'metadata': {
                'variant': self._variant_to_dict(variant),
                'exon': self._exon_to_dict(row, exon, overhang)
            }
        }


class ExonDataset(Dataset):

    def __init__(self, exon_file, fasta_file, split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None, exon_file_kwargs=None):
        """
        """
        self.exon_file = exon_file
        self.fasta_file = fasta_file
        self.split_seq = split_seq
        self.encode = encode
        self.spliter = seq_spliter or SeqSpliter()

        self.exons = self._read_exon_file(exon_file, **exon_file_kwargs)
        self.vseq_extractor = ExonSeqVcfSeqExtrator(fasta_file)
        self.fasta = self.vseq_extractor.fasta

    def _read_exon_file(self, exon_file, **kwargs):
        return pd.read_csv(exon_file, **kwargs)

    def __getitem__(self, idx):
        row = self.exons.iloc[idx]
        exon = Interval(row['chrom'], row['start'], row['end'], row['strand'])
        variant = Variant(row['chrom'], row['pos'], row['ref'], [row['alt']])

    def __len__(self):
        pass
