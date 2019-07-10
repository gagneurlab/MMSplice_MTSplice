import logging
from itertools import islice
from pkg_resources import resource_filename

import pandas as pd
import pyranges
from pybedtools import Interval
from kipoi.data import SampleIterator
from concise.preprocessing import encodeDNA
from kipoiseq.extractors import VariantSeqExtractor, MultiSampleVCF
from mmsplice.utils import pyrange_remove_chr_from_chrom_annotation


logger = logging.getLogger('mmsplice')
logger.addHandler(logging.NullHandler())

GRCH37 = resource_filename('mmsplice', 'models/grch37_exons.csv.gz')
GRCH38 = resource_filename('mmsplice', 'models/grch38_exons.csv.gz')


def read_exon_pyranges(gtf_file, overhang=(100, 100), first_last=True):
    '''
    Read exon as pyranges from gtf_file

    Args:
      gtf_file: gtf file from ensembl/gencode.
      overhang: padding of exon to match variants.
      first_last: set overhang of first and last exon of the gene to zero
        so seq intergenic region will not be processed.
    '''
    df_gtf = pyranges.read_gtf(gtf_file).df
    df_exons = df_gtf[df_gtf['Feature'] == 'exon']
    df_exons = df_exons[['Chromosome', 'Start', 'End', 'Strand',
                         'exon_id', 'gene_id', 'gene_name', 'transcript_id']]

    if first_last:
        df_genes = df_gtf[df_gtf['Feature'] == 'transcript']
        df_genes.set_index('transcript_id', inplace=True)
        df_genes = df_genes.loc[df_exons['transcript_id']]
        df_genes.set_index(df_exons.index, inplace=True)

        starting = df_exons['Start'] == df_genes['Start']
        ending = df_exons['End'] == df_genes['End']

        df_exons.loc[:, 'left_overhang'] = ~starting * 100
        df_exons.loc[:, 'right_overhang'] = ~ending * 100

        df_exons.loc[:, 'Start'] -= df_exons['left_overhang']
        df_exons.loc[:, 'End'] += df_exons['right_overhang']

    return pyranges.PyRanges(df_exons)


def batch_iter_vcf(vcf_file, batch_size=10000):
    '''
    Iterates variatns in vcf file.

    Args:
      vcf_file: path of vcf file.
      batch_size: size of each batch.
    '''
    variants = MultiSampleVCF(vcf_file)
    batch = list(islice(variants, batch_size))

    while batch:
        yield batch
        batch = list(islice(variants, batch_size))


def variants_to_pyranges(variants):
    '''
    Create pyrange object given list of variant objects.

    Args:
      variants: list of variant objects have CHROM, POS, REF, ALT properties.
    '''
    def _from_variants(variants):
        for v in variants:
            if len(v.ALT) == 1:
                yield v.CHROM, v.POS, v.POS + max(len(v.REF), len(v.ALT[0])), v
            else:
                # Only support one alternative.
                # If multiple alternative, need to split into multiple variants
                logger.warning(
                    '%s has more than one or nan ALT sequence,'
                    'split into mutliple variants with bedtools' % v)

    df = pd.DataFrame(list(_from_variants(variants)),
                      columns=['Chromosome', 'Start', 'End', 'variant'])
    return pyranges.PyRanges(df)


def read_vcf_pyranges(vcf_file, batch_size=10000):
    '''
    Reads vcf and returns batch of pyranges objects.

    Args:
      vcf_file: path of vcf file.
      batch_size: size of each batch.
    '''
    for batch in batch_iter_vcf(vcf_file, batch_size):
        yield variants_to_pyranges(batch)


class ExonSeqVcfSeqExtrator:
    """
    Extracts variant applied sequence of exon with give overhang. If variant
    is in exon extractor extends exon but if variant is in intron it keeps
    intron (overhang) length fixed.
    """

    def __init__(self, fasta_file, overhang=(100, 100)):
        self.variant_seq_extractor = VariantSeqExtractor(fasta_file)
        self.fasta = self.variant_seq_extractor.fasta
        self.overhang = overhang

    def extract(self, interval, variants, sample_id=None):
        """
        Args:
          interval (pybedtools.Interval): zero-based interval of exon
            without overhang.
        """
        down_interval = Interval(
            interval.chrom, interval.start - self.overhang[0],
            interval.start, strand=interval.strand)
        up_interval = Interval(
            interval.chrom, interval.end,
            interval.end + self.overhang[1], strand=interval.strand)

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
                 donor_exon_len=5, donor_intron_len=13):
        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.acceptor_intron_cut = acceptor_intron_cut
        self.donor_intron_cut = donor_intron_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_exon_len = donor_exon_len
        self.donor_intron_len = donor_intron_len

    def split(self, x, overhang, exon_row='', pattern_warning=True):
        """
        Split seqeunce for each module.

        Args:
          seq: seqeunce to split.
          overhang: (acceptor, donor) overhang.
        """
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


class SplicingVCFDataloader(SampleIterator):
    """
    Load genome annotation (gtf) file along with a vcf file,
      return wt sequence and mut sequence.

    Args:
        gtf: gtf file. Can be dowloaded from ensembl/gencode.
          Filter for protein coding genes.
        fasta_file: file path; Genome sequence
        vcf_file: vcf file, each line should contain one
          and only one variant, left-normalized
        spit_seq: whether or not already split the sequence
          when loading the data. Otherwise it can be done in the model class.
        endcode: if split sequence, should it be one-hot-encoded
    """

    def __init__(self, gtf, fasta_file, vcf_file,
                 variant_filter=True, split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None):

        self.gtf_file = gtf
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.split_seq = split_seq
        self.encode = encode
        self.spliter = seq_spliter or SeqSpliter()

        self.pr_exons = self._read_exons(gtf)
        self.vseq_extractor = ExonSeqVcfSeqExtrator(fasta_file, overhang)
        self.fasta = self.vseq_extractor.fasta
        self.variants_batchs = read_vcf_pyranges(vcf_file)
        self.vcf = MultiSampleVCF(vcf_file)

        self._check_chrom_annotation()

        self._generator = self._generate(variant_filter=variant_filter)

    def _check_chrom_annotation(self):
        fasta_chroms = set(self.fasta.fasta.keys())
        vcf_chroms = set(self.vcf.seqnames)

        if not fasta_chroms.intersection(vcf_chroms):
            raise ValueError(
                'Fasta chrom names do not match with vcf chrom names')

        if self.gtf_file == 'grch37' or self.gtf_file == 'grch38':
            chr_annotaion = any(chrom.startswith('chr')
                                for chrom in vcf_chroms)
            if not chr_annotaion:
                self.pr_exons = pyrange_remove_chr_from_chrom_annotation(
                    self.pr_exons)

        gtf_chroms = set(self.pr_exons.Chromosome)
        if not gtf_chroms.intersection(vcf_chroms):
            raise ValueError(
                'GTF chrom names do not match with vcf chrom names')

    def _read_exons(self, gtf):
        if gtf == 'grch37':
            return pyranges.PyRanges(pd.read_csv(GRCH37))
        elif gtf == 'grch38':
            return pyranges.PyRanges(pd.read_csv(GRCH38))
        else:
            return read_exon_pyranges(self.gtf_file)

    def _generate(self, variant_filter=True):
        for pr_variants in self.variants_batchs:

            exon_variant_pairs = pr_variants.join(
                self.pr_exons, suffix="_exon")

            for i, row in exon_variant_pairs.df.iterrows():
                yield row

    def __iter__(self):
        return self

    def __next__(self):
        row = next(self._generator)
        overhang = (row['left_overhang'], row['right_overhang'])
        exon = Interval(row['Chromosome'],
                        row['Start_exon'] + overhang[0] - 1,
                        row['End_exon'] - overhang[1],
                        strand=row['Strand'])
        variant = row['variant']

        seq = self.fasta.extract(Interval(
            exon.chrom, exon.start - overhang[0],
            exon.end + overhang[1], strand=exon.strand))
        mut_seq = self.vseq_extractor.extract(exon, [variant])

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

    def batch_iter(self, batch_size=32):
        encode = self.encode
        self.encode = False

        for batch in super().batch_iter(batch_size):
            if encode:
                batch['inputs']['seq'] = self._encode_batch_seq(
                    batch['inputs']['seq'])
                batch['inputs']['mut_seq'] = self._encode_batch_seq(
                    batch['inputs']['mut_seq'])

            yield batch

        self.encode = encode

    def _encode_seq(self, seq):
        return {k: encodeDNA([v.upper()]) for k, v in seq.items()}

    def _encode_batch_seq(self, batch):
        return {k: encodeDNA([i.upper() for i in v.tolist()])
                for k, v in batch.items()}

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
