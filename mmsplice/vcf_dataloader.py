import logging
from pkg_resources import resource_filename
import pandas as pd
import pyranges
from kipoi.data import SampleIterator
from kipoiseq.extractors import MultiSampleVCF, SingleVariantMatcher
from mmsplice.utils import pyrange_remove_chr_from_chrom_annotation,\
    pyrange_add_chr_from_chrom_annotation
from mmsplice.exon_dataloader import ExonSplicingMixin

logger = logging.getLogger('mmsplice')

prebuild_annotation = {
    'grch37': resource_filename('mmsplice', 'models/grch37_exons.csv.gz'),
    'grch38': resource_filename('mmsplice', 'models/grch38_exons.csv.gz')
}


def read_exon_pyranges(gtf_file, overhang=(100, 100), first_last=True):
    '''
    Read exon as pyranges from gtf_file

    Args:
      gtf_file: gtf file from ensembl/gencode.
      overhang: padding of exon to match variants.
      first_last: set overhang of first and last exon of the gene to zero
        so seq intergenic region will not be processed.
    '''

    def _exon_filter(df, overhang=overhang, first_last=first_last):
        df_exons = df[df['Feature'] == 'exon']
        df_exons = df_exons[['Chromosome', 'Start', 'End', 'Strand',
                            'exon_id', 'gene_id', 'gene_name', 'transcript_id']]

        if first_last:
            df_genes = df[df['Feature'] == 'transcript']
            df_genes.set_index('transcript_id', inplace=True)
            df_genes = df_genes.loc[df_exons['transcript_id']]
            df_genes.set_index(df_exons.index, inplace=True)

            starting = df_exons['Start'] == df_genes['Start']
            ending = df_exons['End'] == df_genes['End']

            df_exons.loc[:, 'left_overhang'] = ~starting * overhang[0]
            df_exons.loc[:, 'right_overhang'] = ~ending * overhang[1]

            df_exons.loc[:, 'Start'] -= df_exons['left_overhang']
            df_exons.loc[:, 'End'] += df_exons['right_overhang']
            df_exons['Start'] = df_exons['Start'].astype('int32')
            df_exons['End'] = df_exons['End'].astype('int32')

        return df_exons

    df_gtf = pyranges.read_gtf(gtf_file)
    df_exons = df_gtf.apply(_exon_filter)

    return df_exons


class SplicingVCFMixin(ExonSplicingMixin):

    def __init__(self, pr_exons, annotation, fasta_file, vcf_file,
                 split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None,
                 tissue_specific=False, tissue_overhang=(300, 300),
                 interval_attrs=tuple()):
        super().__init__(fasta_file, split_seq, encode, overhang, seq_spliter,
                         tissue_specific, tissue_overhang)
        self.pr_exons = pr_exons
        self.annotation = annotation
        self.vcf_file = vcf_file
        self.vcf = MultiSampleVCF(vcf_file)
        self._check_chrom_annotation()
        self.matcher = SingleVariantMatcher(
            vcf_file, pranges=self.pr_exons,
            interval_attrs=interval_attrs
        )
        self._generator = iter(self.matcher)

    def _check_chrom_annotation(self):
        fasta_chroms = set(self.fasta.fasta.keys())
        vcf_chroms = set(self.vcf.seqnames)

        if not fasta_chroms.intersection(vcf_chroms):
            raise ValueError(
                'Fasta chrom names do not match with vcf chrom names')

        gtf_chroms = set(self.pr_exons.Chromosome)
        if not gtf_chroms.intersection(vcf_chroms):
            chr_annotaion = any(chrom.startswith('chr')
                                for chrom in vcf_chroms)
            if not chr_annotaion:
                self.pr_exons = pyrange_remove_chr_from_chrom_annotation(
                    self.pr_exons)
            else:
                self.pr_exons = pyrange_add_chr_from_chrom_annotation(
                    self.pr_exons)

        gtf_chroms = set(self.pr_exons.Chromosome)
        if not gtf_chroms.intersection(vcf_chroms):
            raise ValueError(
                'GTF chrom names do not match with vcf chrom names')


class SplicingVCFDataloader(SplicingVCFMixin, SampleIterator):
    """
    Load genome annotation (gtf) file along with a vcf file,
      return reference sequence and alternative sequence.

    Args:
      gtf: gtf file. Can be dowloaded from ensembl/gencode.
        Filter for protein coding genes.
      fasta_file: file path; Genome sequence
      vcf_file: vcf file, each line should contain one
        and only one variant, left-normalized
      split_seq: whether or not already split the sequence
        when loading the data. Otherwise it can be done in the model class.
      endcode: if split sequence, should it be one-hot-encoded.
      overhang: overhang of exon to fetch flanking sequence of exon.
      seq_spliter: SeqSpliter class instance specific how to split seqs.
         if None, use the default arguments of SeqSpliter
      tissue_specific: tissue specific predicts
      tissue_overhang: overhang of exon to fetch flanking sequence of
        tissue specific model.
    """

    def __init__(self, gtf, fasta_file, vcf_file,
                 split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None,
                 tissue_specific=False, tissue_overhang=(300, 300)):
        pr_exons = self._read_exons(gtf, overhang)
        super().__init__(pr_exons, gtf, fasta_file, vcf_file,
                         split_seq, encode, overhang, seq_spliter,
                         tissue_specific, tissue_overhang,
                         interval_attrs=('left_overhang', 'right_overhang',
                                         'exon_id', 'gene_id',
                                         'gene_name', 'transcript_id'))

    def _read_exons(self, gtf, overhang=(100, 100)):
        if gtf in prebuild_annotation:
            if overhang != (100, 100):
                logger.warning('Overhang argument will be ignored'
                               ' for prebuild annotation.')
            df = pd.read_csv(prebuild_annotation[gtf])
            df['Start'] -= 1  # convert prebuild annotation to 0-based
            return pyranges.PyRanges(df)
        else:
            return read_exon_pyranges(gtf, overhang=overhang)

    def __next__(self):
        exon, variant = next(self._generator)
        overhang = (exon.attrs['left_overhang'], exon.attrs['right_overhang'])
        exon._start += overhang[0]
        exon._end -= overhang[1]
        return self._next(exon, variant, overhang)

    def __iter__(self):
        return self
