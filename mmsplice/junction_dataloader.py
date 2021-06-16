import logging
import numpy as np
import pandas as pd
import pyranges
from kipoi.data import SampleIterator
from mmsplice.exon_dataloader import ExonDataset
from mmsplice.vcf_dataloader import SplicingVCFMixin
from kipoiseq.dataclasses import Interval, Variant

logger = logging.getLogger('mmsplice')


def junction_df_junction_str(df):
    return df['Chromosome'] + \
        ':' + df['Junction_Start'].astype('str') + \
        '-' + df['Junction_End'].astype('str') + \
        ':' + df['Strand']


class _JunctionDataset(ExonDataset):
    required_cols = ('Chromosome', 'Junction_Start', 'Junction_End',
                     'Strand', 'pos', 'ref', 'alt')

    def __init__(self, exon_file, fasta_file, event_type, split_seq=True, encode=True,
                 overhang=100, seq_spliter=None,
                 exon_len=100, **kwargs):
        self.event_type = event_type
        super().__init__(exon_file, fasta_file, split_seq=split_seq,
                         encode=encode, overhang=self._get_overhang(overhang),
                         seq_spliter=seq_spliter, **kwargs)
        self.exon_len = exon_len
        self.exons = self._junction_to_acceptor_exons(
            self.exons, exon_len, event_type)

    def _get_overhang(self, overhang):
        if self.event_type == 'psi5':
            return (overhang, 0)
        elif self.event_type == 'psi3':
            return (0, overhang)
        else:
            raise ValueError('event_type should be "psi5" or "psi3"')

    @staticmethod
    def _junction_to_acceptor_exons(df, exon_len, event_type):
        # intron to 0-based position
        df['Junction_Start'] -= 1

        # calculates acceptor exon start end based
        if event_type == 'psi5':
            df['Exon_Start'] = np.where(df['Strand'] == '-',
                                        df['Junction_Start'] - exon_len,
                                        df['Junction_End'])
        elif event_type == 'psi3':
            df['Exon_Start'] = np.where(df['Strand'] == '-',
                                        df['Junction_End'],
                                        df['Junction_Start'] - exon_len)

        df['Exon_End'] = df['Exon_Start'] + exon_len
        df['junction'] = junction_df_junction_str(df)
        return df.rename(columns=ExonDataset.exon_cols_mapping)

    def __getitem__(self, idx):
        row = self.exons.iloc[idx]
        exon_attrs = {k: row[k] for k in self.optional_metadata if k in row}
        exon = Interval(row['Chromosome'], row['Exon_Start'],
                        row['Exon_End'], strand=row['Strand'],
                        attrs=exon_attrs)
        variant = Variant(row['Chromosome'], row['pos'],
                          row['ref'], row['alt'])

        overhang = self.overhang
        if exon.strand == '-':
            overhang = (overhang[1], overhang[0])

        if self.event_type == 'psi5':
            mask = ['donor', 'donor_intron']
        else:
            mask = ['acceptor', 'acceptor_intron']

        row = self._next(exon, variant, overhang, mask)
        return row


class JunctionPSI5Dataset(_JunctionDataset):

    def __init__(self, exon_file, fasta_file, split_seq=True, encode=True,
                 overhang=100, seq_spliter=None, exon_len=100, **kwargs):
        super().__init__(exon_file, fasta_file, 'psi5', split_seq=split_seq,
                         encode=encode, overhang=overhang,
                         seq_spliter=seq_spliter, exon_len=exon_len, **kwargs)


class JunctionPSI3Dataset(_JunctionDataset):

    def __init__(self, exon_file, fasta_file, split_seq=True, encode=True,
                 overhang=100, seq_spliter=None, exon_len=100, **kwargs):
        super().__init__(exon_file, fasta_file, 'psi3', split_seq=split_seq,
                         encode=encode, overhang=overhang,
                         seq_spliter=seq_spliter, exon_len=exon_len, **kwargs)


class _JunctionVCFDataloader(SplicingVCFMixin, SampleIterator):
    """
    Load intron annotation (bed) file along with a vcf file,
      return reference sequence and alternative sequence.

    Args:
      intron_annotation: 'grch37' or 'grch38' or
        path of tabular file contains intron (junction) annotation with
        colunms of `'Chromosome', 'Start', 'End', 'Strand'`. (0-based)
      fasta_file: file path; Genome sequence
      vcf_file: vcf file, each line should contain one
        and only one variant, left-normalized
      split_seq: whether or not already split the sequence
        when loading the data. Otherwise it can be done in the model class.
      endcode: if split sequence, should it be one-hot-encoded.
      overhang: overhang of exon to fetch flanking sequence of exon.
      seq_spliter: SeqSpliter class instance specific how to split seqs.
         if None, use the default arguments of SeqSpliter
    """

    def __init__(self, intron_annotation, fasta_file, vcf_file,
                 event_type, split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None, exon_len=100):
        self.event_type = event_type
        pr_exons = self._read_junction(intron_annotation, event_type,
                                       overhang, exon_len)
        super().__init__(pr_exons, intron_annotation, fasta_file, vcf_file,
                         split_seq, encode, overhang, seq_spliter,
                         interval_attrs=('junction',))

    @staticmethod
    def _read_junction(intron_annotation, event_type, overhang=(100, 100), exon_len=100):
        if type(intron_annotation) == str:
            df = pd.read_csv(intron_annotation, dtype={'Chromosome': str})
        else:
            df = intron_annotation
            df['Chromosome'] = df['Chromosome'].astype(str)

        # TODO: rearrange overhang
        # TO REFACTOR: refactor with junction_loader classes.
        df = pd.concat([df, _JunctionVCFDataloader._junction_to_exon(
            df, overhang, exon_len)], axis=1)

        # split donor-acceptor from horizontal dataframe to two df
        # and merge this two df as one vertical dataframes
        df = df.rename(columns={"Start": "Junction_Start",
                                "End": "Junction_End"})

        if event_type == 'psi5':
            df_exons = df[['Chromosome', 'Acceptor_Start', 'Acceptor_End',
                           'Strand', 'Junction_Start', 'Junction_End']] \
                .rename(columns={'Acceptor_Start': 'Start', 'Acceptor_End': 'End'})
        elif event_type == 'psi3':
            df_exons = df[['Chromosome', 'Donor_Start', 'Donor_End',
                           'Strand', 'Junction_Start', 'Junction_End']] \
                .rename(columns={'Donor_Start': 'Start', 'Donor_End': 'End'})
        else:
            raise ValueError('event_type should be "psi5" or "psi3"')

        df_exons['junction'] = junction_df_junction_str(df_exons)

        del df_exons['Junction_Start']
        del df_exons['Junction_End']
        return pyranges.PyRanges(df_exons)

    @staticmethod
    def _junction_to_exon(df, overhang=(100, 100), exon_len=100):
        # calculates donor-acceptor exon start end based
        # on given fixed exon lenght parameter and overhang (1-based)
        # --...-- represents junction (- exon), (. intron), (/ cut)
        mat = np.where(df['Strand'] == '-',
                       # --a./..d--/
                       (df['End'] - overhang[0],  df['End'] + exon_len,
                        # /--a../.d--
                        df['Start'] - exon_len, df['Start'] + overhang[1]),
                       # /--d../.a--
                       (df['Start'] - exon_len, df['Start'] + overhang[1],
                        # --d./..a--/
                        df['End'] - overhang[0], df['End'] + exon_len)
                       )
        return pd.DataFrame(
            mat.T,
            columns=['Donor_Start', 'Donor_End',
                     'Acceptor_Start', 'Acceptor_End'],
            index=df.index
        )

    def __next__(self):
        exon, variant = next(self._generator)

        # TODO: fix overhang will be reversed based on strand
        #   in `ExonSplicingMixin`
        # ---...*---
        if (self.event_type == 'psi3' and exon.strand == '-') \
           or (self.event_type == 'psi5' and exon.strand == '+'):
            overhang = (self.overhang[0], 0)
        # ---*...---
        elif (self.event_type == 'psi5' and exon.strand == '-') \
                or (self.event_type == 'psi3' and exon.strand == '+'):
            overhang = (0, self.overhang[1])

        exon._start += overhang[0]
        exon._end -= overhang[1]

        if self.event_type == 'psi5':
            mask = ['donor', 'donor_intron']
        else:
            mask = ['acceptor', 'acceptor_intron']

        row = self._next(exon, variant, overhang, mask)

        return row

    def __iter__(self):
        return self


class JunctionPSI5VCFDataloader(_JunctionVCFDataloader):
    def __init__(self, intron_annotation, fasta_file, vcf_file,
                 split_seq=True, encode=True, overhang=(100, 100),
                 seq_spliter=None, exon_len=100, **kwargs):
        super().__init__(intron_annotation, fasta_file, vcf_file,
                         'psi5', split_seq=split_seq, encode=encode,
                         overhang=overhang, seq_spliter=seq_spliter,
                         exon_len=exon_len, **kwargs)


class JunctionPSI3VCFDataloader(_JunctionVCFDataloader):
    def __init__(self, intron_annotation, fasta_file, vcf_file,
                 split_seq=True, encode=True, overhang=(100, 100),
                 seq_spliter=None, exon_len=100, **kwargs):
        super().__init__(intron_annotation, fasta_file, vcf_file, 'psi3',
                         split_seq=split_seq, encode=encode,
                         overhang=overhang, seq_spliter=seq_spliter,
                         exon_len=exon_len, **kwargs)
