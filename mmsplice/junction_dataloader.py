import logging
import numpy as np
import pandas as pd
import pyranges
from kipoi.data import SampleIterator
from mmsplice.exon_dataloader import ExonDataset
from mmsplice.vcf_dataloader import SplicingVCFMixin

logger = logging.getLogger('mmsplice')


def junction_df_junction_str(df):
    return df['Chromosome'] + \
        ':' + df['Junction_Start'].astype('str') + \
        '-' + df['Junction_End'].astype('str') + \
        ':' + df['Strand']


class JunctionPSI5Dataset(ExonDataset):
    required_cols = ('Chromosome', 'Junction_Start', 'Junction_End',
                     'Strand', 'pos', 'ref', 'alt')

    def __init__(self, exon_file, fasta_file, split_seq=True, encode=True,
                 overhang=100, seq_spliter=None,
                 exon_len=100, **kwargs):
        super().__init__(exon_file, fasta_file, split_seq=split_seq,
                         encode=encode, overhang=(overhang, 0),
                         seq_spliter=seq_spliter, **kwargs)
        self.exon_len = exon_len
        self.exons = self._junction_to_acceptor_exons(self.exons, exon_len)

    @staticmethod
    def _junction_to_acceptor_exons(df, exon_len):
        # intron start-end position to exon start-end position (1-based)
        df['Junction_Start'] -= 1
        df['Junction_End'] += 1

        # calculates acceptor exon start end based
        # on given fixed exon lenght parameter and overhang (1-based)
        df['Exon_Start'] = np.where(df['Strand'] == '-',
                                    df['Junction_Start'] - exon_len + 1,
                                    df['Junction_End'])
        df['Exon_End'] = df['Exon_Start'] + exon_len - 1

        # convert junction position to 0-based
        df['Junction_End'] -= 1
        df['junction'] = junction_df_junction_str(df)
        return df.rename(columns=ExonDataset.exon_cols_mapping)


class JunctionPSI3Dataset(ExonDataset):
    required_cols = ('Chromosome', 'Junction_Start', 'Junction_End',
                     'Strand', 'pos', 'ref', 'alt')

    def __init__(self, exon_file, fasta_file, split_seq=True, encode=True,
                 overhang=100, seq_spliter=None,
                 exon_len=100, **kwargs):
        super().__init__(exon_file, fasta_file, split_seq=split_seq,
                         encode=encode, overhang=(0, overhang),
                         seq_spliter=seq_spliter, **kwargs)
        self.exon_len = exon_len
        self.exons = self._junction_to_donor_exons(self.exons, exon_len)

    @staticmethod
    def _junction_to_donor_exons(df, exon_len):
        # intron start-end position to exon start-end position (1-based)
        df['Junction_Start'] -= 1
        df['Junction_End'] += 1

        # calculates acceptor exon start end based
        # on given fixed exon lenght parameter and overhang (1-based)
        df['Exon_Start'] = np.where(df['Strand'] == '-',
                                    df['Junction_End'],
                                    df['Junction_Start'] - exon_len + 1)
        df['Exon_End'] = df['Exon_Start'] + exon_len - 1

        # convert junction position to 0-based
        df['Junction_End'] -= 1
        df['junction'] = junction_df_junction_str(df)
        return df


class JunctionVCFDataloader(SplicingVCFMixin, SampleIterator):
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
                 split_seq=True, encode=True,
                 overhang=(100, 100), seq_spliter=None, exon_len=100):
        pr_exons = self._read_junction(intron_annotation, overhang, exon_len)
        super().__init__(pr_exons, intron_annotation, fasta_file, vcf_file,
                         split_seq, encode, overhang, seq_spliter,
                         interval_attrs=('side', 'junction'))

    @staticmethod
    def _read_junction(intron_annotation, overhang=(100, 100), exon_len=100):
        df = pd.read_csv(intron_annotation, dtype={'Chromosome': str})

        # TODO: rearrange overhang
        # TO REFACTOR: refactor with junction_loader classes.
        df = pd.concat([df, JunctionVCFDataloader._junction_to_exon(
            df, overhang, exon_len)], axis=1)

        # split donor-acceptor from horizontal dataframe to two df
        # and merge this two df as one vertical dataframes
        df = df.rename(columns={"Start": "Junction_Start",
                                "End": "Junction_End"})
        df_donor = df[['Chromosome', 'Donor_Start', 'Donor_End',
                       'Strand', 'Junction_Start', 'Junction_End']] \
            .rename(columns={'Donor_Start': 'Start', 'Donor_End': 'End'})
        df_acceptor = df[['Chromosome', 'Acceptor_Start', 'Acceptor_End',
                          'Strand', 'Junction_Start', 'Junction_End']] \
            .rename(columns={'Acceptor_Start': 'Start', 'Acceptor_End': 'End'})
        df_donor['side'] = 'donor'
        df_acceptor['side'] = 'acceptor'
        df_exons = pd.concat([df_acceptor, df_donor])

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
        if (exon.attrs['side'] == 'donor' and exon.strand == '-') \
           or (exon.attrs['side'] == 'acceptor' and exon.strand == '+'):
            overhang = (self.overhang[0], 0)
        # ---*...---
        elif (exon.attrs['side'] == 'acceptor' and exon.strand == '-') \
                or (exon.attrs['side'] == 'donor' and exon.strand == '+'):
            overhang = (0, self.overhang[1])

        exon._start += overhang[0]
        exon._end -= overhang[1]
        return self._next(exon, variant, overhang)

    def __iter__(self):
        return self
