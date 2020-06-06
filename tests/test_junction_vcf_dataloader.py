import pandas as pd
from kipoiseq import Interval
from conftest import junction_file, fasta_file, junction_psi5_file,  \
    junction_psi3_file
from mmsplice.junction_dataloader import JunctionVCFDataloader, \
    JunctionPSI5Dataset, JunctionPSI3Dataset
from mmsplice.utils import Variant
from mmsplice import MMSplice
from mmsplice import predict_all_table


def test_JunctionPSI5Dataset_junction_to_acceptor_exons():
    df_junc = JunctionPSI5Dataset(junction_psi5_file, fasta_file).exons
    assert df_junc.iloc[0].Exon_Start == 18159912
    assert df_junc.iloc[0].Exon_End == 18160011
    assert df_junc.iloc[1].Exon_Start == 18159740
    assert df_junc.iloc[1].Exon_End == 18159839


def test_JunctionPSI5Dataset__getitem__():
    dl = JunctionPSI5Dataset(junction_psi5_file, fasta_file, encode=False)
    d = dl[0]

    assert d['inputs']['seq'] == {
        'acceptor_intron': 'GCTCGGTGGCCTCCCGGCCAGCACTCACGCGCAGCGATGGCAGGCT'
        'GGACAGCTCCCCATGAAGCGTGGTCAGGTTGTTGTGGCTCACAGACAA',
        'acceptor': 'AGCTCCCCATGAAGCGTGGTCAGGTTGTTGTGGCTCACAGACAAGTGTTCCTG',
        'exon': 'CTGGGCCGGAATGGGGAGCACAGGCTCAGCAGGGTGGGTGGGGTGTCCAGAGTCCCT'
        'GCCTGCCTCCACAACTCTGCACCATCCAGGAAGGCCTCAGCCA',
        'donor': 'AGCCANNNNNNNNNNNNN',
        'donor_intron': 'NNNNNNNN'
    }
    assert d['inputs']['mut_seq'] == {
        'acceptor_intron': 'GCTCGGTGGCCTCCCGGCCAGCACTCACGCGCAGCGATGGCAGGCT'
        'GGACAGCTCCCCATGAAGCGTGGTCAGGTTGTTGTGGCTCACAGACAA',
        'acceptor': 'AGCTCCCCATGAAGCGTGGTCAGGTTGTTGTGGCTCACAGACAAGTGTTCCTG',
        'exon': 'CTGGGCCGGAATGGGGAGCACAGGCTCAGCAGGGTGGGTGGGGTGTCCAGAGTCCCT'
        'GCCTGCCTCCACAACTCTGCACCATCCAGGAAGGCCTCAGCCA',
        'donor': 'AGCCANNNNNNNNNNNNN',
        'donor_intron': 'NNNNNNNN'
    }
    assert d['metadata'] == {
        'variant': {
            'chrom': '17',
            'pos': 18159808,
            'ref': 'C',
            'alt': 'G',
            'annotation': '17:18159808:C>G'
        },
        'exon': {
            'chrom': '17',
            'start': 18159911,
            'end': 18160011,
            'strand': '+',
            'junction': '17:18159839-18159911:+',
            'left_overhang': 100,
            'right_overhang': 0,
            'annotation': '17:18159911-18160011:+'
        }
    }


def test_JunctionPSI3Dataset_junction_to_donor_exons():
    df_junc = JunctionPSI3Dataset(junction_psi3_file, fasta_file).exons
    assert df_junc.iloc[0].Exon_Start == 18159840 - 100
    assert df_junc.iloc[0].Exon_End == 18159840 - 1
    assert df_junc.iloc[1].Exon_Start == 18159911 + 1
    assert df_junc.iloc[1].Exon_End == 18159911 + 100


def test_JunctionPSI3Dataset__getitem__():
    dl = JunctionPSI3Dataset(junction_psi5_file, fasta_file, encode=False)
    d = dl[0]

    assert d['inputs']['seq'] == {
        'donor_intron': 'CGATGGCAGGCTGGACAGCTCCCCATGAAGCGTGGTCAGGTTGTTGTGG'
        'CTCACAGACAAGTGTTCCTGGGCCGGAATGGGGAGCACAGGCTCA',
        'donor': 'CTCACGCGCAGCGATGGC',
        'exon': 'CAGAAATGCTGGCCCAGCTTGGTGTCTAGCCCAGGCTTCACGTCCCTGACCGGACCT'
        'TGGCCCCAACCCCAAGCTCGGTGGCCTCCCGGCCAGCACTCAC',
        'acceptor': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCAG',
        'acceptor_intron': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    }
    assert d['inputs']['mut_seq'] == {
        'acceptor_intron': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN',
        'acceptor': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCAG',
        'exon': 'CAGAAATGCTGGCCCAGCTTGGTGTCTAGCCCAGGCTTCACGTCCCTGACCGGACCT'
        'TGGCCCCAACCGCAAGCTCGGTGGCCTCCCGGCCAGCACTCAC',
        'donor': 'CTCACGCGCAGCGATGGC',
        'donor_intron': 'CGATGGCAGGCTGGACAGCTCCCCATGAAGCGTGGTCAGGTTGTTGTGG'
        'CTCACAGACAAGTGTTCCTGGGCCGGAATGGGGAGCACAGGCTCA'
    }

    assert d['metadata'] == {
        'variant': {
            'annotation': '17:18159808:C>G',
            'ref': 'C',
            'pos': 18159808,
            'chrom': '17',
            'alt': 'G'
        },
        'exon': {
            'chrom': '17',
            'start': 18159739,
            'end': 18159839,
            'strand': '+',
            'left_overhang': 0,
            'right_overhang': 100,
            'annotation': '17:18159739-18159839:+',
            'junction': '17:18159839-18159911:+'
        }
    }


def test_JunctionVCFDataloader__init__(vcf_path):
    junc_dl = JunctionVCFDataloader(junction_file, fasta_file, vcf_path)
    assert junc_dl

    rows = list(junc_dl._generator)
    assert len(rows) != 0


def test_JunctionVCFDataloader_read_junction(vcf_path):
    df = pd.read_csv(junction_file)
    df_junc = JunctionVCFDataloader._read_junction(
        junction_file, overhang=(100, 100), exon_len=50).df

    assert df_junc.shape[0] == df.shape[0] * 2

    row = df_junc[(df_junc['side'] == 'donor') &
                  (df_junc['Strand'] == '-')].iloc[0]
    assert row['Start'] == 41279742 - 100
    assert row['End'] == 41279742 + 50

    row = df_junc[(df_junc['side'] == 'acceptor') &
                  (df_junc['Strand'] == '-')].iloc[0]
    assert row['Start'] == 41276032 - 50
    assert row['End'] == 41276032 + 100

    row = df_junc[(df_junc['side'] == 'donor') &
                  (df_junc['Strand'] == '+')].iloc[0]
    assert row['Start'] == 41276002 - 50
    assert row['End'] == 41276002 + 100

    row = df_junc[(df_junc['side'] == 'acceptor') &
                  (df_junc['Strand'] == '+')].iloc[0]
    assert row['Start'] == 41279042 - 100
    assert row['End'] == 41279042 + 50


def test_JunctionVCFDataloader__next__(vcf_path):
    junc_dl = JunctionVCFDataloader(
        junction_file, fasta_file, vcf_path,
        encode=False, split_seq=False)

    junc_dl._generator = iter([
        (
            Interval('17', 41279042 - 100, 41279042 + 100, strand='-',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 0,
                         'junction': '17:41276003-41279042:-',
                         'side': 'donor'
                     }),
            Variant('17', 41279035, 'C', 'G')
        )
    ])
    d = next(junc_dl)

    assert d['metadata']['exon'] == {
        'chrom': '17',
        'start': 41279042,
        'end': 41279042 + 100,
        'strand': '-',
        'right_overhang': 0,
        'left_overhang': 100,
        'annotation': '17:41279042-41279142:-',
        'junction': '17:41276003-41279042:-',
        'side': 'donor'
    }

    assert d['inputs']['seq'] == (
        'CACGCTGGAGCTCACAGGGGTTGGGAGGGGGCTCGGGCATGGCGGGCTGCAGGTCCTGAGCCTTGCCC'
        'TGTGCAGGGCGGCTGGGGCCCGGTGAGAATTCAAGCGGGGTGCAGGCGGGCCGGCAGTGCTGGGGGACC'
        'CGGCGCACCCTCTGCAGCTGCTGGCCCGGGTGCTAGGCCCCTGACTGCCCGGGGCCGGGGGTG')
    assert len(d['inputs']['seq']) == 200

    assert d['inputs']['mut_seq'] == (
        'CACGCTGGAGCTCACAGGGGTTGGGAGGGGGCTCGGGCATGGCGGGCTGCAGGTCCTGAGCCTTGCCC'
        'TGTGCAGGGCGGCTGGGGCCCGGTGAGAATTCAAGCGGGCTGCAGGCGGGCCGGCAGTGCTGGGGGACC'
        'CGGCGCACCCTCTGCAGCTGCTGGCCCGGGTGCTAGGCCCCTGACTGCCCGGGGCCGGGGGTG')
    assert len(d['inputs']['mut_seq']) == 200


def test_JunctionVCFDataloader__next__split(vcf_path):
    junc_dl = JunctionVCFDataloader(
        junction_file, fasta_file, vcf_path,
        encode=False, split_seq=True)

    junc_dl._generator = iter([
        (
            Interval('17', 41279042 - 100, 41279042 + 100, strand='-',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 0,
                         'junction': '17:41276003-41279042:-',
                         'side': 'donor'
                     }),
            Variant('17', 41279035, 'C', 'G')
        )
    ])
    d = next(junc_dl)
    assert d['inputs']['seq']['acceptor_intron'] == 'N' * 45
    assert d['inputs']['mut_seq']['acceptor_intron'] == 'N' * 45


def test_predict_all_table(vcf_path):
    model = MMSplice()
    dl = JunctionVCFDataloader(junction_file, fasta_file, vcf_path)
    df = predict_all_table(model, dl, pathogenicity=True,
                           splicing_efficiency=True)
    assert 'junction' in df.columns
    assert 'side' in df.columns
    assert len(df['delta_logit_psi']) > 0

    dl = JunctionPSI5Dataset(junction_psi5_file, fasta_file)
    df = predict_all_table(model, dl, pathogenicity=True,
                           splicing_efficiency=True)
    assert 'junction' in df.columns
    assert len(df['delta_logit_psi']) > 0

    dl = JunctionPSI3Dataset(junction_psi3_file, fasta_file)
    df = predict_all_table(model, dl, pathogenicity=True,
                           splicing_efficiency=True)
    assert 'junction' in df.columns
    assert len(df['delta_logit_psi']) > 0
