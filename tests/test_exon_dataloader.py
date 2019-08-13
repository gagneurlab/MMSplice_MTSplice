import pandas as pd
from conftest import fasta_file, exon_file
from mmsplice.exon_dataloader import ExonDataset


def test_ExonDataset():
    dl = ExonDataset(exon_file, fasta_file)
    required_cols = ('CHROM', 'Exon_Start', 'Exon_End',
                     'strand', 'POS', 'REF', 'ALT')
    assert all(c in dl.exons.columns for c in required_cols)


def test_ExonDataset__getitem__():
    dl = ExonDataset(exon_file, fasta_file, encode=False, split_seq=False)
    for i in dl:
        assert i

    seq = (
        'TGGCTGAGGCCTTCCTGGATGGTGCAGAGTTGTGGAGGCAGGCAGGGACTCTGGACACCCCACCCAC'
        'CCTGCTGAGCCTGTGCTCCCCATTCCGGCCCAGGAACACTTGTCTGTGAGCCACAACAACCTGACCA'
        'CGCTTCATGGGGAGCTGTCCAGCCTGCCATCGCTGCGCGTGAGTGCTGGCCGGGAGGCCACCGAGCT'
        'TGGGGTTGGGGCCAAGGTCCGGTCAGGGACGTGAAGCCTGGGCTAGACACCAAGCTGGGCCAGCATT'
        'TCTG'
    )
    assert dl[10]['inputs']['seq'] == seq

    mut_seq = (
        'TGGCTGAGGCCTTCCTGGATGGTGCAGAGTTGTGGAGGCAGGCAGGGACTCTGGACACCCCACCCAC'
        'CCTGCTGAGCCTGTGCTCCCCATTCCGGCCCAGGAACACTTGTCTGTGAGCCATAACAACCTGACCA'
        'CGCTTCATGGGGAGCTGTCCAGCCTGCCATCGCTGCGCGTGAGTGCTGGCCGGGAGGCCACCGAGCT'
        'TGGGGTTGGGGCCAAGGTCCGGTCAGGGACGTGAAGCCTGGGCTAGACACCAAGCTGGGCCAGCATT'
        'TCTG'
    )
    assert dl[10]['inputs']['mut_seq'] == mut_seq

    assert dl[10]['metadata']['variant']['CHROM'] == '17'
    assert dl[10]['metadata']['variant']['POS'] == 18159891
    assert dl[10]['metadata']['variant']['REF'] == 'G'
    assert dl[10]['metadata']['variant']['ALT'][0] == 'A'
    assert dl[10]['metadata']['variant']['STR'] == "17:18159891:G:['A']"

    assert dl[10]['metadata']['exon']['chrom'] == '17'
    assert dl[10]['metadata']['exon']['start'] == 18159839
    assert dl[10]['metadata']['exon']['end'] == 18159911
    assert dl[10]['metadata']['exon']['strand'] == '-'
    assert dl[10]['metadata']['exon']['left_overhang'] == 100
    assert dl[10]['metadata']['exon']['right_overhang'] == 100
    assert dl[10]['metadata']['exon']['annotation'] == '17:18159839-18159911:-'


def test_ExonDataset__len__():
    dl = ExonDataset(exon_file, fasta_file)
    df = pd.read_csv(exon_file)
    assert len(dl) == df.shape[0]
