import pandas as pd
from conftest import fasta_file, exon_file
from mmsplice.exon_dataloader import ExonDataset, ExonSplicingMixin
from kipoiseq.dataclasses import Interval, Variant


def test_ExonDataset():
    dl = ExonDataset(exon_file, fasta_file)
    required_cols = ('Chromosome', 'Exon_Start', 'Exon_End',
                     'Strand', 'pos', 'ref', 'alt')
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

    assert dl[10]['metadata']['variant']['chrom'] == '17'
    assert dl[10]['metadata']['variant']['pos'] == 18159891
    assert dl[10]['metadata']['variant']['ref'] == 'G'
    assert dl[10]['metadata']['variant']['alt'][0] == 'A'
    assert dl[10]['metadata']['variant']['annotation'] == "17:18159891:G>A"

    assert dl[10]['metadata']['exon']['chrom'] == '17'
    assert dl[10]['metadata']['exon']['start'] == 18159839
    assert dl[10]['metadata']['exon']['end'] == 18159911
    assert dl[10]['metadata']['exon']['strand'] == '-'
    assert dl[10]['metadata']['exon']['left_overhang'] == 100
    assert dl[10]['metadata']['exon']['right_overhang'] == 100
    assert dl[10]['metadata']['exon']['annotation'] == '17:18159839-18159911:-'


def test_ExonDataset_tissue():
    dl = ExonDataset(exon_file, fasta_file,
                     overhang=(300, 300),
                     encode=False,
                     split_seq=False,
                     tissue_specific=True)
    for x in dl:
        assert x['inputs']['mut_seq'] == x['inputs']['tissue_seq']


def test_ExonDataset__len__():
    dl = ExonDataset(exon_file, fasta_file)
    df = pd.read_csv(exon_file)
    assert len(dl) == df.shape[0]
    
    
def test_ExonSplicingMixin_extract_seq_in_mmsplice(vcf_path):
    import numpy as np
    from conftest import gtf_file, fasta_file
    from mmsplice.vcf_dataloader import SplicingVCFDataloader
    
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)

    rows = list(dl) 
    row = rows[0]
    assert 41203228 == row['metadata']['variant']['pos']

    mutated_seq = False   
    for module in ['acceptor_intron', 'acceptor', 'exon', 'donor', 'donor_intron']:
        # for snvs there should be only one mutated nucleotide -> sum over different entries in one-hot encoded vector == 2
        assert (pd.DataFrame(row['inputs'])['seq'][module] != pd.DataFrame(row['inputs'])['mut_seq'][module]).sum() <= 2
        if (pd.DataFrame(row['inputs'])['seq'][module] != pd.DataFrame(row['inputs'])['mut_seq'][module]).sum() == 2:
            mutated_seq = True

    assert mutated_seq == True