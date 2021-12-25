import numpy as np
from kipoiseq.dataclasses import Interval, Variant
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice.exon_dataloader import SeqSpliter
from conftest import gtf_file, fasta_file, variants, vcf_file


def test_SplicingVCFDataloader__check_chrom_annotation():
    dl = SplicingVCFDataloader('grch37', fasta_file, vcf_file)
    chroms = {str(i) for i in range(1, 22)}.union(['X', 'Y', 'M'])
    assert len(chroms.difference(set(dl.pr_exons.Chromosome))) == 0
    # assert sum(1 for i in dl) > 0


def test_SplicingVCFDataloader__encode_seq(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    encoded = dl._encode_seq({'acceptor': 'ATT'})
    np.testing.assert_array_equal(
        encoded['acceptor'][0],
        np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1.],
                  [0., 0., 0., 1.]])
    )


def test_SplicingVCFDataloader__encode_batch_seq(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    encoded = dl._encode_batch_seq({'acceptor': np.array(['ATT'])})
    np.testing.assert_array_equal(
        encoded['acceptor'][0],
        np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1.],
                  [0., 0., 0., 1.]])
    )


def test_SplicingVCFDataloader__read_exons(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    df_exons = dl._read_exons(gtf_file, overhang=(10, 20)).df
    row = df_exons[df_exons['exon_id'] == 'ENSE00003559512'].iloc[0]
    assert row['left_overhang'] == 10
    assert row['right_overhang'] == 20


def test_benchmark_SplicingVCFDataloader(benchmark, vcf_path):
    benchmark(SplicingVCFDataloader, gtf_file, fasta_file, vcf_path)


def test_splicing_vcf_dataloader_prebuild_grch37(vcf_path):
    dl = SplicingVCFDataloader('grch37', fasta_file, vcf_path)


def test_splicing_vcf_dataloader_prebuild_grch38(vcf_path):
    dl = SplicingVCFDataloader('grch38', fasta_file, vcf_path)


def test_splicing_vcf_loads_all(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    assert sum(1 for i in dl) == len(variants) - 1

    dl = SplicingVCFDataloader('grch37', fasta_file, vcf_file)
    assert sum(1 for i in dl) > 0


def test_SepSpliter_split_tissue_seq():
    spliter = SeqSpliter(tissue_acceptor_intron=9, tissue_acceptor_exon=3,
                         tissue_donor_intron=9, tissue_donor_exon=3)
    overhang = (9, 9)
    seq = 'ATCATCATC' + 'GGGAAA' + 'CGTGCTCGT'
    d = spliter.split_tissue_seq(seq, overhang)
    assert d['acceptor'] == 'ATCATCATC' + 'GGG'
    assert d['donor'] == 'AAA' + 'CGTGCTCGT'

    overhang = (7, 7)
    seq = 'CATCATC' + 'GGGAAA' + 'CGTGCTC'
    d = spliter.split_tissue_seq(seq, overhang)
    assert d['acceptor'] == 'NNCATCATC' + 'GGG'
    assert d['donor'] == 'AAA' + 'CGTGCTCNN'

    overhang = (11, 11)
    seq = 'ggATCATCATC' + 'GGGAAA' + 'CGTGCTCGTtt'
    d = spliter.split_tissue_seq(seq, overhang)
    assert d['acceptor'] == 'ATCATCATC' + 'GGG'
    assert d['donor'] == 'AAA' + 'CGTGCTCGT'

    overhang = (0, 0)
    seq = '' + 'GGGAAA' + ''
    d = spliter.split_tissue_seq(seq, overhang)
    assert d['acceptor'] == 'NNNNNNNNN' + 'GGG'
    assert d['donor'] == 'AAA' + 'NNNNNNNNN'


def test_SplicingVCFDataloader__next__(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False,
                               tissue_specific=True)
    dl._generator = iter([
        (
            Interval('17', 41275933, 41276232, strand='-',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 100,
                         'exon_id': 'exon_id',
                         'gene_id': 'gene_id',
                         'gene_name': 'gene_name',
                         'transcript_id': 'transcript_id',
                     }),
            Variant('17', 41276033, 'C', 'G')
        )
    ])

    expected_snps_seq = {
        'seq':
        'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTT'
        'TATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAG'
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGG'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'alt_seq':
        'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTT'
        'TATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAG'
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'tissue_seq': 'AGGCTACCACCACCTACCCGGTCAGTCACTCCTCTG'
        'TAGCTTTCTCTTTCTTGGAGAAAGGAAAAGACCCAAGGGGTTGGCAGCAA'
        'TATGTGAAAAAATTCAGAATTTATGTTGTCTAATTACAAAAAGCAACTTC'
        'TAGAATCTTTAAAAATAAAGGACGTTGTCATTAGTTCTTTGGTTTGTATT'
        'ATTCTAAAACCTTCCAAATCTTAAATTTACTTTATTTTAAAATGATAAAA'
        'TGAAGTTGTCATTTTATAAACCTTTTAAAAAGATATATATATATGTTTTT'
        'CTAATGTGTTAAAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCT'
        'TCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAG'
        'AGTGTCCCATCTGCTAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCT'
        'ATGATTATCTCCTATGCAAATGAACAGAATTGACCTTACATACTAGGGAA'
        'GAAAAGACATGTCTAGTAAGATTAGGCTATTGTAATTGCTGATTTTCTTA'
        'ACTGAAGAACTTTAAAAATATAGAAAATGATTCCTTGTTCTCCATCCACT'
        'CTGCCTCTCCCACTCCTCTCCTTTTCAACACAAATCCTGTGGTCCGGGAA'
        'AGACAGGGACTCTGTCTTGATTGGTTCTGCACTGGGGCAGGAATCTAGTT'
        'TAGATTAACTGGC'
    }

    d = next(dl)
    assert d['inputs']['seq'] == expected_snps_seq['seq']
    assert d['inputs']['mut_seq'] == expected_snps_seq['alt_seq']
    assert d['inputs']['tissue_seq'] == expected_snps_seq['tissue_seq']

    dl._generator = iter([
        (
            Interval('17', 41275933, 41276132, strand='-',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 0,
                         'exon_id': 'exon_id',
                         'gene_id': 'gene_id',
                         'gene_name': 'gene_name',
                         'transcript_id': 'transcript_id',
                     }),
            Variant('17', 41276033, 'C', 'G')
        )
    ])

    expected_snps_seq = {
        'seq':
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGG'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'alt_seq':
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'tissue_seq':
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTCT'
        'AGTAAGATTAGGCTATTGTAATTGCTGATTTTCTTAACTGAAGAACTTTA'
        'AAAATATAGAAAATGATTCCTTGTTCTCCATCCACTCTGCCTCTCCCACT'
        'CCTCTCCTTTTCAACACAAATCCTGTGGTCCGGGAAAGACAGGGACTCTG'
        'TCTTGATTGGTTCTGCACTGGGGCAGGAATCTAGTTTAGATTAACTGGC'
    }

    d = next(dl)
    assert d['inputs']['seq'] == expected_snps_seq['seq']
    assert d['inputs']['mut_seq'] == expected_snps_seq['alt_seq']
    assert d['inputs']['tissue_seq'] == expected_snps_seq['tissue_seq']


def test_SplicingVCFDataloader__next__split(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=True, encode=False,
                               tissue_specific=True)
    dl._generator = iter([
        (
            Interval('17', 41275933, 41276232, strand='-',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 100,
                         'exon_id': 'exon_id',
                         'gene_id': 'gene_id',
                         'gene_name': 'gene_name',
                         'transcript_id': 'transcript_id',
                     }),
            Variant('17', 41276033, 'C', 'G')
        )
    ])

    expected_snps_seq = {
        'seq':
        'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTT'
        'TATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAG'
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGG'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'alt_seq':
        'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTT'
        'TATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAG'
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'tissue_seq': {
            'acceptor': 'AGGCTACCACCACCTACCCGGTCAGTCACTCCTC'
            'TGTAGCTTTCTCTTTCTTGGAGAAAGGAAAAGACCCAAGGGGTTGG'
            'CAGCAATATGTGAAAAAATTCAGAATTTATGTTGTCTAATTACAAA'
            'AAGCAACTTCTAGAATCTTTAAAAATAAAGGACGTTGTCATTAGTT'
            'CTTTGGTTTGTATTATTCTAAAACCTTCCAAATCTTAAATTTACTT'
            'TATTTTAAAATGATAAAATGAAGTTGTCATTTTATAAACCTTTTAA'
            'AAAGATATATATATATGTTTTTCTAATGTGTTAAAGTTCATTGGAA'
            'CAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAA'
            'ATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC',
            'donor': 'GTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCT'
            'TCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATC'
            'TTAGAGTGTCCCATCTGCTAAGTCAGCACAAGAGTGTATTAATTTG'
            'GGATTCCTATGATTATCTCCTATGCAAATGAACAGAATTGACCTTA'
            'CATACTAGGGAAGAAAAGACATGTCTAGTAAGATTAGGCTATTGTA'
            'ATTGCTGATTTTCTTAACTGAAGAACTTTAAAAATATAGAAAATGA'
            'TTCCTTGTTCTCCATCCACTCTGCCTCTCCCACTCCTCTCCTTTTC'
            'AACACAAATCCTGTGGTCCGGGAAAGACAGGGACTCTGTCTTGATT'
            'GGTTCTGCACTGGGGCAGGAATCTAGTTTAGATTAACTGGC'
        }
    }

    d = next(dl)
    # TO TEST:
    # assert d['inputs']['seq'] == expected_snps_seq['seq']
    # assert d['inputs']['mut_seq'] == expected_snps_seq['alt_seq']
    assert d['inputs']['tissue_seq'] == expected_snps_seq['tissue_seq']

    dl._generator = iter([
        (
            Interval('17', 41275933, 41276132, strand='-',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 0,
                         'exon_id': 'exon_id',
                         'gene_id': 'gene_id',
                         'gene_name': 'gene_name',
                         'transcript_id': 'transcript_id',
                     }),
            Variant('17', 41276033, 'C', 'G')
        )
    ])

    expected_snps_seq = {
        'seq':
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGG'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'alt_seq':
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'tissue_seq': {
            'acceptor': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTCATTGGAA'
            'CAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAA'
            'ATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC',
            'donor': 'NTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCT'
            'TCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATC'
            'TTAGAGTGTCCCATCTGCTAAGTCAGCACAAGAGTGTATTAATTTG'
            'GGATTCCTATGATTATCTCCTATGCAAATGAACAGAATTGACCTTA'
            'CATACTAGGGAAGAAAAGACATGTCTAGTAAGATTAGGCTATTGTA'
            'ATTGCTGATTTTCTTAACTGAAGAACTTTAAAAATATAGAAAATGA'
            'TTCCTTGTTCTCCATCCACTCTGCCTCTCCCACTCCTCTCCTTTTC'
            'AACACAAATCCTGTGGTCCGGGAAAGACAGGGACTCTGTCTTGATT'
            'GGTTCTGCACTGGGGCAGGAATCTAGTTTAGATTAACTGGC'
        }
    }

    d = next(dl)
    assert d['inputs']['tissue_seq'] == expected_snps_seq['tissue_seq']
    assert len(d['inputs']['tissue_seq']['acceptor']) == 400
    assert len(d['inputs']['tissue_seq']['donor']) == 400


def test_SplicingVCFDataloader__next__split_seq_True(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=True, encode=False)
    dl._generator = iter([
        (
            Interval('17', 61278131 - 100, 61278318 + 100, strand='+',
                     attrs={
                         'left_overhang': 100,
                         'right_overhang': 100,
                         'exon_id': 'exon_id',
                         'gene_id': 'gene_id',
                         'gene_name': 'gene_name',
                         'transcript_id': 'transcript_id',
                     }),
            Variant('17', 61278319, 'G', 'A')
        )
    ])

    expected_snps_seq = {
        'seq': {
            'acceptor_intron':
            'GCATGTAGTTTTTTTTTTCATCAAAGCATTTATTTTATCTTAAAATATACTTT'
            'AACAGCTGATCAGGTTATCTTACTTATTCATGATTCCAATT',
            'acceptor':
            'TTTAACAGCTGATCAGGTTATCTTACTTATTCATGATTCCAATTTTTCAGACC',
            'exon':
            'ACCTCAGCAATCACCCAGCGGATAAGTCCTTGTTCCACTCTGACTAGCAGCAC'
            'TGCCTCTCCACCAGCCAGTAGCCCCTGCTCTACACTCCCACCCATCAGTACAA'
            'ATGCAACTGCCAAGGACTGCAGCTATGGGGCTGTTACTAGTCCAACCTCTACC'
            'CTTGAAAGCAGAGATAGTGGCATCATTG',
            'donor':
            'CATTGGTGAGTTGGTTTT',
            'donor_intron':
            'TGGTTTTTATATTGATAATTTTGTGTCCTTTTTTTCCTTTTTAAAATAATTAC'
            'ACAAGCTTAAGGTTTTTAAAAATTCTGCTTTTGAATTGTTG'
        },
        'alt_seq': {
            'acceptor_intron':
            'GCATGTAGTTTTTTTTTTCATCAAAGCATTTATTTTATCTTAAAATATACTTT'
            'AACAGCTGATCAGGTTATCTTACTTATTCATGATTCCAATT',
            'acceptor':
            'TTTAACAGCTGATCAGGTTATCTTACTTATTCATGATTCCAATTTTTCAGACC',
            'exon': 'ACCTCAGCAATCACCCAGCGGATAAGTCCTTGTTCCACTCTGACT'
            'AGCAGCACTGCCTCTCCACCAGCCAGTAGCCCCTGCTCTACACTCCCACCCAT'
            'CAGTACAAATGCAACTGCCAAGGACTGCAGCTATGGGGCTGTTACTAGTCCAA'
            'CCTCTACCCTTGAAAGCAGAGATAGTGGCATCATTG',
            'donor': 'CATTGATGAGTTGGTTTT',
            'donor_intron':
            'TGGTTTTTATATTGATAATTTTGTGTCCTTTTTTTCCTTTTTAAAATAATTAC'
            'ACAAGCTTAAGGTTTTTAAAAATTCTGCTTTTGAATTGTTG'
        }
    }

    d = next(dl)
    assert d['inputs']['seq'] == expected_snps_seq['seq']
    assert d['inputs']['mut_seq'] == expected_snps_seq['alt_seq']


def test_splicing_vcf_loads_snps(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False)

    expected_snps_seq = {
        'seq':
        'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTT'
        'TATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAG'
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGG'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
        'alt_seq':
        'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTT'
        'TATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAG'
        'TTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAG'
        'TACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGC'
        'TAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTA'
        'TGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC'
    }

    rows = list(dl)
    # d = rows[11]

    for i in rows:
        if '17:41276033:C>G' == i['metadata']['variant']['annotation']:
            d = i

    assert d['inputs']['seq'] == expected_snps_seq['seq']
    assert d['inputs']['mut_seq'] == expected_snps_seq['alt_seq']


def test_splicing_vcf_loads_deletions(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False)

    expected_snps_seq = {
        '17:41267740:TACTT>A': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCT'
            'AGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTTCA'
            'TAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAATTA'
        },
        '17:41267741:ACTT>AC': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCG'
            'TAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTT'
            'CATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41267742:CTTGCAAAATATGTGGTCACACTTTGTGGAGACAGGTTCCTTGATCAACTCCAGA>C': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGGT'
            'AAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTTC'
            'ATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41267742:CTT>C': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCG'
            'TAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTT'
            'CATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41267795:GAC>GA': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'ATAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAAT'
            'AAATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTATC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41267795:GA>G': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGCT'
            'GGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAA'
            'GTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCT'
            'TCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41267796:ACT>A': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATA'
            'AATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'CATAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAA'
            'TAAATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTTC'
            'TGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCA'
            'AGTAAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCC'
            'TTCATAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        }
    }

    rows = list(dl)

    for d in rows:
        variant = d['metadata']['variant']['annotation']

        if variant in expected_snps_seq:
            assert d['inputs']['seq'] == expected_snps_seq[variant]['seq']
            assert d['inputs']['mut_seq'] == expected_snps_seq[variant]['alt_seq']


def test_splicing_vcf_loads_insertions(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False)

    expected_snps_seq = {
        '17:41267795:G>GAA': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATAA'
            'ATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTCTG'
            'GAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAGT'
            'AAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTTCA'
            'TAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATAAAT'
            'TATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTTTCTG'
            'GAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAGT'
            'AAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTTCA'
            'TAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41267797:CT>CTAA': {
            'seq':
            'TAACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATAA'
            'ATTATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTAGTCTG'
            'GAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAGT'
            'AAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTTCA'
            'TAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT',
            'alt_seq':
            'ACAGCTCAAAGTTGAACTTATTCACTAAGAATAGCTTTATTTTTAAATAAAT'
            'TATTGAGCCTCATTTATTTTCTTTTTCTCCCCCCCTACCCTGCTTTAGTCTG'
            'GAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAGT'
            'AAGTTTGAATGTGTTATGTGGCTCCATTATTAGCTTTTGTTTTTGTCCTTCA'
            'TAACCCAGGAAACACCTAACTTTATAGAAGCTTTACTTTCTTCAAT'
        },
        '17:41276032:AC>ACCA': {
            'seq':
            'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTTTA'
            'TAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAGTTCA'
            'TTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAA'
            'ATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGGTAAGTCAG'
            'CACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTATGCAAATGAA'
            'CAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
            'alt_seq':
            'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTTTA'
            'TAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAGTTCA'
            'TTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAA'
            'ATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTGGTAAGTC'
            'AGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTATGCAAATG'
            'AACAGAATTGACCTTACATACTAGGGAAGAAAAGACATG'
        },
        '17:41276033:C>CCAGATG': {
            'seq':
            'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTTTA'
            'TAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAGTTCA'
            'TTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAA'
            'ATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGGTAAGTCAG'
            'CACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTATGCAAATGAA'
            'CAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
            'alt_seq':
            'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGTCATTTTA'
            'TAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTAAAGTTCA'
            'TTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAA'
            'ATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGCATCTGGTA'
            'AGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTATGCA'
            'AATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGA'
        },
        '17:41276132:A>ACT': {
            'seq': 'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAGTTGT'
            'CATTTTATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGTGTTA'
            'AAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTTGAAGAA'
            'GTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGGT'
            'AAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTCCTATGC'
            'AAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC',
            'alt_seq': 'CAAATCTTAAATTTACTTTATTTTAAAATGATAAAATGAAG'
            'TTGTCATTTTATAAACCTTTTAAAAAGATATATATATATGTTTTTCTAATGT'
            'GTTAAAGAGTTCATTGGAACAGAAAGAAATGGATTTATCTGCTCTTCGCGTT'
            'GAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCA'
            'TCTGGTAAGTCAGCACAAGAGTGTATTAATTTGGGATTCCTATGATTATCTC'
            'CTATGCAAATGAACAGAATTGACCTTACATACTAGGGAAGAAAAGACATGTC'
        }
    }

    rows = list(dl)

    for d in rows:
        variant = d['metadata']['variant']['annotation']

        if variant in expected_snps_seq:
            assert d['inputs']['seq'] == expected_snps_seq[variant]['seq']
            assert d['inputs']['mut_seq'] == expected_snps_seq[variant]['alt_seq']
