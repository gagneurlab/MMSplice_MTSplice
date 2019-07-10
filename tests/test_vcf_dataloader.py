import numpy as np
from kipoiseq.extractors import MultiSampleVCF
from mmsplice.vcf_dataloader import SplicingVCFDataloader, \
    read_exon_pyranges, batch_iter_vcf, variants_to_pyranges, \
    read_vcf_pyranges

from conftest import gtf_file, fasta_file, snps, deletions, \
    insertions, variants, vcf_file


def test_read_exon_pyranges():
    df_exons = read_exon_pyranges(gtf_file).df

    assert df_exons.shape[0] == 68
    exon = df_exons[df_exons['exon_id'] == 'ENSE00001831829'].iloc[0]
    assert exon['Start'] == 41196822
    assert exon['End'] == 41197819 + 100


def test_batch_iter_vcf(vcf_path):
    batchs = list(batch_iter_vcf(vcf_path, 10))
    assert sum(len(i) for i in batchs) == len(variants)


def test_variants_to_pyranges(vcf_path):
    variants = list(MultiSampleVCF(vcf_path))
    df = variants_to_pyranges(variants).df
    assert df.shape[0] == len(variants)

    v = df.iloc[0]
    assert v.Chromosome == '13'
    assert v.Start == 32953886
    assert v.End == 32953889
    assert v.variant.REF == 'GTT'
    assert v.variant.ALT[0] == 'AA'


def test_read_vcf_pyranges(vcf_path):
    batchs = list(read_vcf_pyranges(vcf_path, batch_size=10))
    assert sum(i.df.shape[0] for i in batchs) == len(variants)


def test__chech_chrom_annotation():
    dl = SplicingVCFDataloader('grch37', fasta_file, vcf_file)
    chroms = {str(i) for i in range(1, 22)}.union(['X', 'Y', 'M'])
    assert len(chroms.difference(set(dl.pr_exons.Chromosome))) == 0
    assert sum(1 for i in dl) > 0


def test__encode_seq(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    encoded = dl._encode_seq({'acceptor': 'ATT'})
    np.testing.assert_array_equal(
        encoded['acceptor'][0],
        np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1.],
                  [0., 0., 0., 1.]])
    )
    encoded = dl._encode_seq({'acceptor': 'aTT'})
    np.testing.assert_array_equal(
        encoded['acceptor'][0],
        np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1.],
                  [0., 0., 0., 1.]])
    )


def test__encode_batch_seq(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    encoded = dl._encode_batch_seq({'acceptor': np.array(['ATT'])})
    np.testing.assert_array_equal(
        encoded['acceptor'][0],
        np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1.],
                  [0., 0., 0., 1.]])
    )
    encoded = dl._encode_batch_seq({'acceptor': np.array(['aTT'])})
    np.testing.assert_array_equal(
        encoded['acceptor'][0],
        np.array([[1., 0., 0., 0.],
                  [0., 0., 0., 1.],
                  [0., 0., 0., 1.]])
    )


def test_benchmark_SplicingVCFDataloader(benchmark, vcf_path):
    benchmark(SplicingVCFDataloader, gtf_file, fasta_file, vcf_path)


def test_splicing_vcf_dataloader_prebuild_grch37(vcf_path):
    dl = SplicingVCFDataloader('grch37', fasta_file, vcf_path)


def test_splicing_vcf_dataloader_prebuild_grch38(vcf_path):
    dl = SplicingVCFDataloader('grch38', fasta_file, vcf_path)


def test_splicing_vcf_loads_all(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    assert sum(1 for i in dl) == len(variants) - 1


def test_splicing_vcf_loads_snps(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False)

    expected_snps_seq = [
        {
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
    ]

    for i in range(len(snps)):
        d = next(dl)
        print(d)
        print(d['metadata']['exon']['start'])
        print(d['metadata']['exon']['end'])
        assert d['inputs']['seq'] == expected_snps_seq[i]['seq']
        assert d['inputs']['mut_seq'] == expected_snps_seq[i]['alt_seq']


def test_splicing_vcf_loads_deletions(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False)

    expected_snps_seq = [
        {
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
        {
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
        {
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
        {
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
        {
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
        },
        {
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
        {
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
        }
    ]

    for i in range(len(snps)):
        d = next(dl)

    for i in range(len(deletions) - 1):
        d = next(dl)

        print('Variant position:', d['metadata']['variant']['POS'])
        print('Interval:',
              d['metadata']['exon']['start'],
              '-',
              d['metadata']['exon']['end'])
        print(d)
        assert d['inputs']['seq'] == expected_snps_seq[i]['seq']
        assert d['inputs']['mut_seq'] == expected_snps_seq[i]['alt_seq']


def test_splicing_vcf_loads_insertions(vcf_path):
    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path,
                               split_seq=False, encode=False)

    for i in range(len(snps) + len(deletions) - 1):
        d = next(dl)

    expected_snps_seq = [
        {
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
        {
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
        {
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
        {
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
        {
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
    ]

    for i in range(len(insertions)):

        d = next(dl)

        print(d)
        print(d['metadata']['exon']['start'])
        print(d['metadata']['exon']['end'])
        assert d['inputs']['seq'] == expected_snps_seq[i]['seq']
        assert d['inputs']['mut_seq'] == expected_snps_seq[i]['alt_seq']
