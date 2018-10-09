
# -*- coding: utf-8 -*-

"""Tests for `mmsplice` package."""
from mmsplice import MMSplice
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import predict_all_table

from conftest import gtf_file, fasta_file, variants


def test_mmsplice():
    x = {'seq': "ATGCGACGTACCCAGTAAAT",
         'intronl_len': 4,
         'intronr_len': 4}
    model = MMSplice()
    pred = model.predict(x)
    assert len(pred) == 5


def test_predict_all_table(vcf_path):
    model = MMSplice(
        exon_cut_l=0,
        exon_cut_r=0,
        acceptor_intron_cut=6,
        donor_intron_cut=6,
        acceptor_intron_len=50,
        acceptor_exon_len=3,
        donor_exon_len=5,
        donor_intron_len=13)

    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    df = predict_all_table(model, dl, batch_size=1024,
                           split_seq=False, assembly=True,
                           pathogenicity=True, splicing_efficiency=True)

    assert len(df['mmsplice_dlogitPsi']) == len(variants) - 1

    # dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    # df = predict_all_table(model, dl, batch_size=1024,
    #                        split_seq=True, assembly=True)
