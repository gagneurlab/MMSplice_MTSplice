
# -*- coding: utf-8 -*-

"""Tests for `mmsplice` package."""
from concise.preprocessing import encodeDNA
from mmsplice import MMSplice
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import predict_all_table

from conftest import gtf_file, fasta_file, variants


def test_mmsplice():
    seq = 'ATGCGACGTACCCAGTAAAT'
    overhang = (4, 4)
    model = MMSplice()
    pred = model.predict(seq, overhang)
    assert len(pred) == 5


def test_predict_save(vcf_path):
    pass


def test_predict_all_table(vcf_path):
    model = MMSplice()

    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    df = predict_all_table(model, dl, pathogenicity=True,
                           splicing_efficiency=True)

    assert len(df['delta_logit_psi']) == len(variants) - 1


def test_exon_model_masking():
    model = MMSplice()

    preds = [
        model.exonM.predict(encodeDNA(['AAA']))[0][0],
        model.exonM.predict(encodeDNA(['AAA', 'CATACA']))[0][0],
        model.exonM.predict(encodeDNA(['AAA', 'CATACAGGAA']))[0][0]
    ]

    for i in preds:
        assert abs(preds[0] - i) < 1e-6
