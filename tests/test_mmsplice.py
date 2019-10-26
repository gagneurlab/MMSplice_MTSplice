"""Tests for `mmsplice` package."""
import pandas as pd
from mmsplice import MMSplice
from mmsplice.utils import encodeDNA
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice.exon_dataloader import ExonDataset
from mmsplice import predict_all_table
from conftest import gtf_file, fasta_file, variants, exon_file


def test_mmsplice():
    seq = 'ATGCGACGTACCCAGTAAAT'
    overhang = (4, 4)
    model = MMSplice()
    pred = model.predict_on_seq(seq, overhang)
    assert len(pred) == 5


def test_predict_save(vcf_path):
    pass


def test_predict_all_table(vcf_path):
    model = MMSplice()

    dl = SplicingVCFDataloader(gtf_file, fasta_file, vcf_path)
    df = predict_all_table(model, dl, pathogenicity=True,
                           splicing_efficiency=True)

    assert len(df['delta_logit_psi']) == len(variants) - 1


def test_predict_all_table_tissue_specific(vcf_path):
    model = MMSplice()
    dl = SplicingVCFDataloader(
        gtf_file, fasta_file, vcf_path, tissue_specific=True)
    df = predict_all_table(model, dl)
    assert len(df['delta_logit_psi']) == len(variants) - 1
    assert df.shape[1] == 7 + 10 + 56


def test_predict_all_table_exon_dataloader(vcf_path):
    model = MMSplice()
    df_exons = pd.read_csv(exon_file)
    dl = ExonDataset(exon_file, fasta_file)
    df = predict_all_table(model, dl, pathogenicity=True,
                           splicing_efficiency=True)

    assert len(df['delta_logit_psi']) == df_exons.shape[0]


def test_predict_all_table_tissue_specific_exon_dataloader():
    model = MMSplice()
    df_exons = pd.read_csv(exon_file)
    dl = ExonDataset(exon_file, fasta_file, tissue_specific=True)
    df = predict_all_table(model, dl)

    assert len(df['delta_logit_psi']) == df_exons.shape[0]
    assert df.shape[1] == 3 + 10 + 56


def test_exon_model_masking():
    model = MMSplice()

    preds = [
        model.exonM.predict(encodeDNA(['AAA']))[0][0],
        model.exonM.predict(encodeDNA(['AAA', 'CATACA']))[0][0],
        model.exonM.predict(encodeDNA(['AAA', 'CATACAGGAA']))[0][0]
    ]

    for i in preds:
        assert abs(preds[0] - i) < 1e-6
