# -*- coding: utf-8 -*-

"""Tests for `mmsplice` package."""
import json
import pytest

from mmsplice import MMSplice, predict_all_table
from mmsplice.vcf_dataloader import SplicingVCFDataloader, GenerateExonIntervalTree
from mmsplice.api import app


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
    pass


def test_mmsplice():
    x = {'seq': "ATGCGACGTACCCAGTAAAT",
         'intronl_len': 4,
         'intronr_len': 4}
    model = MMSplice()
    pred = model.predict(x)
    assert len(pred) == 5

# test vcf dataloader
# clinvar_20180429.filtered.BRCA1.vcf.gz
# clinvar_20180429.filtered.BRCA1.vcf.gz.tbi
# Homo_sapiens.GRCh37.75.uniq_exon.BRCA1.gtf

# def test_GenerateExonIntervalTree():
#     gtf = 'data/test.gtf'
#     # fasta = 'data/hg19.nochr.chr17.fa'
#     exonTree = GenerateExonIntervalTree(gtf,
#         out_file = 'data/gtfIntervalTree.pkl')
#     return exonTree
# 17:41267742:CTT:['C']


def test_api():
    with app.test_client() as c:

        resp = c.post('/create-model', json={})
        assert resp.status_code == 200

        resp = c.post('/psi-score', json={
            'seq': "ATGCGACGTACCCAGTAAAT",
            'intronl_len': 4,
            'intronr_len': 4
        })
        assert len(resp.data.decode("utf-8").split(',')) == 5
