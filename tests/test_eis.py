#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `eis` package."""

import pytest

from eis import Eis


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

     
def test_eis():
    x = {'seq': "ATGCGACGTACCCAGTAAAT",
     'intronl_len': 4,
    'intronr_len': 4}
    model = Eis()
    pred = model.predict(x)
    assert len(pred) == 5
    