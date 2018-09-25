# -*- coding: utf-8 -*-

"""Tests for `mmsplice` package."""
import pytest

from mmsplice import MMSplice

def test_mmsplice():
    x = {'seq': "ATGCGACGTACCCAGTAAAT",
         'intronl_len': 4,
         'intronr_len': 4}
    model = MMSplice()
    pred = model.predict(x)
    assert len(pred) == 5
