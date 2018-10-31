# -*- coding: utf-8 -*-

"""Top-level package for mmsplice."""

__author__ = """Jun Cheng"""
__email__ = 'chengju@in.tum.de'
__version__ = '0.2.6'

from keras.models import load_model
import mmsplice.generic
import mmsplice.interval_tree
import mmsplice.vcf_dataloader
import mmsplice.exon_dataloader
from .utils import postproc
from mmsplice.mmsplice import MMSplice, \
    predict_all, \
    writeVCF, \
    predict_all_table, \
    ACCEPTOR_INTRON, \
    ACCEPTOR, \
    EXON, \
    EXON3,\
    DONOR, \
    DONOR_INTRON, \
    LINEAR_MODEL, \
    LOGISTIC_MODEL,\
    EFFICIENCY_MODEL

__all__ = [
    'load_model',
    'MMSplice',
    'predict_all',
    'writeVCF',
    'predict_all_table',
    'ACCEPTOR_INTRON',
    'ACCEPTOR',
    'EXON',
    'EXON3',
    'DONOR',
    'DONOR_INTRON',
    'LINEAR_MODEL',
    'LOGISTIC_MODEL',
    'EFFICIENCY_MODEL'
]
