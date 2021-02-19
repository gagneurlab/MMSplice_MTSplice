"""Top-level package for mmsplice."""

__author__ = """Jun Cheng & M.Hasan Celik"""
__email__ = 's6juncheng@gmail.com'
__version__ = '2.1.1'

from keras.models import load_model
from mmsplice.mmsplice import MMSplice, \
    writeVCF, \
    predict_save, \
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

from mmsplice.mtsplice import MTSPLICE, MTSplice

__all__ = [
    'load_model',
    'MMSplice',
    'writeVCF',
    'predict_save',
    'predict_all_table',
    'ACCEPTOR_INTRON',
    'ACCEPTOR',
    'EXON',
    'EXON3',
    'DONOR',
    'DONOR_INTRON',
    'LINEAR_MODEL',
    'LOGISTIC_MODEL',
    'EFFICIENCY_MODEL',
    'MTSPLICE',
    'MTSplice'
]
