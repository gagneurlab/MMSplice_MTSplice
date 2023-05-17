"""Top-level package for mmsplice."""

__author__ = """Jun Cheng & M.Hasan Celik"""
__email__ = 'chengju@in.tum.de'
__version__ = '2.4.0'

from tensorflow.keras.models import load_model
from mmsplice.mmsplice import MMSplice, \
    writeVCF, \
    predict_save, \
    predict_all_table, \
    ACCEPTOR_INTRON, \
    ACCEPTOR, \
    EXON, \
    EXON3,\
    DONOR, \
    DONOR_INTRON


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
    'MTSPLICE',
    'MTSplice'
]
