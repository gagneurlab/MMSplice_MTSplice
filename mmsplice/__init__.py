# -*- coding: utf-8 -*-

"""Top-level package for mmsplice."""

__author__ = """Jun Cheng"""
__email__ = 'chengju@in.tum.de'
__version__ = '0.2.0'

from mmsplice.mmsplice import MMSplice, predict_all, writeVCF, predict_all_table
import mmsplice.generic
import mmsplice.IntervalTree
import mmsplice.vcf_dataloader
import mmsplice.exon_dataloader
from .utils import postproc

# export modules
from mmsplice.mmsplice import ACCEPTOR_INTRON, ACCEPTOR, EXON, DONOR, DONOR_INTRON, LINEAR_MODEL, LOGISTIC_MODEL
from keras.models import load_model

ACCEPTOR_INTRON_MODULE = load_model(ACCEPTOR_INTRON)
ACCEPTOR_MODULE = load_model(ACCEPTOR)
EXON_MODULE = load_model(EXON)
DONOR_MODULE = load_model(DONOR)
DONOR_INTRON_MODULE = load_model(DONOR_INTRON)