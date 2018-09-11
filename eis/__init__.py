# -*- coding: utf-8 -*-

"""Top-level package for eis."""

__author__ = """Jun Cheng"""
__email__ = 'chengju@in.tum.de'
__version__ = '0.2.0'

from eis.eis import Eis, predict_all, writeEIS, predict_all_table
import eis.generic
import eis.IntervalTree
import eis.vcf_dataloader
import eis.exon_dataloader
from .utils import postproc

# export modules
from eis.eis import ACCEPTOR_INTRON, ACCEPTOR, EXON, DONOR, DONOR_INTRON 
from keras.models import load_model

ACCEPTOR_INTRON_MODULE = load_model(ACCEPTOR_INTRON)
ACCEPTOR_MODULE = load_model(ACCEPTOR)
EXON_MODULE = load_model(EXON)
DONOR_MODULE = load_model(DONOR)
DONOR_INTRON_MODULE = load_model(DONOR_INTRON)