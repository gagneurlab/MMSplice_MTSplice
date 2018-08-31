# -*- coding: utf-8 -*-

"""Top-level package for eis."""

__author__ = """Jun Cheng"""
__email__ = 'chengju@in.tum.de'
__version__ = '0.2.0'

from eis.eis import Eis, predict_all, writeEIS, predict_all_table
import eis.generic
import eis.IntervalTree
import eis.vcf_dataloader
from .utils import postproc