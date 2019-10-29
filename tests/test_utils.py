import numpy as np
import pyranges
from kipoiseq.dataclasses import Interval, Variant
from mmsplice.utils import pyrange_remove_chr_from_chrom_annotation, \
    left_normalized, get_var_side, encodeDNA


def test_pyrange_remove_chr_to_chrom_annotation():
    pr_exons = pyranges.data.exons()
    pyrange_remove_chr_from_chrom_annotation(pr_exons)
    assert set(pr_exons.Chromosome) == {'chrX', 'chrY'}


def test_Variant_left_normalized():
    v = left_normalized(Variant('1', 10, 'CA', 'CAGC'))

    assert v.chrom == '1'
    assert v.pos == 12
    assert v.ref == ''
    assert v.alt == 'GC'
    assert v.start == 11


def test_get_var_side():
    exon = Interval('chr1', 11, 20, strand='+')
    variant = Variant('chr1', 10, 'A', 'AGG')
    assert get_var_side(variant, exon) == 'left'

    variant = Variant('chr1', 11, 'A', 'AGG')
    assert get_var_side(variant, exon) == "exon"

    variant = Variant('chr1', 20, 'A', 'AGG')
    assert get_var_side(variant, exon) == 'right'

    variant = Variant('chr1', 13, 'A', 'AGG')
    assert get_var_side(variant, exon) == "exon"

    exon = Interval('chr1', 11, 20, strand='-')
    variant = Variant('chr1', 8, 'A', 'AGG')
    assert get_var_side(variant, exon) == "right"

    variant = Variant('chr1', 11, 'A', 'AGG')
    assert get_var_side(variant, exon) == "exon"

    variant = Variant('chr1', 20, 'A', 'AGG')
    assert get_var_side(variant, exon) == 'left'


def test_encodeDNA():
    seq_vec = ['ACGTN']
    np.testing.assert_array_equal(
        encodeDNA(seq_vec),
        np.array([[[1., 0., 0., 0.],
                   [0., 1., 0., 0.],
                   [0., 0., 1., 0.],
                   [0., 0., 0., 1.],
                   [0., 0., 0., 0.]]])
    )

    seq_vec = ['AA', 'ATT', 'ATTCGG']
    np.testing.assert_array_equal(
        encodeDNA(seq_vec),
        np.array([[[1., 0., 0., 0.],
                   [1., 0., 0., 0.],
                   [0., 0., 0., 0.],
                   [0., 0., 0., 0.],
                   [0., 0., 0., 0.],
                   [0., 0., 0., 0.]],

                  [[1., 0., 0., 0.],
                   [0., 0., 0., 1.],
                   [0., 0., 0., 1.],
                   [0., 0., 0., 0.],
                   [0., 0., 0., 0.],
                   [0., 0., 0., 0.]],

                  [[1., 0., 0., 0.],
                   [0., 0., 0., 1.],
                   [0., 0., 0., 1.],
                   [0., 1., 0., 0.],
                   [0., 0., 1., 0.],
                   [0., 0., 1., 0.]]])
    )
