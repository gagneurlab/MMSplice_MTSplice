import pyranges
from pybedtools import Interval
from mmsplice.utils import pyrange_remove_chr_from_chrom_annotation, Variant, \
    left_normalized, get_var_side


def test_pyrange_remove_chr_to_chrom_annotation():
    pr_exons = pyranges.data.exons()
    pyrange_remove_chr_from_chrom_annotation(pr_exons)
    assert set(pr_exons.Chromosome) == {'chrX', 'chrY'}


def test_Variant():
    v = Variant('1', 10, 'CA', ['CAGC'])
    assert v.CHROM == '1'
    assert v.POS == 10
    assert v.REF == 'CA'
    assert v.ALT[0] == 'CAGC'
    assert v.start == 9


def test_Variant_left_normalized():
    v = left_normalized(Variant('1', 10, 'CA', ['CAGC']))
    assert v.CHROM == '1'
    assert v.POS == 12
    assert v.REF == ''
    assert v.ALT[0] == 'GC'
    assert v.start == 11


def test_get_var_side():
    exon = Interval('chr1', 11, 20, strand='+')
    variant = Variant('chr1', 10, 'A', ['AGG'])
    assert get_var_side(variant, exon) == 'left'

    variant = Variant('chr1', 11, 'A', ['AGG'])
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
