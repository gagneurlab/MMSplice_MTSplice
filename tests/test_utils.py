import pyranges
from mmsplice.utils import pyrange_remove_chr_from_chrom_annotation


def test_pyrange_remove_chr_to_chrom_annotation():
    pr_exons = pyranges.data.exons()
    pyrange_remove_chr_from_chrom_annotation(pr_exons)
    assert set(pr_exons.Chromosome) == {'chrX', 'chrY'}
