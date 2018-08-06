===
eis
===


.. image:: https://img.shields.io/pypi/v/eis.svg
        :target: https://pypi.python.org/pypi/eis

.. image:: https://img.shields.io/travis/s6juncheng/eis.svg
        :target: https://travis-ci.org/s6juncheng/eis

.. image:: https://readthedocs.org/projects/eis/badge/?version=latest
        :target: https://eis.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


.. image:: https://pyup.io/repos/github/s6juncheng/eis/shield.svg
     :target: https://pyup.io/repos/github/s6juncheng/eis/
     :alt: Updates



Predict splicing variant effect from VCF


* Free software: MIT license
* Documentation: https://eis.readthedocs.io.


## Usage example
------

Check notebooks/example.ipynb

```python
# Import
from eis.vcf_dataloader import SplicingVCFDataloader
from eis import Eis, predict_all_table
from eis.utils import max_varEff

# example files
gtf = 'tests/data/test.gtf'
vcf = 'tests/data/test.vcf.gz'
fasta = 'tests/data/hg19.nochr.chr17.fa'
gtfIntervalTree = '../tests/data/test.pkl' # pickle exon interval Tree

# dataloader to load variants from vcf
dl = SplicingVCFDataloader(gtf, 
                          fasta,
                          vcf,
                          out_file=gtfIntervalTree,
                          split_seq=False)

# Specify model
model = Eis(
    exon_cut_l=0,
    exon_cut_r=0,
    acceptor_intron_cut=6,
    donor_intron_cut=6,
    acceptor_intron_len=50,
    acceptor_exon_len=3,
    donor_exon_len=5,
    donor_intron_len=13)

 # Do prediction
 predictions = predict_all_table(model, dl, batch_size=1024, split_seq=False, assembly=False)

 # Summerize with maximum effect size
 predictionsMax = max_varEff(predictions)
```