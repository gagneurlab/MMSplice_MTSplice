# mmsplice

[![pypi](https://img.shields.io/pypi/v/mmsplice.svg)](https://pypi.python.org/pypi/mmsplice)

Predict splicing variant effect from VCF

Paper: Cheng et al. https://doi.org/10.1101/438986

![MMSplice](https://raw.githubusercontent.com/kipoi/models/master/MMSplice/Model.png)


## Usage example
------

```bash
pip install mmsplice
```

### Preparation
------

#### 1. Prepare annotation (gtf) file
Standard human gene annotation file in GTF format can be downloaded from ensembl or gencode.
`MMSplice` can work directly with those files, however, some filtering is higly recommended.

- Filter for protein coding genes.
- Filter out duplicated exons. The same exon can be annotated multiple times if it appears in multiple transcripts. 
  This will cause duplicated predictions.

We provide a filtered version [here](https://raw.githubusercontent.com/gagneurlab/MMSplice_paper/master/data/shared/Homo_sapiens.GRCh37.75.chr.uniq_exon.gtf.gz). 
Note this version has chromosome names in the format `chr*`. You may need to remove them to match the chromosome names in your fasta file.

#### 2. Prepare variant (VCF) file
A correctly formatted VCF file with work with `MMSplice`, however the following steps will make it less prone to false positives:

- Quality filtering. Low quality variants leads to unreliable predictions.
- Avoid presenting multiple variants in one line by splitting them into multiple lines. Example code to do it:
  ```bash
  bcftools norm -m-both -o out.vcf in.vcf.gz
  ```
- Left-normalization. For instance, GGCA-->GG is not left-normalized while GCA-->G is. Details for unified representation of genetic variants see [Tan et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481842/)
  ```bash
  bcftools norm -f reference.fasta -o out.vcf in.vcf
  ```
  
#### 3. Prepare reference genome (fasta) file
Human reference fasta file can be downloaded from ensembl/gencode. Make sure the chromosome name matches with GTF annotation file you use.


### Example code
------

Check [notebooks/example.ipynb](https://github.com/gagneurlab/MMSplice/blob/master/notebooks/example.ipynb)

To score variants (including indels), we suggest to use primarily the `deltaLogitPSI` predictions, which is the default output. The differential splicing efficiency (dse) model was trained from MMSplice modules and exonic variants from MaPSy, thus only the predictions for exonic variants are calibrated.

```python
# Import
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
from mmsplice.utils import max_varEff

# example files
gtf = 'tests/data/test.gtf'
vcf = 'tests/data/test.vcf.gz'
fasta = 'tests/data/hg19.nochr.chr17.fa'
gtfIntervalTree = 'tests/data/test.pkl' # pickle exon interval Tree

# dataloader to load variants from vcf
dl = SplicingVCFDataloader(gtf,
                          fasta,
                          vcf,
                          out_file=gtfIntervalTree, # same pikled gtf IntervalTree
                          split_seq=False)

# Specify model
model = MMSplice(
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

## VEP Plugin

The VEP plugin wraps the prediction function from `mmsplice` python package. Please check documentation of vep plugin [under VEP_plugin/README.md](VEP_plugin/README.md).
