# MMSplice & MTSplice
[![CircleCI](https://circleci.com/gh/gagneurlab/MMSplice_MTSplice.svg?style=svg)](https://circleci.com/gh/gagneurlab/MMSplice_MTSplice)
[![pypi](https://img.shields.io/pypi/v/mmsplice.svg)](https://pypi.python.org/pypi/mmsplice)

Predict (tissue-specific) splicing variant effect from VCF. MTSplice is integrated into MMSplice with the same API. 

Paper: Cheng et al. https://doi.org/10.1101/438986, https://www.biorxiv.org/content/10.1101/2020.06.07.138453v1

![MMSplice](https://raw.githubusercontent.com/kipoi/models/master/MMSplice/Model.png)
![MTSplice](https://raw.githubusercontent.com/s6juncheng/figshare/master/MTSplice.JPG)

## Installation
-----------------

External dependencies:
```bash
pip install cyvcf2 cython
```

Conda installation is recommended:
```bash
conda install cyvcf2 cython -y
```

```bash
pip install mmsplice
```

## Run MMSplice Online

You can run mmsplice with following google colab notebooks online:

- [run on vcf file](https://colab.research.google.com/drive/1Kw5rHMXaxXXsmE3WecxbXyGQJma80Eq6)

### Preparation
-----------------

#### 1. Prepare annotation (gtf) file
Standard human gene annotation file in GTF format can be downloaded from ensembl or gencode.
`MMSplice` can work directly with those files, however, some filtering is higly recommended.

- Filter for protein coding genes.

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
-------------------

Check [notebooks/example.ipynb](https://github.com/gagneurlab/MMSplice/blob/master/notebooks/example.ipynb)

To score variants (including indels), we suggest to use primarily the `deltaLogitPSI` predictions, which is the default output. The differential splicing efficiency (dse) model was trained from MMSplice modules and exonic variants from MaPSy, thus only the predictions for exonic variants are calibrated.

**MTSplice** To predict tissue-specific variant effect with MTSplice, specify `tissue_specific=True` in `SplicingVCFDataloader`. 

```python
# Import
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
from mmsplice.utils import max_varEff

# example files
gtf = 'tests/data/test.gtf'
vcf = 'tests/data/test.vcf.gz'
fasta = 'tests/data/hg19.nochr.chr17.fa'
csv = 'pred.csv'

# dataloader to load variants from vcf
dl = SplicingVCFDataloader(gtf, fasta, vcf, tissue_specific=False)

# Specify model
model = MMSplice()

# predict and save to csv file
predict_save(model, dl, csv, pathogenicity=True, splicing_efficiency=True)

# Or predict and return as df
predictions = predict_all_table(model, dl, pathogenicity=True, splicing_efficiency=True)

# Summerize with maximum effect size
predictionsMax = max_varEff(predictions)
```

### Output

Output of MMSplice is an tabular data which contains following described columns:

* `ID`: id string of the variant
* `delta_logit_psi`: The main score is predicted by MMSplice, which shows the effect of the variant on the inclusion level (PSI percent spliced in) of the exon. The score is on a logit scale.  If the score is positive, it shows that variant leads higher inclusion rate for the exon. If the score is negative, it shows that variant leads higher exclusion rate for the exon. If delta_logit_psi is bigger than 2 or smaller than -2, the effect of variant can be considered strong.
* `exons`: Genetics location of exon whose inclusion rate is effected by variant
* `exon_id`: Genetic id of exon whose inclusion rate is effected by variant
* `gene_id`: Genetic id of the gene which the exon belongs to. 
* `gene_name`:  Name of the gene which the exon belongs to. 
* `transcript_id`: Genetic id of the transcript which the exon belongs to. 
* `ref_acceptorIntron`: acceptor intron score of the reference sequence
* `ref_acceptor`: acceptor score of the reference sequence
* `ref_exon`: exon score of the reference sequence
* `ref_donor`: donor score of the reference sequence
* `ref_donorIntron`: donor intron score of the reference sequence
* `alt_acceptorIntron`: acceptor intron score of variant sequence
* `alt_acceptor`: acceptor score of the sequence with variant
* `alt_exon`: exon score of the sequence with variant
* `alt_donor`: donor score of the sequence with variant
* `alt_donorIntron`: donor intron score of the sequence with variant
* `pathogenicity`: Potential pathogenic effect of the variant.
* `efficiency`: The effect of the variant on the splicing efficiency of the exon.


## VEP Plugin

The VEP plugin wraps the prediction function from `mmsplice` python package. Please check documentation of vep plugin [under VEP_plugin/README.md](VEP_plugin/README.md).
