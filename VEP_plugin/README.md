# VEP MMSplice plugin user tutorial

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that runs MMSplice (modular modeling of splicing) which performs a set of prediction on splicing. The plugin requires MMSplice python package as an external dependency since it wraps mmsplice package as vep plugin. Thus, MMSplice package should be installed (check installation). Then, it automatically runs python server in background and analysis variant with python server.

For main MMSpilice documentation, check [main README.md](../README.md).

## Installation

Install ensemble vep, if you didn't already.

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```

For more [detail about vep installation](https://github.com/Ensembl/ensembl-vep)

Then, download mmsplice vep plugin. Alternatively, you can just install `MMSplice.pm` file under `VEP_plugin`.

```bash
git clone https://github.com/gagneurlab/MMSplice
cd VEP_plugin
```

Copy for mmsplice.pm file to your ensemble ensemble vep plugin directory.
```bash
cp MMSplice.pm ~/.vep/Plugins/
```

If your vep cache directory different then default, you need to copy MMSplice to your cachedir directory accordingly. For example:
```bash
cp MMSplice.pm /ensembl-vep/92/cachedir/
```

Lastly, install mmsplice python package to your machine or virtualenv.
```bash
pip install mmsplice
```

## Usage

if you are not already familiar with the usage of VEP plugins, please check [this documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) first.

Now, You can analyze your vcf with following comments using default mmsplice configuration.

```bash
./vep -i vcf_file.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --cache --port 3337
```

If your vep configurations are different then default, you need to add them as parameter.
For example, if your cache dir different then default please speficy it as follow:

```bash
./vep -i vcf_file.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --port 3337 --cache --dir /ensembl-vep/92/cachedir/
```

For further details about VEP plugin parameters, please check [this documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_plugin).

You may want to run mmsplice different then default parameter. The full list of parameters of mmsplice with their default is given below.
Please be cautious changing parameters. If the corresponding mmsplice model don't support the parameters it will throw errors.

```bash
 ./vep -i vcf_file.vcf --plugin MMSplice,[port_of_mmsplice_server=5000],[intronl_len=100],[intronr_len=80],[exon_cut_l=0],[exon_cut_r=0],[acceptor_intron_cut=6],[donor_intron_cut=3],[acceptor_intron_len=20],[acceptor_exon_len=3],[donor_exon_len=3],[donor_intron_len=6],[acceptor_intronM],[acceptorModelFile],[exonModelFile],[donorModelFile],[donor_intronModelFile]
```

## Results

If you are able to run the code above, it will produce an vcf file with following columns are added to INFO section:

```
mmsplice_ref_acceptor_intron => "acceptor intron score of reference sequence",
mmsplice_ref_acceptor => "acceptor score of reference sequence",
mmsplice_ref_exon => "exon score of reference sequence",
mmsplice_ref_donor => "donor score of reference sequence",
mmsplice_ref_donor_intron => "donor intron score of reference sequence",
mmsplice_alt_acceptor_intron => "acceptor intron score of variant sequence ",
mmsplice_alt_acceptor => "acceptor score of variant sequence",
mmsplice_alt_exon => "exon score of variant sequence",
mmsplice_alt_donor => "donor score of variant sequence",
mmsplice_alt_donor_intron => "donor intron score of variant sequence",
mmsplice_delta_logit_psi => "delta logit psi score of variant",
mmsplice_pathogenicity => "pathogenicity effect of variant"
```

The plugin don't filters any variant. Some of the variants may not have prediction because they are not matched. In this case, emtpy values are returned.

## Troubleshoot

### Gziped Vcf

VEP don't support gzip files. So if you get `gzip: stdout: Broken pipe` this error, or tring to run gziped vcf file, please unzip your file first. Then, use unzipped vcf version.

### Thread Safety

python server is not tread safe due to some technical limition in deep learning models. So, if you want to analysis two vcf file some time, please start two different python server with specifying different port number (check parameter of mmsplice plugin). For example:

```bash
./vep -i vcf_file1.vcf --plugin MMSplice,5000 ... & ./vep -i vcf_file1.vcf --plugin MMSplice,5001 ...
```
