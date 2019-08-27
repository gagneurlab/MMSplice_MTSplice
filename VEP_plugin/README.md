# VEP MMSplice plugin user tutorial

 This is a plugin for the Ensembl Variant Effect Predictor (VEP) that runs MMSplice (modular modeling of splicing) which performs a set of prediction on splicing. The plugin requires MMSplice python package as an external dependency since it wraps mmsplice package as vep plugin. Thus, MMSplice package should be installed (check [MMSplice installation](../README.md)).

For main MMSplice documentation, check [main README.md](../README.md).

## Installation

Install ensemble vep, if you didn't already.

```bash
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```

For more [detail about vep installation](https://github.com/Ensembl/ensembl-vep)

Then, download mmsplice vep plugin. Alternatively, you can just copy `MMSplice.pm` file under `VEP_plugin`.

```bash
git clone https://github.com/gagneurlab/MMSplice
cd VEP_plugin
# Copy for mmsplice.pm file to your ensemble ensemble vep plugin directory.
cp MMSplice.pm ~/.vep/Plugins/
```

If your vep cache directory differs from the default, you need to copy `MMSplice.pm` to your cachedir directory accordingly. For example:
```bash
cp MMSplice.pm /ensembl-vep/92/cachedir/Plugins/
```

Lastly, install mmsplice python package.
```bash
pip install mmsplice
```

## Usage

if you are not already familiar with the usage of VEP plugins, please check [this documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) first.

Now, You can analyze your vcf with following command using the default mmsplice configuration.

```bash
vep -i vcf_file.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --cache --port 3337
```

If your vep configurations are different than default, you need to add them as parameter.
For example, if your cache dir different than default please speficy it as follow:

```bash
vep -i vcf_file.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --port 3337 --cache --dir /ensembl-vep/92/cachedir/
```

For further details about VEP plugin parameters, please check [this documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_plugin).

You may want to run mmsplice with parameters differ from the default. The full list of parameters of mmsplice with their default is given below.
Please be cautious changing parameters. If the corresponding mmsplice model don't support the parameters it will throw errors.

```bash
vep -i vcf_file.vcf --plugin MMSplice,[port_of_mmsplice_server=5000],[intronl_len=100],[intronr_len=80],[exon_cut_l=0],[exon_cut_r=0],[acceptor_intron_cut=6],[donor_intron_cut=3],[acceptor_intron_len=20],[acceptor_exon_len=3],[donor_exon_len=3],[donor_intron_len=6],[acceptor_intronM],[acceptorModelFile],[exonModelFile],[donorModelFile],[donor_intronModelFile]
```

## Docker Usage

Also, dockerfile which contains VEP and mmsplice is provided with the repository. Docker can be built:
```
docker build -t mmsplice .
```

Now, you can attach the docker container and analyze your data with following usage section:

```bash
docker run -t -i mmsplice /bin/bash
```

If you already have .vep cache files, you may mount it to docker:
```bash
docker run -t -i -v $HOME/your_vep_data:/root/.vep ensemblorg/ensembl-vep /bin/bash

```

If you do not want to attach to docker container then you may mount all the necessary directories and analyze your vcf with:
```bash
sudo docker run -t -i -v ~/.vep:/root/.vep -v ~/Projects/MMSplice/tests/data:/data -v ~/Desktop/outputs:/opt/vep/src/ensembl-vep/outputs mmsplice vep -i /data/test.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --cache --port 3337 -o outputs/results.txt
```

If you want to use docker without dealing with cache files, you can use vep database and stdin-out:
```
cat tests/data/test.vcf | sudo docker run -i mmsplice vep --plugin MMSplice --format vcf --assembly GRCh37 --database --port 3337 --vcf -o STDOUT | tee variant_effect_output.txt
sed -n '/##fileformat=VCFv4.0/,$p' variant_effect_output.txt > variant_effect_output.vcf
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

### Forking and Thread-safety

MMSplice VEP plugin is not thread-safe so avoid running with `fork` parameter to prevent unintended behaviors.
