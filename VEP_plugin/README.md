# VEP MMSplice plugin user tutorial

## Installation

Install ensemble vep, if you didn't already.

```
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```

For more [detail about vep installation](https://github.com/Ensembl/ensembl-vep)

Then, download mmsplice vep plugin. Alternatively, you can just install `MMSplice.pm` file under `VEP_plugin`.

```
git clone https://github.com/gagneurlab/MMSplice
cd VEP_plugin
```

Copy for mmsplice.pm file to your ensemble ensemble vep plugin directory.

```
cp MMSplice.pm ~/.vep/Plugins/
```

Lastly, install mmsplice python package to your machine or virtualenv.
```
pip install mmsplice
```

## Usage

Now, You can analyze your vcf with following comments using default mmsplice configuration.

```
./vep -i ~/dirOf/vcf_file.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --cache --port 3337
```

If your vep configurations are different then default, you need to add them as parameter.
For example, if your cache dir different then default please speficy it as follow:

```
vep -i test.vcf --plugin MMSplice --vcf --force --assembly GRCh37 --port 3337 --cache --dir /ensembl-vep/92/cachedir/
```

You may want to run mmsplice different then default parameter. The full list of parameters of mmsplice with their default is given below.
Please be cautious changing parameters. If the corresponding mmsplice model don't support the parameters it will throw errors.

```
 ./vep -i variants.vcf --plugin MMSplice,[port_of_mmsplice_server=5000],[intronl_len=100],[intronr_len=80],[exon_cut_l=0],[exon_cut_r=0],[acceptor_intron_cut=6],[donor_intron_cut=3],[acceptor_intron_len=20],[acceptor_exon_len=3],[donor_exon_len=3],[donor_intron_len=6],[acceptor_intronM],[acceptorModelFile],[exonModelFile],[donorModelFile],[donor_intronModelFile]
```
