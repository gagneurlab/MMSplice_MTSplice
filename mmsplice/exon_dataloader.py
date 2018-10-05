from kipoi.data import Dataset
from kipoi.metadata import GenomicRanges

from concise.preprocessing import encodeDNA
import pandas as pd
import warnings
from pyfaidx import Fasta
import six

from .generic import Variant, get_var_side
from .vcf_dataloader import ExonInterval


def readExon(exonfile, **kwargs):
    exons = pd.read_csv(exonfile, **kwargs)
    exons = exons.rename(columns={
        "hg19_variant_position": "POS",
        "variant_position": "POS",
        "reference": "REF",
        "variant": "ALT",
        "start": "Exon_Start",
        "end": "Exon_End",
        "exon_start": "Exon_Start",
        "exon_end": "Exon_End",
        "chr": "CHROM",
        "seqnames": "CHROM",
        "chromosome": "CHROM"
    })
    return exons


class ExonDataset(Dataset):
    ''' This class work with files that is a list of alternative exons
    exonfile: need to have [chr, Exon_Start, Exon_End, strand, ID, REF, ALT, POS], best also have transcript_id, gene_id
    split_to: order of spliting and return feature sequence
    sep: exonfile delimiter
    maxExonLen: exon length to consider, shorter sequence will be padded, langer one will be trimmed at the start end
    '''

    def __init__(self,
                 exonfile,
                 fasta_file,
                 overhang=(20, 80),  # overhang from the exon
                 exon_cut_l=0,
                 exon_cut_r=0,
                 acceptor_intron_cut=6,
                 donor_intron_cut=6,
                 acceptor_intron_len=50,
                 acceptor_exon_len=3,
                 donor_exon_len=5,
                 donor_intron_len=13,
                 sep=',',
                 split_seq=False,
                 encode=False
                 ):
        if isinstance(fasta_file, six.string_types):
            fasta = Fasta(fasta_file, as_raw=False)
        self.fasta = fasta
        self.overhang = overhang
        self.sep = sep
        self.split_seq = split_seq
        exons = readExon(exonfile, sep=self.sep)

        self.encode = encode
        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.acceptor_intron_cut = acceptor_intron_cut
        self.donor_intron_cut = donor_intron_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_exon_len = donor_exon_len
        self.donor_intron_len = donor_intron_len

        ########
        # Specific part for Vex-seq
        # variant side
        exons['side'] = list(map(get_var_side,
                                 zip(exons['POS'],
                                     exons['REF'],
                                     exons['ALT'],
                                     exons['Exon_Start'],
                                     exons['Exon_End'],
                                     exons['strand'])))
        ########
        self.exons = exons

    def __len__(self):
        return len(self.exons)

    def __getitem__(self, idx):
        exon = self.exons.iloc[idx]
        ########
        # Specific part for Vex-seq
        # variant side
        variant = Variant(CHROM=exon.CHROM,
                          POS=exon.POS,
                          REF=exon.REF,
                          ALT=exon.ALT,
                          ID=exon.ID,
                          strand=exon.strand,
                          side=exon.side)
        ########

        attributes = {}
        try:
            attributes['transcript_id'] = [exon.transcript_id]
        except:
            attributes['transcript_id'] = [""]
        try:
            attributes['gene_id'] = [exon.gene_id]
        except:
            attributes['gene_id'] = [""]
        try:
            attributes['exon_id'] = [exon.CHROM + ':' +
                                     exon.Exon_Start + '-' + exon.Exon_End]
        except:
            attributes['exon_id'] = [""]
        try:
            attributes['biotype'] = exon.biotype
        except:
            attributes['biotype'] = ""
        try:
            attributes['order'] = exon.order
        except:
            attributes['order'] = ""

        exon = ExonInterval.from_exonfile(
            exon, attributes, overhang=self.overhang)

        out = {}
        out['inputs'] = {}
        out['mut_inputs'] = {}
        if self.split_seq:
            out['inputs']['seq'] = self.split(exon.get_seq(self.fasta))
            out['mut_inputs']['seq'] = self.split(
                exon.get_mut_seq(self.fasta, variant).upper())
        else:
            out['inputs']['seq'] = exon.get_seq(self.fasta)
            out['mut_inputs']['seq'] = exon.get_mut_seq(
                self.fasta, variant).upper()
        out['inputs']['intronl_len'] = self.overhang[0]
        out['inputs']['intronr_len'] = self.overhang[1]
        out['mut_inputs']['intronl_len'] = self.overhang[0]
        out['mut_inputs']['intronr_len'] = self.overhang[1]

        out['metadata'] = {}
        out['metadata']['gene_id'] = exon.gene_id
        out['metadata']['transcript_id'] = exon.transcript_id
        out['metadata']['biotype'] = attributes['biotype']
        out['metadata']['order'] = exon.order
        out['metadata']['ranges'] = GenomicRanges(exon.seqid,  # exon is now object of class ExonInterval
                                                  exon.start - 1,  # for kipoi 0-base
                                                  exon.end,  # actual got sequence coordinates
                                                  exon.gene_id,
                                                  exon.strand)
        return out

    def split(self, x, overhang):
        ''' x: a sequence to split
        '''
        intronl_len, intronr_len = overhang
        # need to pad N if left seq not enough long
        lackl = self.acceptor_intron_len - intronl_len
        if lackl >= 0:
            x = "N" * (lackl + 1) + x
            intronl_len += lackl + 1
        lackr = self.donor_intron_len - intronr_len
        if lackr >= 0:
            x = x + "N" * (lackr + 1)
            intronr_len += lackr + 1
        acceptor_intron = x[:intronl_len - self.acceptor_intron_cut]
        acceptor = x[(intronl_len - self.acceptor_intron_len): (intronl_len + self.acceptor_exon_len)]
        exon = x[(intronl_len + self.exon_cut_l): (-intronr_len - self.exon_cut_r)]
        donor = x[(-intronr_len - self.donor_exon_len): (-intronr_len + self.donor_intron_len)]
        donor_intron = x[-intronr_len + self.donor_intron_cut:]
        if donor[self.donor_exon_len:self.donor_exon_len + 2] != "GT":
            warnings.warn("None GT donor", UserWarning)
        if acceptor[self.acceptor_intron_len - 2:self.acceptor_intron_len] != "AG":
            warnings.warn("None AG donor", UserWarning)

        if self.encode:
            return {
                "acceptor_intron": encodeDNA([acceptor_intron]),
                "acceptor": encodeDNA([acceptor]),
                "exon": encodeDNA([exon]),
                "donor": encodeDNA([donor]),
                "donor_intron": encodeDNA([donor_intron])
            }
        else:
            return {
                "acceptor_intron": acceptor_intron,
                "acceptor": acceptor,
                "exon": exon,
                "donor": donor,
                "donor_intron": donor_intron
            }
