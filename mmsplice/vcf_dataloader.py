import pickle
import warnings
import functools

import gffutils
from cyvcf2 import VCF
from pyfaidx import Fasta
from kipoi.data import SampleIterator
from kipoi.metadata import GenomicRanges
from concise.preprocessing import encodeDNA

from .generic import Variant, get_var_side
from .interval_tree import IntervalTree, Interval


class VariantInterval(Interval):

    def __init__(self, REF=None, ALT=None, **kwargs):
        super().__init__(**kwargs)
        self.REF = REF
        self.ALT = ALT

    @classmethod
    def from_Variant(cls, variant):
        # Only support one alternative.
        # If multiple alternative, need to split into multiple variants
        if len(variant.ALT) != 1:
            warnings.warn('%s has more than one alternative sequence,'
                         'only the first one was taken,'
                         'split into mutliple variants with bedtools' % variant.__repr__(), UserWarning)
        if variant.ID is None:
            ID = '.'
        else:
            ID = variant.ID
        return cls(chrom=variant.CHROM,
                   start=variant.POS,
                   end=variant.POS,
                   strand="*",
                   name=ID,
                   REF=variant.REF,
                   ALT=variant.ALT[0])

    def to_Variant(self, strand, side):
        # convert to my Variant class, which needs strand and side
        assert self.start == self.end
        var = Variant(CHROM=self.chrom,
                      POS=self.start,
                      REF=self.REF,
                      ALT=self.ALT,
                      strand=strand,
                      ID=self.name)
        var.side = side
        return var


class ExonInterval(gffutils.Feature):

    def __init__(self, order=-1, **kwargs):
        super().__init__(**kwargs)
        self.order = order
        self.name = self.attributes["exon_id"][0]
        self.transcript_id = self.attributes["transcript_id"][0]
        self.gene_id = self.attributes["gene_id"][0]
        # self.biotype = self.attributes["gene_biotype"][0]
        self.Exon_Start = self.start  # with overhang, will be set later
        self.Exon_End = self.end  # with overhang
        self.overhang = (0, 0)
        self._isLast = False

    @property
    def isLast(self):
        return self._isLast

    @isLast.setter
    def isLast(self, value):
        self._isLast = value

    @property
    def isFirst(self):
        return self.order == 1

    @property
    def grange(self):
        return GenomicRanges(self.chrom,
                             self.start,
                             self.end,
                             self.transcript_id,
                             self.strand)

    def __str__(self):
        return '{0}_{1}_{2}:{3}'.format(
            self.chrom, self.Exon_Start, self.Exon_End, self.strand)

    @property
    def to_dict(self):
        return {'isLast': self.isLast,
                'isFirst': self.isFirst,
                'order': self.order,
                'name': self.name,
                'gene_id': self.gene_id,
                'Exon_Start': self.Exon_Start,
                'Exon_End': self.Exon_End,
                'intronl_len': self.overhang[0],
                'intronr_len': self.overhang[1],
                'seqid': self.seqid,
                'strand': self.strand,
                'start': self.start,
                'end': self.end}

    @classmethod
    def from_Feature(cls,
                     feature,
                     overhang=(0, 0)):
        # convert gffutils.Feature to ExonInterval
        iv = cls(seqid=feature.chrom,
                 source=feature.source,
                 start=feature.start,  # exon start
                 end=feature.end,  # exon end
                 strand=feature.strand,
                 frame=feature.frame,
                 attributes=feature.attributes,
                 order=int(feature.attributes['exon_number'][0]))
        if iv.strand == "+":
            iv.start = iv.Exon_Start - overhang[0]
            iv.end = iv.Exon_End + overhang[1]
        else:
            iv.start = iv.Exon_Start - overhang[1]
            iv.end = iv.Exon_End + overhang[0]
        iv.overhang = overhang
        return iv

    @classmethod
    def from_exonfile(cls, exon, attributes, overhang=(0, 0)):
        iv = cls(seqid=exon.CHROM,
                 source='',
                 start=exon.Exon_Start,  # exon start
                 end=exon.Exon_End,  # exon end
                 strand=exon.strand,
                 frame='',
                 attributes=attributes,
                 order=attributes['order'])
        if iv.strand == "+":
            iv.start = iv.Exon_Start - overhang[0]
            iv.end = iv.Exon_End + overhang[1]
        else:
            iv.start = iv.Exon_Start - overhang[1]
            iv.end = iv.Exon_End + overhang[0]
        iv.overhang = overhang
        return iv

    def get_seq(self, fasta, use_strand=True):
        seq = self.sequence(fasta, use_strand=use_strand)
        seq = seq.upper()
        return seq

    def _ref_check(self, fasta, variant):
        '''Checks that ref seq matchs with fasta'''
        ref_check = fasta.get_seq(
            self.chrom,
            variant.POS, variant.POS +
            len(variant.REF) - 1,  # -1 because 1 based
            self.strand == '-').seq  # last option corresponds to rc= in Fasta.get_seq function from pyfaidx
        return ref_check.upper() != variant.REF

    def _get_snp_seq(self, seq, variant, position):
        '''Seq of SNP, or more than one nt equal length substitution'''
        mut_seq = seq[:max(0, position)] + variant.ALT + \
            seq[position + len(variant.ALT):]
        if len(mut_seq) > len(seq):
            mut_seq = mut_seq[:len(seq)]
            assert len(mut_seq) == len(seq), variant
        return mut_seq

    def _get_insertion_seq(self, seq, variant, position):
        # insertion
        if len(variant.REF) + position > len(seq):
            # The actual variant position exceeded the retrieved sequence
            return seq
        mut_seq = seq[:max(0, position)] + variant.ALT + \
            seq[position + len(variant.REF):]
        assert len(seq) - len(mut_seq) == variant.len_diff, variant
        if variant.side is None:
            # only substitute
            return mut_seq
        elif variant.side == 'left':
            return mut_seq[-variant.len_diff:]
        else:
            return mut_seq[:variant.len_diff]

    def _get_del_seq(self, seq, variant, position, fasta):
            # deletion
        if len(variant.REF) + position > len(seq):
            # The actual variant position exceeded the retrieved sequence
            return seq
        if variant.side is None:
            mut_seq = seq[:max(0, position)] + variant.ALT + \
                seq[position + len(variant.REF):]
            assert len(seq) - len(mut_seq) == variant.len_diff, variant
            return mut_seq
        elif variant.side == 'left':
            if self.strand == "+":
                mut_seq = fasta.get_seq(
                    self.chrom,
                    self.start - variant.len_diff, self.end,
                    self.strand == '-').seq
            else:
                mut_seq = fasta.get_seq(
                    self.chrom,
                    self.start, self.end + variant.len_diff,
                    self.strand == '-').seq
            mut_seq = mut_seq[:position + variant.len_diff] + \
                variant.ALT + seq[position + len(variant.REF):]
            assert len(mut_seq) == len(seq), variant
            return mut_seq
        else:
            if self.strand == "+":
                mut_seq = fasta.get_seq(
                    self.chrom,
                    self.start, self.end + variant.len_diff,
                    self.strand == '-').seq
            else:
                mut_seq = fasta.get_seq(
                    self.chrom,
                    self.start - variant.len_diff, self.end,
                    self.strand == '-').seq
            mut_seq = seq[:position] + variant.ALT + \
                mut_seq[position + len(variant.REF):]
            assert len(mut_seq) == len(seq), variant
            return mut_seq

    def get_mut_seq(self, fasta, variant):
        assert variant.side in (None, "left", "right")

        if self._ref_check(fasta, variant):
            warnings.warn(
                "Reference not match, cannot mutate, return original sequence.", UserWarning)
            print(variant)
            return self.get_seq(fasta)

        seq = self.get_seq(fasta)
        position = self._var_pos(variant)

        # position is 1 based, len(seq) 0 based
        if position < 0 or position >= len(seq):
            return seq
        elif variant.len_diff == 0:
            return self._get_snp_seq(seq, variant, position)
        elif variant.len_diff < 0:
            return self._get_insertion_seq(seq, variant, position)
        elif variant.len_diff > 0:
            return self._get_del_seq(seq, variant, position, fasta)

    def _var_pos(self, variant):
        """ Get variant relative position in the sequence
        """
        if self.strand == "+":
            return variant.POS - self.start
        else:
            return self.end - variant.POS - len(variant.REF) + 1


class FastaSeq(Fasta):
    ''' Implement a getSeq method that return upper string and take strand
    '''

    def getSeq(self, iv):
        seq = self.get_seq(iv.chrom, iv.start, iv.end, iv.strand == '-')
        return seq.seq.upper()


@functools.lru_cache(maxsize=1)
def GenerateExonIntervalTree(gtf_file,
                             overhang=(100, 100),  # overhang from the exon
                             gtf_db_path=":memory:",
                             out_file=None,
                             disable_infer_transcripts=True,
                             disable_infer_genes=True,
                             firstLastNoExtend=True,
                             source_filter=None):
    """
    Build IntervalTree object from gtf file for one feature unit (e.g. gene, exon). If give out_file, pickle it.
    Args:
        gtf_file: gtf format file or pickled Intervaltree object.
        overhang: flanking intron length to take along with exon. Corresponding to left (acceptor side) and right (donor side)
        gtf_db_path: (optional) gtf database path. Database for one gtf file only need to be created once
        out_file: (optional) file path to store the pickled Intervaltree obejct. Next time run it can be given to `gtf_file`
        disable_infer_transcripts: option to disable infering transcripts. Can be True if the gtf file has transcripts annotated.
        disable_infer_genes: option to disable infering genes. Can be True if the gtf file has genes annotated.
        firstLastNoExtend: if True, overhang is not taken for 5' of the first exon, or 3' of the last exon of a gene.
        source_filter: gene source filters, such as "protein_coding" filter for protein coding genes
    """
    try:
        gtf_db = gffutils.interface.FeatureDB(gtf_db_path)
    except ValueError:
        gtf_db = gffutils.create_db(
            gtf_file,
            gtf_db_path,
            disable_infer_transcripts=disable_infer_transcripts,
            disable_infer_genes=disable_infer_genes)

    genes = gtf_db.features_of_type('gene')
    exonTree = IntervalTree()
    default_overhang = overhang
    for gene in genes:
        if source_filter is not None:
            if gene.source != source_filter:
                continue
        for exon in gtf_db.children(gene, featuretype='exon'):
            isLast = False  # track whether is last exon
            if firstLastNoExtend:
                if (gene.strand == "+" and exon.end == gene.end) or (gene.strand == "-" and exon.start == gene.start):
                    overhang = (overhang[0], 0)
                    isLast = True
                elif (gene.strand == "+" and exon.start == gene.start) or (gene.strand == "-" and exon.end == gene.end):
                    # int(exon.attributes['exon_number'][0]) == 1:
                    overhang = (0, overhang[1])
            iv = ExonInterval.from_Feature(exon, overhang)
            iv.isLast = isLast
            overhang = default_overhang
            exonTree.insert(iv)
    if out_file is not None:
        with open(out_file, 'wb') as f:
            pickle.dump(exonTree, f)
    return exonTree


class SplicingVCFDataloader(SampleIterator):
    """
    Load genome annotation (gtf) file along with a vcf file, return wt sequence and mut sequence.
    Args:
        gtf: gtf file or pickled gtf IntervalTree. Can be dowloaded from ensembl/gencode. Filter for protein coding genes.
        fasta_file: file path; Genome sequence
        vcf_file: vcf file, each line should contain one and only one variant, left-normalized
        spit_seq: whether or not already split the sequence when loading the data. Otherwise it can be done in the model class.
        endcode: if split sequence, should it be one-hot-encoded
        exon_cut_l: when extract exon feature, how many base pair to cut out at the begining of an exon
        exon_cut_r: when extract exon feature, how many base pair to cut out at the end of an exon
           (cut out the part that is considered as acceptor site or donor site)
        acceptor_intron_cut: how many bp to cut out at the end of acceptor intron that consider as acceptor site
        donor_intron_cut: how many bp to cut out at the end of donor intron that consider as donor site
        acceptor_intron_len: what length in acceptor intron to consider for acceptor site model
        acceptor_exon_len: what length in acceptor exon to consider for acceptor site model
        donor_intron_len: what length in donor intron to consider for donor site model
        donor_exon_len: what length in donor exon to consider for donor site model
        out_file: file path to save pickle file for IntervalTree object derived from GTF annotation. 
            Once the IntervalTree object is generated and saved as a pickle file, next run it can be directly provided to the `gtf` argument,
            the program will save time by not generating it again.
        variant_filter: if set True (default), variants with `FILTER` field other than `PASS` will be filtered out.
        **kwargs: kwargs for `GenerateExonIntervalTree` object
    """

    def __init__(self,
                 gtf,
                 fasta_file,
                 vcf_file=None,
                 split_seq=False,
                 exon_cut_l=0,
                 exon_cut_r=0,
                 acceptor_intron_cut=6,
                 donor_intron_cut=6,
                 acceptor_intron_len=50,
                 acceptor_exon_len=3,
                 donor_exon_len=5,
                 donor_intron_len=13,
                 variant_filter=True,
                 encode=True,
                 **kwargs):
        try:
            with open(gtf, 'rb') as f:
                self.exons = pickle.load(f)
        except (FileNotFoundError, pickle.UnpicklingError):
            self.exons = GenerateExonIntervalTree(gtf, **kwargs)
        import six
        if isinstance(fasta_file, six.string_types):
            fasta = Fasta(fasta_file, as_raw=False)
        self.fasta = fasta
        self.ssGenerator = self.spliceSiteGenerator(vcf_file, self.exons, variant_filter)

        self.encode = encode
        self.split_seq = split_seq
        self.exon_cut_l = exon_cut_l
        self.exon_cut_r = exon_cut_r
        self.acceptor_intron_cut = acceptor_intron_cut
        self.donor_intron_cut = donor_intron_cut
        self.acceptor_intron_len = acceptor_intron_len
        self.acceptor_exon_len = acceptor_exon_len
        self.donor_exon_len = donor_exon_len
        self.donor_intron_len = donor_intron_len

    @staticmethod
    def spliceSiteGenerator(vcf_file, exonTree, variant_filter=True):
        variants = VCF(vcf_file)
        for var in variants:
            if variant_filter and var.FILTER:
                next
            iv = VariantInterval.from_Variant(var)

            matches = map(lambda x: x.interval,
                          exonTree.intersect(iv, ignore_strand=True))

            for match in matches:
                side = get_var_side((
                    var.POS,
                    var.REF,
                    var.ALT,
                    match.Exon_Start,
                    match.Exon_End,
                    match.strand
                ))
                var = iv.to_Variant(match.strand, side)  # to my Variant class
                yield match, var

    def __iter__(self):
        return self

    def __next__(self):
        ss, var = next(self.ssGenerator)
        seq = ss.get_seq(self.fasta).upper()
        mut_seq = ss.get_mut_seq(self.fasta, var).upper()

        if self.split_seq:
            seq = self.split(seq, ss.overhang)
            mut_seq = self.split(mut_seq, ss.overhang)

        return {
            'inputs': {
                'seq': seq,
                'intronl_len': ss.overhang[0],
                'intronr_len': ss.overhang[1],

            },
            'inputs_mut': {
                'seq': mut_seq,
                'intronl_len': ss.overhang[0],
                'intronr_len': ss.overhang[1]
            },
            'metadata': {
                'ranges': ss.grange,
                'variant': var.to_dict,
                'ExonInterval': ss.to_dict,
                'annotation': str(ss)
            }
        }

    def batch_predict_iter(self, **kwargs):
        """Returns samples directly useful for prediction x["inputs"]
        Args:
          **kwargs: Arguments passed to self.batch_iter(**kwargs)
        """
        return (x for x in self.batch_iter(**kwargs))

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

        acceptor_start = intronl_len - self.acceptor_intron_len
        acceptor_end = intronl_len + self.acceptor_exon_len
        acceptor = x[acceptor_start: acceptor_end]

        exon_start = intronl_len + self.exon_cut_l
        exon_end = -intronr_len - self.exon_cut_r
        exon = x[exon_start: exon_end]

        donor_start = -intronr_len - self.donor_exon_len
        donor_end = -intronr_len + self.donor_intron_len
        donor = x[donor_start: donor_end]

        donor_intron = x[-intronr_len + self.donor_intron_cut:]

        if donor[self.donor_exon_len:self.donor_exon_len + 2] != "GT":
            warnings.warn("None GT donor", UserWarning)
        if acceptor[self.acceptor_intron_len - 2:self.acceptor_intron_len] != "AG":
            warnings.warn("None AG donor", UserWarning)

        splits = {
            "acceptor_intron": acceptor_intron,
            "acceptor": acceptor,
            "exon": exon,
            "donor": donor,
            "donor_intron": donor_intron
        }

        if self.encode:
            return {k: encodeDNA([v]) for k, v in splits.items()}

        return splits
