import tempfile
import pytest

gtf_file = 'tests/data/test.gtf'
fasta_file = 'tests/data/hg19.nochr.chr17.fa'
vcf_file = 'tests/data/test.vcf.gz'
exon_file = 'tests/data/test_exons.csv'
junction_file = 'tests/data/test_junction.csv'
junction_psi5_file = 'tests/data/test_junction_psi5_file.csv'
junction_psi3_file = 'tests/data/test_junction_psi3_file.csv'
multi_vcf = 'tests/data/multi_test.vcf.gz'


snps = [
    "17:41276033:C:['G']"
]

deletions = [
    "13:32953886:GTT:['AA']",
    "17:41267740:TACTT:['A']",
    "17:41267741:ACTT:['AC']",
    "17:41267742:CTTGCAAAATATGTGGTCACACTTTGTGGAGACAGGTTCCTTGATCAACTCCAGA:['C']",
    "17:41267742:CTT:['C']",
    "17:41267795:GAC:['GA']",
    "17:41267795:GA:['G']",
    "17:41267796:ACT:['A']",
]

insertions = [
    "17:41267795:G:['GAA']",
    "17:41267797:CT:['CTAA']",
    "17:41276032:AC:['ACCA']",
    "17:41276033:C:['CCAGATG']",
    "17:41276132:A:['ACT']"
]

variants = [
    *snps,
    *deletions,
    *insertions
]

outside = [
    "17:41275833:A:['G']"
]


def parse_vcf_id(vcf_id):
    return vcf_id.replace("'", '').replace('[', '').replace(']', '').split(':')


@pytest.fixture
def vcf_path():
    with tempfile.NamedTemporaryFile('w') as temp_vcf:
        temp_vcf.write('##fileformat=VCFv4.0\n')
        temp_vcf.write('##contig=<ID=13,length=115169878>\n')
        temp_vcf.write('##contig=<ID=17,length=81195210>\n')
        temp_vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        for v in variants:
            temp_vcf.write('%s\t%s\t1\t%s\t%s\t.\t.\t.\n'
                           % tuple(parse_vcf_id(v)))

        temp_vcf.flush()
        yield temp_vcf.name
