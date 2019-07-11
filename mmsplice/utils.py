from collections import namedtuple
import pandas as pd
import numpy as np
import pyranges
from sklearn.externals import joblib
from pkg_resources import resource_filename

LINEAR_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/linear_model.pkl'))
LOGISTIC_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/Pathogenicity.pkl'))
EFFICIENCY_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/splicing_efficiency.pkl'))


class Variant(namedtuple('Variant', ['CHROM', 'POS', 'REF', 'ALT'])):

    @property
    def start(self):
        return self.POS - 1


def clip(x):
    return np.clip(x, 0.00001, 0.99999)


def logit(x):
    x = clip(x)
    return np.log(x) - np.log(1 - x)


def expit(x):
    return 1. / (1. + np.exp(-x))


def pyrange_remove_chr_from_chrom_annotation(pr):
    df = pr.df
    df['Chromosome'] = df['Chromosome'].str.replace('chr', '')
    return pyranges.PyRanges(df)


def max_varEff(df):
    """ Summarize largest absolute effect per variant across all affected exons.
    Args:
        df: result of `predict_all_table`
    """
    if isinstance(df, str):
        df = pd.read_csv(df, index_col=0)

    df_max = df.groupby(['ID'], as_index=False).agg(
        {'delta_logit_psi': lambda x: max(x, key=abs)})

    df_max = df_max.merge(df, how='left', on=['ID', 'delta_logit_psi'])
    df_max = df_max.drop_duplicates(subset=['ID', 'delta_logit_psi'])
    return df_max


def _not_close0(arr):
    return ~np.isclose(arr, 0)


def _and_not_close0(x, y):
    return np.logical_and(_not_close0(x), _not_close0(y))


def transform(X, region_only=False):
    ''' Make interaction terms for the overlapping prediction region
    Args:
        X: modular prediction. Shape (, 5)
        region_only: only interaction terms with indicator function on overlapping
    '''
    exon_overlap = np.logical_or(
        _and_not_close0(X[:, 1], X[:, 2]),
        _and_not_close0(X[:, 2], X[:, 3])
    )
    acceptor_intron_overlap = _and_not_close0(X[:, 0], X[:, 1])
    donor_intron_overlap = _and_not_close0(X[:, 3], X[:, 4])

    if not region_only:
        exon_overlap = X[:, 2] * exon_overlap
        donor_intron_overlap = X[:, 4] * donor_intron_overlap
        acceptor_intron_overlap = X[:, 0] * acceptor_intron_overlap

    return np.hstack([
        X,
        exon_overlap.reshape(-1, 1),
        donor_intron_overlap.reshape(-1, 1),
        acceptor_intron_overlap.reshape(-1, 1)
    ])


def predict_delta_logit_psi(X_ref, X_alt):
    return LINEAR_MODEL.predict(transform(X_alt - X_ref, region_only=False))


def predict_pathogenicity(X_ref, X_alt):
    X = transform(X_alt - X_ref, region_only=True)
    X = np.concatenate([X_ref, X_alt, X[:, -3:]], axis=-1)
    return LOGISTIC_MODEL.predict_proba(X)[:, 1]


def predict_splicing_efficiency(X_ref, X_alt):
    X = transform(X_alt - X_ref, region_only=False)
    X = X[:, [1, 2, 3, 5]]  # no intronic modules
    return EFFICIENCY_MODEL.predict(X)


def read_vep(vep_result_path,
             max_per_var=False):
    ''' Read MMSplice VEP plugin output. Only support vcf type output.

    Args:
        vep_result_path: file path to the returned result of VEP plugin.
        max_per_var: return maximum absolute effect size per variant.
    '''

    from cyvcf2 import VCF
    from collections import defaultdict

    score_pred = []

    keys = [
        'mmsplice_alt_acceptor',
        'mmsplice_alt_acceptorIntron',
        'mmsplice_alt_donor',
        'mmsplice_alt_donorIntron',
        'mmsplice_alt_exon',
        'mmsplice_delta_logit_psi',
        'mmsplice_pathogenicity',
        'mmsplice_ref_acceptor',
        'mmsplice_ref_acceptorIntron',
        'mmsplice_ref_donor',
        'mmsplice_ref_donorIntron',
        'mmsplice_ref_exon'
    ]

    alt_seqs = defaultdict(list)
    ref_seqs = defaultdict(list)

    for l in VCF(vep_result_path):
        csq = l.INFO['CSQ'].split(',')
        predictions = map(lambda x: tuple(x.split('|')[-len(keys):]), csq)

        for pred in predictions:
            if pred != ('',) * len(keys):
                x = dict(
                    zip(keys, map(float, (i if i != '' else 0 for i in pred))))
                x['ID'] = "%s:%d:%s:%s" % (
                    l.CHROM, int(l.start) + 1, l.REF, l.ALT)
                score_pred.append(x)

    df_plugin = pd.DataFrame(score_pred)

    if max_per_var:
        df_plugin = max_varEff(df_plugin).set_index('ID')

    return df_plugin
