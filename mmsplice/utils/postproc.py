
import pandas as pd
import numpy as np
from sklearn.externals import joblib
from pkg_resources import resource_filename
LINEAR_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/linear_model.pkl'))
LOGISTIC_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/Pathogenicity.pkl'))


def max_varEff(df):
    """ Summarize largest absolute effect per variant across all affected exons.
    Args:
        df: result of `predict_all_table`
    """
    if isinstance(df, str):
        df = pd.read_csv(df, index_col=0)

    ref_list = ['mmsplice_ref_acceptorIntron',
                'mmsplice_ref_acceptor',
                'mmsplice_ref_exon',
                'mmsplice_ref_donor',
                'mmsplice_ref_donorIntron']
    alt_list = ['mmsplice_alt_acceptorIntron',
                'mmsplice_alt_acceptor',
                'mmsplice_alt_exon',
                'mmsplice_alt_donor',
                'mmsplice_alt_donorIntron']

    if 'mmsplice_dlogitPsi' not in df.columns:
        X = df[alt_list].values - df[ref_list].values
        X = transform(X)
        df['mmsplice_dlogitPsi'] = LINEAR_MODEL.predict(X)

    dfMax = df.groupby(['ID'], as_index=False).agg(
        {'mmsplice_dlogitPsi': lambda x: max(x, key=abs)})

    dfMax = dfMax.merge(df, how='left', on=['ID', 'mmsplice_dlogitPsi'])
    dfMax = dfMax.drop_duplicates(subset=['ID', 'mmsplice_dlogitPsi'])
    # dfMax = dfMax.drop("mmsplice_dlogitPsi", axis=1)
    return dfMax


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


def predict_deltaLogitPsi(X_ref, X_alt):
    return LINEAR_MODEL.predict(transform(X_alt - X_ref, region_only=False))


def predict_pathogenicity(X_ref, X_alt):
    X = transform(X_alt - X_ref, region_only=True)
    X = np.concatenate([X_ref, X_alt, X[:, -3:]], axis=-1)
    return LOGISTIC_MODEL.predict(X)
