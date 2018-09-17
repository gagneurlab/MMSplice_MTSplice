import pandas as pd
import numpy as np
from sklearn.externals import joblib
from pkg_resources import resource_filename
LINEAR_MODEL = joblib.load(resource_filename('mmsplice', 'models/linear_model.pkl'))

def max_varEff(df):
    """ Summarize largest absolute effect per variant across all affected exons.
    Args:
        df: result of `predict_all_table`
    """
    if isinstance(df, str):
        df = pd.read_csv(df, index_col=0)
    ref_list = ['mmsplice_ref_acceptorIntron', 'mmsplice_ref_acceptor', 'mmsplice_ref_exon', 'mmsplice_ref_donor', 'mmsplice_ref_donorIntron']
    alt_list = ['mmsplice_alt_acceptorIntron', 'mmsplice_alt_acceptor', 'mmsplice_alt_exon', 'mmsplice_alt_donor', 'mmsplice_alt_donorIntron']
    if 'mmsplice_diff' not in df.columns:
        X = df[alt_list].values - df[ref_list].values
        X = transform(X)
        df['mmsplice_diff'] = LINEAR_MODEL.predict(X)
    dfMax = df.groupby(['ID'], as_index=False).agg({'mmsplice_diff': lambda x: max(x, key=abs)})
    dfMax = dfMax.merge(df, how='left', on = ['ID', 'mmsplice_diff'])
    dfMax = dfMax.drop_duplicates(subset=['ID', 'mmsplice_diff'])
    # dfMax = dfMax.drop("mmsplice_diff", axis=1)
    return dfMax

def not_close0(arr):
    return ~np.isclose(arr, 0)

def transform(X, region_only=False):
    ''' Make interaction terms for the overlapping prediction region
    Args:
        X: modular prediction. Shape (, 5)
        region_only: only interaction terms with indicator function on overlapping
    '''
    exon_overlap = np.logical_or(np.logical_and(not_close0(X[:,1]), not_close0(X[:,2])), np.logical_and(not_close0(X[:,2]), not_close0(X[:,3])))
    acceptor_intron_overlap = np.logical_and(not_close0(X[:,0]), not_close0(X[:,1]))
    donor_intron_overlap = np.logical_and(not_close0(X[:,3]), not_close0(X[:,4]))

    if region_only:
        X = np.hstack([X, exon_overlap.reshape(-1,1)])
        X = np.hstack([X, donor_intron_overlap.reshape(-1,1)]) 
        X = np.hstack([X, acceptor_intron_overlap.reshape(-1,1)])
    else:
        X = np.hstack([X, (X[:,2]*exon_overlap).reshape(-1,1)])
        X = np.hstack([X, (X[:,4]*donor_intron_overlap).reshape(-1,1)]) 
        X = np.hstack([X, (X[:,0]*acceptor_intron_overlap).reshape(-1,1)])
    
    return X