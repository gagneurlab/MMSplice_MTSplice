import pandas as pd
import numpy as np

def max_varEff(df, combine_fn=lambda x: np.sum(x, axis=1)):
	""" Summarize largest absolute effect per variant across all affected exons.
	Args:
		df: result of `predict_all_table`
		combine_fn: Maximum effect size is calculated by calculating 
		the combined effect size of all 5 modules.
	"""
	ref_list = ['EIS_ref_acceptorIntron', 'EIS_ref_acceptor', 'EIS_ref_exon', 'EIS_ref_donor', 'EIS_ref_donorIntron']
	alt_list = ['EIS_alt_acceptorIntron', 'EIS_alt_acceptor', 'EIS_alt_exon', 'EIS_alt_donor', 'EIS_alt_donorIntron']
	EIS_diff = df[alt_list].as_matrix() - df[ref_list].as_matrix()
	df['EIS_diff'] = combine_fn(EIS_diff)
	dfMax = df.groupby(['ID'], as_index=False).agg({'EIS_diff': lambda x: max(x, key=abs)})
	dfMax = dfMax.merge(df, how='left', on = ['ID', 'EIS_diff'])
	dfMax = dfMax.drop_duplicates(subset=['ID', 'EIS_diff'])
	# dfMax = dfMax.drop("EIS_diff", axis=1)
	return dfMax