import os.path
import pytest
from mmsplice.vcf_dataloader import SplicingVCFDataloader
from mmsplice import MMSplice, predict_all_table
from mmsplice.utils import read_vep, max_varEff
from scipy.stats import pearsonr


vep_output = 'variant_effect_output.txt'

@pytest.mark.skipif(not os.path.isfile(vep_output),
                    reason="No vep result file to test")
def test_vep_plugin():
    gtf = 'tests/data/test.gtf'
    vcf = 'tests/data/test.vcf.gz'
    fasta = 'tests/data/hg19.nochr.chr17.fa'

    dl = SplicingVCFDataloader(gtf, fasta, vcf)
    model = MMSplice()

    df_python = predict_all_table(
        model, dl, pathogenicity=True, splicing_efficiency=True)
    df_python_predictionsMax = max_varEff(df_python).set_index('ID')

    df_plugin = read_vep(vep_output)
    df_plugin_predictionsMax = max_varEff(df_plugin).set_index('ID')

    indexes = list(set(df_plugin_predictionsMax.index) &
                   set(df_python_predictionsMax.index))

    vep_plugin_dlogitPsi = df_plugin_predictionsMax.loc[
        indexes, 'delta_logit_psi']
    python_package = df_python_predictionsMax.loc[
        indexes, 'delta_logit_psi']

    assert pearsonr(vep_plugin_dlogitPsi, python_package)[0] >= 0.95
