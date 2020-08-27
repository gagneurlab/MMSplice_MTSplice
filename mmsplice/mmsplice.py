import warnings
from pkg_resources import resource_filename
from tqdm import tqdm
import numpy as np
import pandas as pd
from keras.models import load_model
from sklearn.externals import joblib
import concise
from mmsplice.utils import logit, predict_deltaLogitPsi, \
    predict_pathogenicity, predict_splicing_efficiency, encodeDNA, \
    read_ref_psi_annotation, delta_logit_PSI_to_delta_PSI, \
    mmsplice_ref_modules, mmsplice_alt_modules, df_batch_writer
from mmsplice.exon_dataloader import SeqSpliter
from mmsplice.mtsplice import MTSplice, tissue_names
from mmsplice.layers import GlobalAveragePooling1D_Mask0


ACCEPTOR_INTRON = resource_filename('mmsplice', 'models/Intron3.h5')
DONOR = resource_filename('mmsplice', 'models/Donor.h5')
EXON = resource_filename('mmsplice', 'models/Exon.h5')
EXON3 = resource_filename('mmsplice', 'models/Exon_prime3.h5')
ACCEPTOR = resource_filename('mmsplice', 'models/Acceptor.h5')
DONOR_INTRON = resource_filename('mmsplice', 'models/Intron5.h5')
LINEAR_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/linear_model.pkl'))
LOGISTIC_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/Pathogenicity.pkl'))
EFFICIENCY_MODEL = joblib.load(resource_filename(
    'mmsplice', 'models/splicing_efficiency.pkl'))


class MMSplice(object):
    """
    Load modules of mmsplice model, perform prediction on batch of dataloader.

    Args:
      acceptor_intronM: acceptor intron model,
        score acceptor intron sequence.
      acceptorM: accetpor splice site model. Score acceptor sequence
        with 50bp from intron, 3bp from exon.
      exonM: exon model, score exon sequence.
      donorM: donor splice site model, score donor sequence
        with 13bp in the intron, 5bp in the exon.
      donor_intronM: donor intron model, score donor intron sequence.
    """

    def __init__(self,
                 acceptor_intronM=ACCEPTOR_INTRON,
                 acceptorM=ACCEPTOR,
                 exonM=EXON,
                 donorM=DONOR,
                 donor_intronM=DONOR_INTRON,
                 seq_spliter=None):
        self.spliter = seq_spliter or SeqSpliter()
        self.acceptor_intronM = load_model(acceptor_intronM, compile=False)
        self.acceptorM = load_model(acceptorM, compile=False)
        self.exonM = load_model(exonM,
                                custom_objects={"GlobalAveragePooling1D_Mask0":
                                                GlobalAveragePooling1D_Mask0},
                                compile=False)
        self.donorM = load_model(donorM, compile=False)
        self.donor_intronM = load_model(donor_intronM, compile=False)

    def predict_on_batch(self, batch):
        warnings.warn(
            "`self.predict_on_batch` is deprecated,"
            " use `self.predict_modular_scores_on_batch instead`",
            DeprecationWarning
        )
        return self.predict_modular_scores_on_batch(batch)

    def predict_modular_scores_on_batch(self, batch):
        '''
        Perform prediction on batch of dataloader.

        Args:
          batch: batch of dataloader.

        Returns:
          np.matrix of modular predictions
          as [[acceptor_intronM, acceptor, exon, donor, donor_intron]]

        '''
        score = np.concatenate([
            self.acceptor_intronM.predict(batch['acceptor_intron']),
            logit(self.acceptorM.predict(batch['acceptor'])),
            self.exonM.predict(batch['exon']),
            logit(self.donorM.predict(batch['donor'])),
            self.donor_intronM.predict(batch['donor_intron'])
        ], axis=1)
        return score

    def predict(self, *args, **kwargs):
        warnings.warn(
            "self.predict is deprecated, use self.predict_on_seq instead",
            DeprecationWarning
        )
        return self.predict_on_seq(*args, **kwargs)

    def predict_on_seq(self, seq, overhang=(100, 100)):
        """
        Performe prediction of overhanged exon sequence string.

        Args:
          seq (str):  sequence of overhanged exon.
          overhang (Tuple[int, int]): overhang of seqeunce.

        Returns:
          np.array of modular predictions
          as [[acceptor_intronM, acceptor, exon, donor, donor_intron]].
        """
        batch = self.spliter.split(seq, overhang)
        batch = {k: encodeDNA([v]) for k, v in batch.items()}
        return self.predict_modular_scores_on_batch(batch)[0]

    def _predict_batch(self, batch, optional_metadata=None):
        optional_metadata = optional_metadata or []

        X_ref = self.predict_modular_scores_on_batch(
            batch['inputs']['seq'])
        X_alt = self.predict_modular_scores_on_batch(
            batch['inputs']['mut_seq'])
        ref_pred = pd.DataFrame(X_ref, columns=mmsplice_ref_modules)
        alt_pred = pd.DataFrame(X_alt, columns=mmsplice_alt_modules)

        df = pd.DataFrame({
            'ID': batch['metadata']['variant']['annotation'],
            'exons': batch['metadata']['exon']['annotation'],
        })

        for key in optional_metadata:
            for k, v in batch['metadata'].items():
                if key in v:
                    df[key] = v[key]

        df['delta_logit_psi'] = predict_deltaLogitPsi(X_ref, X_alt)
        df = pd.concat([df, ref_pred, alt_pred], axis=1)
        return df

    def _predict_batch_mtsplice(self, batch, df, mtsplice,
                                natural_scale, df_ref):
        X_tissue = mtsplice.predict_on_batch(
            batch['inputs']['tissue_seq'])
        X_tissue += np.expand_dims(
            df['delta_logit_psi'].values, axis=1)
        tissue_pred = pd.DataFrame(X_tissue, columns=tissue_names)
        df = pd.concat([df, tissue_pred], axis=1)

        if natural_scale:
            df_ref = df_ref[df_ref.columns[6:]]
            df = df.join(df_ref, on='exons', rsuffix='_ref')

            ref_tissue_names = ['%s_ref' % i for i in df_ref.columns]

            delta_psi_pred = pd.DataFrame(
                delta_logit_PSI_to_delta_PSI(
                    df[df_ref.columns].values,
                    df[ref_tissue_names].values
                ), columns=['%s_delta_psi' % i
                            for i in df_ref.columns])
            df.reset_index(drop=True, inplace=True)
            delta_psi_pred.reset_index(drop=True, inplace=True)
            df = pd.concat([df, delta_psi_pred], axis=1)
        return df

    def _predict_on_dataloader(self, dataloader, batch_size=512, progress=True,
                               pathogenicity=False, splicing_efficiency=False,
                               natural_scale=False, ref_psi_version=None):
        """
        Make prediction from a dataloader, return results as a table

        Args:
           model: mmsplice model object.
           dataloader: dataloader object.
           progress: show progress bar.
           pathogenicity: adds pathogenicity prediction as column
           splicing_efficiency: adds splicing_efficiency prediction as column

        Returns:
           iterator of pd.DataFrame includes modular prediction,
             delta_logit_psi, splicing_efficiency, pathogenicity.
        """
        from mmsplice.exon_dataloader import ExonSplicingMixin
        assert isinstance(dataloader, ExonSplicingMixin), \
            "Unknown dataloader type"

        if dataloader.tissue_specific:
            mtsplice = MTSplice()
            if natural_scale:
                df_ref = read_ref_psi_annotation(
                    ref_psi_version, set(dataloader.vcf.seqnames))
            else:
                df_ref = None

        dt_iter = dataloader.batch_iter(batch_size=batch_size)
        if progress:
            dt_iter = tqdm(dt_iter)

        for batch in dt_iter:
            df = self._predict_batch(
                batch, dataloader.optional_metadata)
            X_ref = df[mmsplice_ref_modules].values
            X_alt = df[mmsplice_alt_modules].values

            if dataloader.tissue_specific:
                df = self._predict_batch_mtsplice(
                    batch, df, mtsplice, natural_scale, df_ref)

            if pathogenicity:
                df['pathogenicity'] = predict_pathogenicity(
                    X_ref, X_alt)
            if splicing_efficiency:
                df['efficiency'] = predict_splicing_efficiency(
                    X_ref, X_alt)

            yield df

    def predict_on_dataloader(self, dataloader, batch_size=512, progress=True,
                              pathogenicity=False, splicing_efficiency=False,
                              natural_scale=False, ref_psi_version=None):
        """Make prediction from a dataloader, return results as a table
        Args:
           model: mmsplice model object.
           dataloader: dataloader object.
           progress: show progress bar.
           pathogenicity: adds pathogenicity prediction as column
           splicing_efficiency: adds splicing_efficiency prediction as column

        Returns:
           pd.DataFrame includes modular prediction, delta_logit_psi,
           splicing_efficiency, pathogenicity.

        """
        return pd.concat(
            self._predict_on_dataloader(
                dataloader,
                batch_size=batch_size,
                progress=progress,
                pathogenicity=pathogenicity,
                splicing_efficiency=splicing_efficiency,
                natural_scale=natural_scale, ref_psi_version=ref_psi_version)
        )


# TODO: implement prediction methods within MMSplice class,
#   should be more error prone
def predict_save(model, dataloader, output_csv, batch_size=512, progress=True,
                 pathogenicity=False, splicing_efficiency=False):
    from mmsplice import MMSplice
    assert isinstance(model, MMSplice), \
        "model should be a mmsplice.MMSplice class instance"

    df_iter = model._predict_on_dataloader(
        dataloader, progress=progress,
        pathogenicity=pathogenicity,
        splicing_efficiency=splicing_efficiency)

    return df_batch_writer(df_iter, output_csv)


def predict_all_table(model, dataloader, batch_size=512, progress=True,
                      pathogenicity=False, splicing_efficiency=False,
                      natural_scale=False, ref_psi_version=None):
    """
    Return the prediction as a table

    Args:
      model: mmsplice model object.
      dataloader: dataloader object.
      progress: show progress bar.
      pathogenicity: adds pathogenicity prediction as column
      splicing_efficiency: adds  splicing_efficiency prediction as column

    Returns:
      pd.DataFrame of modular prediction, delta_logit_psi, splicing_efficiency,
        pathogenicity.
    """
    from mmsplice import MMSplice
    assert isinstance(model, MMSplice), \
        "model should be a mmsplice.MMSplice class instance"

    return model.predict_on_dataloader(
        dataloader, progress=progress, batch_size=batch_size,
        pathogenicity=pathogenicity, splicing_efficiency=splicing_efficiency,
        natural_scale=natural_scale, ref_psi_version=ref_psi_version)


def writeVCF(vcf_in, vcf_out, predictions):
    from cyvcf2 import Writer, VCF
    with VCF(vcf_in) as vcf:
        vcf.add_info_to_header({
            'ID': 'mmsplice',
            'Description': 'mmsplice splice variant effect',
            'Type': 'Character',
            'Number': '.'
        })

        with Writer(vcf_out, vcf) as w:
            for var in vcf:
                pred = predictions.get(var.ID)
                if pred is not None:
                    var.INFO['mmsplice'] = pred
                w.write_record(var)
