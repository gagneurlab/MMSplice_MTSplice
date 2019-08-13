from pkg_resources import resource_filename
from tqdm import tqdm
import numpy as np
import pandas as pd
from keras import backend as K
from keras.models import load_model
from sklearn.externals import joblib
from concise.preprocessing import encodeDNA

from mmsplice.utils import logit, predict_deltaLogitPsi, \
    predict_pathogenicity, predict_splicing_efficiency
from mmsplice.exon_dataloader import SeqSpliter
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

        K.clear_session()
        self.acceptor_intronM = load_model(acceptor_intronM, compile=False)
        self.acceptorM = load_model(acceptorM, compile=False)
        self.exonM = load_model(exonM,
                                custom_objects={"GlobalAveragePooling1D_Mask0":
                                                GlobalAveragePooling1D_Mask0},
                                compile=False)
        self.donorM = load_model(donorM, compile=False)
        self.donor_intronM = load_model(donor_intronM, compile=False)

    def predict_on_batch(self, batch):
        '''
        Performe prediction on batch of dataloader.

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

    def predict(self, seq, overhang=(100, 100)):
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
        return self.predict_on_batch(batch)[0]


def predict_batch(model, dataloader, batch_size=512, progress=True,
                  pathogenicity=False, splicing_efficiency=False):
    """
    Return the prediction as a table

    Args:
      model: mmsplice model object.
      dataloader: dataloader object.
      progress: show progress bar.
      pathogenicity: adds pathogenicity prediction as column
      splicing_efficiency: adds  splicing_efficiency prediction as column

    Returns:
      iterator of pd.DataFrame of modular prediction, delta_logit_psi,
        splicing_efficiency, pathogenicity.
    """
    dt_iter = dataloader.batch_iter(batch_size=batch_size)
    if progress:
        dt_iter = tqdm(dt_iter)

    ref_cols = ['ref_acceptorIntron', 'ref_acceptor',
                'ref_exon', 'ref_donor', 'ref_donorIntron']
    alt_cols = ['alt_acceptorIntron', 'alt_acceptor',
                'alt_exon', 'alt_donor', 'alt_donorIntron']

    for batch in dt_iter:
        X_ref = model.predict_on_batch(batch['inputs']['seq'])
        X_alt = model.predict_on_batch(batch['inputs']['mut_seq'])
        ref_pred = pd.DataFrame(X_ref, columns=ref_cols)
        alt_pred = pd.DataFrame(X_alt, columns=alt_cols)

        df = pd.DataFrame({
            'ID': batch['metadata']['variant']['STR'],
            'exons': batch['metadata']['exon']['annotation'],
        })
        for k in ['exon_id', 'gene_id', 'gene_name', 'transcript_id']:
            if k in batch['metadata']['exon']:
                df[k] = batch['metadata']['exon'][k]

        df['delta_logit_psi'] = predict_deltaLogitPsi(X_ref, X_alt)
        df = pd.concat([df, ref_pred, alt_pred], axis=1)

        if pathogenicity:
            df['pathogenicity'] = predict_pathogenicity(X_ref, X_alt)
        if splicing_efficiency:
            df['efficiency'] = predict_splicing_efficiency(X_ref, X_alt)

        yield df


def predict_save(model, dataloader, output_csv, batch_size=512, progress=True,
                 pathogenicity=False, splicing_efficiency=False):
    df_iter = predict_batch(model, dataloader, progress=progress,
                            pathogenicity=pathogenicity,
                            splicing_efficiency=splicing_efficiency)

    df = next(df_iter)
    with open(output_csv, 'w') as f:
        df.to_csv(f, index=False)

    for df in df_iter:
        with open(output_csv, 'a') as f:
            df.to_csv(f, index=False, header=False)


def predict_all_table(model,
                      dataloader,
                      batch_size=512,
                      progress=True,
                      pathogenicity=False,
                      splicing_efficiency=False):
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
    return pd.concat(predict_batch(model, dataloader, progress=progress,
                                   pathogenicity=pathogenicity,
                                   splicing_efficiency=splicing_efficiency))


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
