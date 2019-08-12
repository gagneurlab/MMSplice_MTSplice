import sys
import json

import click
import numpy as np
from keras import backend as K

from mmsplice import MMSplice
from mmsplice.exon_dataloader import SeqSpliter
from mmsplice.utils import predict_deltaLogitPsi, predict_pathogenicity


@click.group()
def cli():
    pass


@cli.command(name='run')
def run():
    options = json.loads(sys.stdin.readline().strip())

    K.clear_session()
    psi_model = MMSplice(
        **{k: v for k, v in options.items() if v})
    psi_model.spliter = SeqSpliter(pattern_warning=False)

    # warms up the model
    psi_model.predict("A" * 100, (4, 4))

    sys.stdout.write('MMSPLICE-RESPONSE:' + '1\n')
    sys.stdout.flush()

    while True:
        variant = json.loads(sys.stdin.readline().strip())

        overhang = (variant['intronl_len'], variant['intronr_len'])
        ref_scores = np.matrix(psi_model.predict(variant['ref_seq'], overhang))
        alt_scores = np.matrix(psi_model.predict(variant['alt_seq'], overhang))
        scores = np.hstack([ref_scores, alt_scores]).tolist()[0]

        scores.extend([
            predict_deltaLogitPsi(ref_scores, alt_scores)[0],
            predict_pathogenicity(ref_scores, alt_scores)[0]
        ])

        sys.stdout.write('MMSPLICE-RESPONSE:' +
                         ','.join(map(str, scores)) + '\n')
        sys.stdout.flush()


if __name__ == '__main__':
    cli()
