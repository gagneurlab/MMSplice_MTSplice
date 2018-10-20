import sys
import json

import click
import numpy as np

from keras import backend as K
from .mmsplice import MMSplice
from .utils import predict_deltaLogitPsi, predict_pathogenicity


@click.group()
def cli():
    pass


@cli.command(name='run')
def run():
    options = json.loads(sys.stdin.readline().strip())

    K.clear_session()
    psi_model = MMSplice(
        **{k: v for k, v in options.items() if v},
        pattern_warning=False)

    # warms up the model
    psi_model.predict({
        'seq': "A" * 100,
        'intronl_len': 4,
        'intronr_len': 4
    })

    sys.stdout.write('MMSPLICE-RESPONSE:' + '1\n')
    sys.stdout.flush()

    while True:
        variant = json.loads(sys.stdin.readline().strip())

        ref_scores, alt_scores = [
            np.matrix(
                psi_model.predict({
                    "intronl_len": variant['intronl_len'],
                    "intronr_len": variant['intronr_len'],
                    "seq": variant[seq]
                }).values
            )
            for seq in ['ref_seq', 'alt_seq']
        ]

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
