import numpy as np
from flask import Flask, jsonify, request
from keras import backend as K

from .mmsplice import MMSplice
from .utils import predict_deltaLogitPsi, predict_pathogenicity

app = Flask(__name__)

psi_model = None


@app.route('/create-model', methods=['POST'])
def create_model():
    global psi_model

    K.clear_session()
    psi_model = MMSplice(**{
        k: v
        for k, v in request.json.items()
        if v
    })

    # warms up the model
    psi_model.predict({
        'seq': "A" * 100,
        'intronl_len': 4,
        'intronr_len': 4
    })

    return jsonify(success=True)


@app.route('/psi-score', methods=['POST'])
def psi_score():
    global psi_model

    ref_scores, alt_scores = [
        np.matrix(
            psi_model.predict({
                "intronl_len": request.json['intronl_len'],
                "intronr_len": request.json['intronr_len'],
                "seq": request.json[seq]
            }).values
        )
        for seq in ['ref_seq', 'alt_seq']
    ]

    scores = np.hstack([ref_scores, alt_scores]).tolist()[0]

    scores.extend([
        predict_deltaLogitPsi(ref_scores, alt_scores)[0],
        int(predict_pathogenicity(ref_scores, alt_scores)[0])
    ])

    return ','.join(map(str, scores))
