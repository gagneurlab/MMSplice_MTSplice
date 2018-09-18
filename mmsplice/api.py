from flask import Flask, jsonify, request
from keras import backend as K

from .mmsplice import MMSplice

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

    scores = psi_model.predict(request.json)

    return ','.join(map(str, scores.tolist()))
