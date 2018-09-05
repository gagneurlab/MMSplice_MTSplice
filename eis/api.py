import uuid

from flask import Flask, jsonify, request

from eis import Eis

app = Flask(__name__)

psi_models = dict()


@app.route('/create-model', methods=['POST'])
def create_model():
    global psi_models
    model_id = str(uuid.uuid4())

    psi_models[model_id] = Eis(
        exon_cut_l=int(request.json['exon_cut_l']),
        exon_cut_r=int(request.json['exon_cut_r']),
        acceptor_intron_cut=int(request.json['acceptor_intron_cut']),
        donor_intron_cut=int(request.json['donor_intron_cut']),
        acceptor_intron_len=int(request.json['acceptor_intron_len']),
        acceptor_exon_len=int(request.json['acceptor_exon_len']),
        donor_exon_len=int(request.json['donor_exon_len']),
        donor_intron_len=int(request.json['donor_intron_len'])
    )

    # warms up the model

    x = {'seq': "A" * 100,
         'intronl_len': 4,
         'intronr_len': 4}

    psi_models[model_id].predict(x)

    return jsonify(model_id)


@app.route('/psi-score', methods=['POST'])
def psi_score():
    global psi_models

    model_id = request.json['model_id']
    psi_model = psi_models[model_id]

    x = {
        'seq': request.json['seq'],
        'intronl_len': int(request.json['intronl_len']),
        'intronr_len': int(request.json['intronr_len'])
    }

    scores = psi_model.predict(x)

    return ','.join(map(str, scores.tolist()))
