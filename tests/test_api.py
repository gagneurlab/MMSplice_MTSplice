from mmsplice.api import app

def test_api():
    with app.test_client() as c:

        resp = c.post('/create-model', json={})
        assert resp.status_code == 200

        resp = c.post('/psi-score', json={
            'ref_seq': "ATGCGACGTACCCAGTAAAT",
            'alt_seq': "TTGCGACGTACCCAGTAAAT",
            'intronl_len': 4,
            'intronr_len': 4
        })
        assert len(resp.data.decode("utf-8").split(',')) == 12
