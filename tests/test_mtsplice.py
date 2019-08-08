from mmsplice import MTSplice


def test_mtsplice():
    seq = 'ATGCGACGTACCCAGTAAAT'
    overhang = (4, 4)
    model = MTSplice()
    pred = model.predict(seq, overhang)
    assert pred.shape == (56,)
