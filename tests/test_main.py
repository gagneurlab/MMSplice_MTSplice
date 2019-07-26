import json
from subprocess import Popen, PIPE


def test_cli():
    process = Popen(['mmsplice', 'run'], stdin=PIPE, stdout=PIPE)

    options = {}
    process.stdin.write((json.dumps(options) + '\n').encode())
    process.stdin.flush()
    assert process.stdout.readline().decode() == 'MMSPLICE-RESPONSE:1\n'

    variant = {
        'intronl_len': 4,
        'intronr_len': 4,
        'ref_seq': 'A' * 100,
        'alt_seq': 'T' * 100
    }
    process.stdin.write((json.dumps(variant) + '\n').encode())
    process.stdin.flush()
    out = process.stdout.readline().decode().strip()
    process.terminate()

    head, rest = out.split(':')
    assert head == 'MMSPLICE-RESPONSE'
    pred = list(map(float, rest.split(',')))

    assert len(pred) == 12
    assert pred[10] != 0
