from mmsplice.generic import get_var_side


def test_get_var_side():
    var = (9, 'A', 'AGG', 11, 20, '+')
    assert get_var_side(var) == 'left'

    var = (10, 'A', 'AGG', 11, 20, '+')
    assert get_var_side(var) is "exon"

    var = (19, 'A', 'AGG', 11, 20, '+')
    assert get_var_side(var) is 'right'

    var = (12, 'A', 'AGG', 11, 20, '+')
    assert get_var_side(var) is "exon"

    var = (9, 'A', 'AGG', 11, 20, '-')
    assert get_var_side(var) == 'right'

    var = (10, 'A', 'AGG', 11, 20, '-')
    assert get_var_side(var) is "exon"

    var = (19, 'A', 'AGG', 11, 20, '-')
    assert get_var_side(var) == 'left'

    var = (12, 'A', 'AGG', 11, 20, '-')
    assert get_var_side(var) is "exon"
