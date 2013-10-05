from astrodendro import load_perseus

def test_load_perseus():
    p = load_perseus()
    assert p.data.shape == (313, 217)
