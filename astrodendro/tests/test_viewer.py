import pytest
import numpy as np

import matplotlib.pyplot as plt
from ..dendrogram import Dendrogram
from matplotlib.backend_bases import MouseEvent


DATA = np.array([[1, 3, 4, 4, 1, 4],
                 [1, 2, 3, 2, 1, 3],
                 [2, 1, 1, 3, 1, 2],
                 [1, 1, 1, 1, 1, 1],
                 [2, 3, 2, 1, 1, 2],
                 [2, 3, 5, 3, 1, 1]])


def test_viewer(capsys):

    original_backend = plt.get_backend()

    try:
        plt.switch_backend('qtagg')
    except ImportError:
        pytest.skip("This test requires Qt to be installed")

    d = Dendrogram.compute(DATA)
    viewer = d.viewer()

    plt.show(block=False)

    cb = viewer.fig.canvas.callbacks

    cb.process('button_press_event', MouseEvent('button_press_event', viewer.fig.canvas, 660, 520, 1))
    cb.process('button_press_event', MouseEvent('button_press_event', viewer.fig.canvas, 890, 800, 1))
    cb.process('button_press_event', MouseEvent('button_press_event', viewer.fig.canvas, 700, 700, 1))

    plt.switch_backend(original_backend)

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
