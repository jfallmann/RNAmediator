import os
import pytest
import random
import numpy as np
from RIssmed.RNAtweaks.RNAtweaks import api_rnaplfold, cmd_rnaplfold

TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
PARPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_DIR = os.path.dirname(os.path.dirname(PARPATH))


def random_sequence(seed: int = 1):
    random.seed(seed)
    return "".join(random.choices(["A", "C", "U", "G"], k=350))


@pytest.mark.parametrize(
    "seq,window,span,constraint",
    [
        ("A" * 500, 70, 70, []),
        ("A" * 500, 70, 70, [("p", 20, 30)]),
        (random_sequence(), 70, 70, [("p", 20, 30), ("u", 70, 90)]),
        (random_sequence(), 70, 70, [("u", 20, 30)]),
        (random_sequence(), 70, 70, []),
    ],
)
def test_api_and_cmd_plfold(seq, window, span, constraint):
    cmd_result = cmd_rnaplfold(seq, window, span, constraint=constraint)
    api_result = api_rnaplfold(seq, window, span, constraint=constraint)
    assert api_result == cmd_result


@pytest.mark.parametrize(
    "seq,window,span,constraint",
    [
        (random_sequence(5), 70, 70, [("nonesense", 3, 5)]),
    ],
)
def test_constraint_error(seq, window, span, constraint):
    with pytest.raises(ValueError):
        api_rnaplfold(seq, window, span, constraint=constraint)
    with pytest.raises(ValueError):
        cmd_rnaplfold(seq, window, span, constraint=constraint)


def test_localization():
    result = api_rnaplfold(random_sequence(), 70, 70, region=7, constraint=[("p", 20, 30)])
    result_array = result.numpy_array
    result.localize(30, 50)
    test_array = result.numpy_array
    assert np.array_equal(result_array[30:50], test_array, equal_nan=True)


@pytest.mark.parametrize(
    "seq,window,span,constraint",
    [
        (random_sequence(4), 50, 50, ("p", 200, 207)),
        (random_sequence(4), 50, 50, ("u", 200, 207)),
    ],
)
def test_constraint(seq, window, span, constraint):
    mode = constraint[0]
    api_result = cmd_rnaplfold(seq, window, span, region=7, constraint=[constraint])
    api_result.localize(constraint[1], constraint[2])
    test = api_result.numpy_array[:, 0]
    if mode == "paired" or mode == "p":
        assert np.all(test == 0)
    elif mode == "unpaired" or mode == "u":
        assert np.all(test == 1)
