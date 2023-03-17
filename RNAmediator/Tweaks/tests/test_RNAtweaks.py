import os
import random
from tempfile import TemporaryDirectory
import numpy as np
from RNAmediator.Tweaks.RNAtweaks import (
    api_rnaplfold,
    cmd_rnafold,
    cmd_rnaplfold,
    api_rnafold,
    cmd_rnafold,
    _pl_to_array,
    _read_precalc_plfold,
    printdiff,
    _npprint,
    _calc_gibbs,
    _mutate,
)
import pytest
import RNA

TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
TESTDATAPATH = os.path.join(TESTFOLDER, "testdata")
PARPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_DIR = os.path.dirname(os.path.dirname(PARPATH))

def random_sequence(seed: int = 1):
    random.seed(seed)
    return "".join(random.choices(["A", "C", "U", "G"], k=350))


@pytest.fixture
def sequence_string():
    return "AAATTTTTGGGGUUUUUUUCCCCTTTTTttt"


@pytest.mark.parametrize(
    "seq,window,span,constraint",
    [
        ("A" * 500, 70, 70, []),
        ("A" * 500, 70, 70, [("p", 20, 30)]),
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


@pytest.mark.xfail
def test_api_and_cmd_fold(seq, span, constraint, temp=37, coordinates=None):
    cmd_result = cmd_rnafold(seq, span, temp, constraint=constraint)
    api_result = api_rnafold(seq, span, temp, constraint=constraint)
    assert api_result == cmd_result


@pytest.mark.parametrize(
    "seq,window,span,constraint,regions",
    [
        ("A" * 500, 70, 70, [], (3, 5)),
        ("A" * 500, 70, 70, [("p", 20, 30)], (11, 5)),
        (random_sequence(), 70, 70, [("p", 20, 30), ("u", 70, 90)], (3, 5)),
        (random_sequence(), 70, 70, [("p", 20, 30), ("u", 70, 90)], (11, 5)),
        (random_sequence(), 70, 70, [("p", 20, 30), ("u", 70, 90)], (20, 7)),
        (random_sequence(), 70, 70, [("u", 20, 30)], (3, 5)),
        (random_sequence(), 70, 70, [], (3, 5)),
    ],
)
def test_region_param(seq, window, span, constraint, regions):
    api_result_r1 = api_rnaplfold(seq, window, span, constraint=constraint, region=regions[0])
    api_result_r2 = api_rnaplfold(seq, window, span, constraint=constraint, region=regions[1])
    idx = min(regions)
    api_r1_array = api_result_r1.numpy_array[:, 0:idx]
    api_r2_array = api_result_r2.numpy_array[:, 0:idx]
    assert np.allclose(api_r1_array, api_r2_array, equal_nan=True)


@pytest.mark.parametrize(
    "seq,window,span,constraint",
    [
        (random_sequence(5), 70, 70, [("nonesense", 3, 5)]),
        (random_sequence(5), 70, 70, [("p", -3, 5)]),
        (random_sequence(5), 70, 70, [("p", 350, 360)]),
        (random_sequence(5), 70, 70, [("p", 360, 100)]),
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
    "seq,window,span,constraint,temperature",
    [
        (random_sequence(4), 50, 50, ("p", 200, 207), 37),
        (random_sequence(4), 50, 50, ("u", 200, 207), 40),
    ],
)
def test_constraint(seq, window, span, constraint, temperature):
    mode = constraint[0]
    api_result = cmd_rnaplfold(seq, window, span, region=7, temperature=temperature, constraint=[constraint])
    api_result.localize(constraint[1], constraint[2])
    test = api_result.numpy_array[:, 0]
    if mode == "paired" or mode == "p":
        assert np.all(test == 0)
    elif mode == "unpaired" or mode == "u":
        assert np.all(test == 1)


@pytest.fixture
def pl_fold_file():
    """Returns Path to gziped RNAmediator output (same content as pl_npy_file)"""
    return os.path.join(TESTDATAPATH, "pl_test_output.txt.gz")


@pytest.fixture()
def pl_npy_file():
    """Returns Path to gziped RNAmediator output (same content as pl_npy_file)"""

    return os.path.join(TESTDATAPATH, "pl_test_output.npy")


@pytest.fixture
def pl_fold_test_sequence():
    """return sequence used for pl_fold_file and pl_npy_file"""
    seq = "TTTTTTCTTTATAATTATTCCCCTATTTGAAAAATCAACTTGTATATGAGGCAGCAAACACCTTGCAGAGCAGCATTCCCTTTTAGTTTCAGGACGTGGTGGTGGATGGAACCACTGTAACCTGGCCTCCCTCCATGAGAGGAGGGAATCCAGGTGGCCATGTTGAAATGTGCCTGTGTGCAGCAAGGCTTCTGAAATGACAAGAGAGCCCAGCAGCTTCCAAAGCAGCTGTGACTCTGGATCTCACCCATCATCTCTGCTTCTCACTGTTAGAGGAGTGAATCTGTGCTGCCTTAGGAGGCATGGAACCTGGGACTTTTCTTCCTTGTTTAATGTTTAATTTTATTAAAATAATTTGTAAGTGATAGATGTTGATCTCGTGACAAAAGAGAGATTCCCTCTTTATAAAACTATTCTAACTAAAGATCTTTTGTAAGCCCATGTGTTAGAAATAAAACTTGAATATCCCC"
    return seq[0:135]


@pytest.fixture()
def diff_numpy_array():
    return np.load(os.path.join(TESTDATAPATH, "diffnp.npy"))


@pytest.mark.parametrize(
    "ulim,works",
    [
        (7, True),
        (5, True),
        (100, False),
    ],
)
def test_pl_to_array_from_npy(pl_npy_file, ulim, works, caplog):
    try:
        pl_array = _pl_to_array(pl_npy_file, ulim, fmt="npy")
        if not works:
            assert pl_array is None
        else:
            assert isinstance(pl_array, np.ndarray)
    except IndexError:
        pass


@pytest.mark.parametrize(
    "ulim,works",
    [
        (7, True),
        (5, True),
        (100, False),
    ],
)
def test_pl_to_array_from_txt(pl_fold_file, ulim, works, caplog):
    try:
        pl_array = _pl_to_array(pl_fold_file, ulim, fmt="txt")
        if not works:
            assert pl_array is None
        else:
            assert isinstance(pl_array, np.ndarray)
    except Exception as e:
        assert isinstance(e, ValueError) or isinstance(e, IndexError)


def test_pl_to_array_difference(pl_npy_file, pl_fold_file):
    npy_array = _pl_to_array(pl_npy_file, 7, fmt="npy")
    txt_array = _pl_to_array(pl_fold_file, 7, fmt="txt")
    assert np.array_equal(npy_array, txt_array, equal_nan=True)


def test_read_precalc_plfold(pl_fold_file, pl_fold_test_sequence):
    data = _read_precalc_plfold([], pl_fold_file, pl_fold_test_sequence)
    assert data is not None
    assert type(data) == list
    for line in data:
        assert len(line) != 0
        assert type(line) == list


@pytest.mark.parametrize(
    "sequence, start, end, value",
    [("sequence_string", 5, 6, "A"), ("A" * 500, 100, 107, "UUUUUUU")],
)
def test_mutate(sequence, start, end, value, request, caplog):
    if not type(sequence) is str:
        sequence = request.getfixturevalue(sequence)
    assert isinstance(sequence, str)
    mutseq = _mutate(sequence, start, end, value)
    print(f"{mutseq} {value}")
    assert isinstance(mutseq, str)
    assert len(sequence) == len(mutseq)
    assert mutseq[start:end] == value


@pytest.mark.parametrize(
    "outfile",
    [
        "saved_array.npy",
    ],
)
def test_printdiff(diff_numpy_array, outfile, tmp_path):
    outfile = os.path.join(tmp_path, outfile)
    printdiff(diff_numpy_array, outfile)
    assert os.path.exists(outfile)


def test_printdiff_with_defaults(diff_numpy_array, caplog):
    try:
        printdiff(diff_numpy_array)
        for log in caplog.records:
            assert log.levelname == "ERROR"
    except TypeError:
        pass


@pytest.mark.parametrize(
    "outfile",
    [
        "saved_array_print.npy",
        None,
    ],
)
def test_npprint(diff_numpy_array, outfile, capsys, tmp_path):
    if outfile is not None:
        outfile = os.path.join(tmp_path, outfile)
        file = open(outfile, "wb")
        _npprint(diff_numpy_array, file)
        captured = capsys.readouterr()
        assert captured.out == ""
        assert os.path.exists(outfile)
    else:
        _npprint(diff_numpy_array, outfile)
        captured = capsys.readouterr()
        assert captured.out != ""


@pytest.fixture()
def fold_compound():
    fc = RNA.fold_compound(random_sequence(), RNA.md())
    return fc


@pytest.mark.xfail()
def test_calc_gibbs(fold_compound):
    result = _calc_gibbs(fold_compound)
    assert type(result) == float
    raise NotImplementedError


@pytest.mark.xfail
def test_calc_bpp():
    raise NotImplementedError


@pytest.mark.xfail()
def test_calc_nrg():
    raise NotImplementedError


@pytest.mark.xfail()
def test_calc_ddgs():
    raise NotImplementedError


@pytest.mark.xfail()
def test_get_ddg():
    raise NotImplementedError


@pytest.mark.xfail()
def test_get_bppm():
    raise NotImplementedError
