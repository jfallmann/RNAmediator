import os
from argparse import Namespace
import numpy as np
import pytest

from RIssmed.ConstraintPLFold import main as pl_main


@pytest.fixture()
def default_args():
    args = Namespace(sequence="test_single.fa",
                     window=240,
                     span=60,
                     region=1,
                     multi=2,
                     cutoff=0.01,
                     unconstraint="STDOUT",
                     unpaired="STDOUT",
                     paired="STDOUT",
                     length=100,
                     gc=0,
                     number=1,
                     constrain="sliding",
                     conslength=1,
                     temprange="",
                     alphabet="AUCG",
                     save=0,
                     outdir="",
                     procs=1,
                     vrna="",
                     pattern="",
                     genes="",
                     verbosity=0,
                     loglevel="WARNING",
                     logdir="LOGS"
                     )

    return args


@pytest.fixture()
def sliding_args(default_args):
    default_args.window = 100
    default_args.procs = 7  # TODO: change number  of procs upon deployment
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = "sliding_test"
    return default_args


def test_data_available():
    assert os.path.exists("test.fa.gz"), "Test data not available"
    assert os.path.exists("test_single.fa"), "Test data not available"


def test_sliding_window(sliding_args):
    pl_main(sliding_args)  # TODO: maybe change output directory to tmpdir
    expected_path = "sliding_result/ENSG00000270742"
    assert os.path.exists(expected_path), "expected result does not exist"
    test_path = "sliding_test/ENSG00000270742"
    tests = os.listdir(test_path)
    for element in tests:
        test_el = os.path.join(test_path, element)
        expected_el = os.path.join(expected_path, element)
        assert os.path.exists(test_el), f"test: {test_el} does not exist"
        assert os.path.exists(expected_el), f"expected: {expected_el} does not exist"
        test_el = np.load(test_el)
        expected_el = np.load(expected_el)
        assert np.array_equal(test_el, expected_el, equal_nan=True), f"{element} does not match expected result"


@pytest.fixture()
def single_constraint_args(default_args):
    default_args.window = 100


def test_single_constraint():
    pass


