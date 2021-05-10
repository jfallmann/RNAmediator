import os
from argparse import Namespace
import numpy as np
import pytest
import gzip
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
                     loglevel="INFO",
                     logdir="LOGS"
                     )

    return args


def compare_output_folders(test_path: str, expected_path: str):
    assert os.path.exists(expected_path), "expected result does not exist"
    assert os.path.exists(test_path), "test result was not produced"
    results = os.listdir(test_path)
    for subdir in results:
        test_subdir = os.path.join(test_path, subdir)
        expected_subdir = os.path.join(expected_path, subdir)
        assert os.path.exists(test_subdir), f"test: {test_subdir} does not exist"
        assert os.path.exists(expected_subdir), f"expected: {expected_subdir} does not exist"
        expected_files = os.listdir(expected_subdir)
        for file in expected_files:
            test_file = os.path.join(test_subdir, file)
            expected_file = os.path.join(expected_subdir, file)
            assert os.path.exists(test_file), f"test: {test_file} does not exist"
            assert os.path.exists(expected_file), f"expected: {expected_file} does not exist"
            if ".gz" in test_file:
                with gzip.open(test_file) as test, gzip.open(expected_file) as expected:
                    for test_line in test:
                        expected_line = expected.readline()
                        assert expected_line == test_line
            elif test_file.endswith(".npy"):
                test_file = np.load(test_file)
                expected_file = np.load(expected_file)
                assert np.array_equal(test_file, expected_file, equal_nan=True), f"{test_file} does not match expected result"


@pytest.fixture()
def single_constraint_args(default_args):
    default_args.constrain = "test_single_constraint.bed"
    default_args.window = 100
    default_args.procs = 1  # TODO: change number  of procs upon deployment
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = "single_constraint_test"
    default_args.save = 1
    default_args.logdir = "LOG_SINGLE"
    return default_args


@pytest.fixture()
def multi_constraint_args(default_args):
    default_args.sequence = "test.fa.gz"
    default_args.constrain = "test_constraints.bed"
    default_args.window = 100
    default_args.procs = 1  # TODO: change number  of procs upon deployment
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = "multi_constraint_test"
    default_args.save = 10
    default_args.logdir = "LOG_MULTI"
    return default_args


@pytest.fixture()
def sliding_args(default_args):
    default_args.window = 100
    default_args.procs = 4  # TODO: change number  of procs upon deployment
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
    assert os.path.isfile("test_single_constraint.bed")
    assert os.path.isfile("test_constraints.bed")


def test_single_constraint(single_constraint_args):
    pl_main(single_constraint_args)
    expected_path = "single_constraint_result"
    test_path = single_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)


def test_sliding_window(sliding_args):
    pl_main(sliding_args)  # TODO: maybe change output directory to tmpdir
    expected_path = "sliding_result"
    test_path = sliding_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)


def test_multi_constraint(multi_constraint_args):
    pl_main(multi_constraint_args)
    expected_path = "multi_constraint_result"
    test_path = multi_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)






