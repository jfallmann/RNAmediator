import gzip
import os
import sys
from argparse import Namespace

import numpy as np
import pytest

TESTPATH = os.path.dirname(os.path.abspath(__file__))
PARPATH = os.path.dirname(TESTPATH)
sys.path.append(PARPATH)
from RIssmed.ConstraintPLFold import main as pl_main

EXPECTED_LOGS = os.path.join(TESTPATH, "Expected_Logs")
EXPECTED_RESULTS = os.path.join(TESTPATH, "Expected_Results")


@pytest.fixture()
def default_args():
    args = Namespace(sequence=os.path.join(TESTPATH, "test_single.fa"),
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
                     loglevel="DEBUG",
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
                assert np.array_equal(test_file, expected_file,
                                      equal_nan=True), f"{test_file} does not match expected result"


def compare_logs(test_log: str, expected_log: str):
    expected_lines = set()
    with open(expected_log) as expected_file:
        for line in expected_file:
            expected_line = " ". join(line.split(" ")[4:])
            expected_lines.add(expected_line)
    with open(test_log) as test_file:
        for line in test_file:
            test_line = " ". join(line.split(" ")[4:])
        if "CLI:" and "JetBrains" not in test_line and "Running ConstraintPLFold on" not in test_line:
            assert test_line in expected_lines



@pytest.fixture()
def single_constraint_args(default_args):
    default_args.constrain = os.path.join(TESTPATH, "test_single_constraint.bed")
    default_args.window = 100
    default_args.procs = 1  # TODO: change number  of procs upon deployment
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = os.path.join(TESTPATH, "single_constraint_test")
    default_args.save = 1
    default_args.logdir = os.path.join(TESTPATH, "LOG_SINGLE")
    return default_args


@pytest.fixture()
def multi_constraint_args(default_args):
    default_args.sequence = os.path.join(TESTPATH, "test.fa.gz")
    default_args.constrain = os.path.join(TESTPATH, "test_constraints.bed")
    default_args.window = 100
    default_args.procs = os.cpu_count() - 1 or 1 # TODO: change number  of procs upon deployment
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = os.path.join(TESTPATH, "multi_constraint_test")
    default_args.save = 10
    default_args.logdir = os.path.join(TESTPATH, "LOG_MULTI")
    return default_args


@pytest.fixture()
def sliding_args(default_args):
    default_args.window = 100
    default_args.procs = os.cpu_count() - 1 or 1  # TODO: change number  of procs upon deployment
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = os.path.join(TESTPATH, "sliding_test")
    default_args.logdir = os.path.join(TESTPATH, "LOG_SLIDING")

    return default_args


def test_data_available():
    assert os.path.exists(os.path.join(TESTPATH, "test.fa.gz")), "Test data not available"
    assert os.path.exists(os.path.join(TESTPATH, "test_single.fa")), "Test data not available"
    assert os.path.isfile(os.path.join(TESTPATH, "test_single_constraint.bed"))
    assert os.path.isfile(os.path.join(TESTPATH, "test_constraints.bed"))


def test_single_constraint(single_constraint_args):
    pl_main(single_constraint_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "single_constraint_result")
    test_path = single_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(single_constraint_args.logdir, os.listdir(single_constraint_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Single_Constraint.log")
    compare_logs(test_log=test_log, expected_log=expected_log)
    os.system(f"rm {test_path} -r")
    os.system(f"rm {single_constraint_args.logdir} -r")


def test_sliding_window(sliding_args):
    pl_main(sliding_args)  # TODO: maybe change output directory to tmpdir
    expected_path = os.path.join(EXPECTED_RESULTS, "sliding_result")
    test_path = sliding_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(sliding_args.logdir, os.listdir(sliding_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Sliding_Constraint.log")
    compare_logs(test_log=test_log, expected_log=expected_log)
    os.system(f"rm {test_path} -r")
    os.system(f"rm {sliding_args.logdir} -r")


def test_multi_constraint(multi_constraint_args):
    pl_main(multi_constraint_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "multi_constraint_result")
    test_path = multi_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(multi_constraint_args.logdir, os.listdir(multi_constraint_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Multi_Constraint.log")
    compare_logs(test_log=test_log, expected_log=expected_log)
    os.system(f"rm {test_path} -r")
    os.system(f"rm {multi_constraint_args.logdir} -r")

