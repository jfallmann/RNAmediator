from __future__ import annotations
import gzip
import os
from argparse import Namespace
from tempfile import TemporaryDirectory
import numpy as np
import pytest
from RIssmed.ConstraintPLFold import main as pl_main
from RIssmed.ConstraintPLFold import fold_unconstraint, constrain_seq
from Bio import SeqIO
from RIssmed.RNAtweaks.RNAtweaks import PLFoldOutput, cmd_rnaplfold
from RIssmed.RNAtweaks.RIssmed import expand_window, localize_window


TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
PARPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPECTED_LOGS = os.path.join(TESTFOLDER, "Expected_Logs")
EXPECTED_RESULTS = os.path.join(TESTFOLDER, "Expected_Results")
TESTDATAPATH = os.path.join(TESTFOLDER, "testdata")

TMP_DIR = TemporaryDirectory()
TMP_TEST_DIR = "./"


@pytest.fixture()
def default_args():
    args = Namespace(sequence=os.path.join(TESTDATAPATH, "test_single.fa"),
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
    results = os.listdir(expected_path)
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
                with gzip.open(test_file, "rt") as test, gzip.open(expected_file, "rt") as expected:
                    for x, test_line in enumerate(test):
                        expected_line = expected.readline().rstrip().split("\t")
                        test_line = test_line.rstrip().split("\t")
                        for y in range(len(expected_line)):
                            expected_num = expected_line[y]
                            test_num = test_line[y]
                            if expected_num != "nan" and test_num != "nan":
                                expected_num = float(expected_num)
                                test_num = float(test_num)
                                assert np.allclose(expected_num, test_num, atol=0.0000001), f"files are different at {x}, {y}"
            elif test_file.endswith(".npy"):
                test_file = np.load(test_file)
                expected_file = np.load(expected_file)
                compare_arrays(test_file, expected_file)


def compare_logs(test_log: str, expected_log: str):
    test_lines = set()
    substrings = ["CLI:", "JetBrains", "Running ConstraintPLFold on", "/home/rabsch/Documents/RIssmed/", "/tmp/"]
    with open(test_log) as test_file:
        for line in test_file:
            test_line = " ".join(line.split(" ")[4:])
            test_lines.add(test_line)
    with open(expected_log) as expected_file:
        for line in expected_file:
            expected_line = " ".join(line.split(" ")[4:])
            if not any(substring in expected_line for substring in substrings):
                assert expected_line in test_lines


def compare_arrays(test_array: np.ndarray, expected_array: np.ndarray):
    array_difference = test_array - expected_array
    max_difference = np.max(np.abs(array_difference), where=~np.isnan(array_difference), initial=-1)
    max_diff_idx = np.unravel_index(np.nanargmax(array_difference), array_difference.shape)
    assert np.allclose(test_array, expected_array, equal_nan=True, atol=0.000001), \
        f"detected high difference between RIssmed and command line result " \
        f"with a max of {max_difference} at index: {max_diff_idx}"


@pytest.fixture()
def single_constraint_args(default_args):
    default_args.constrain = os.path.join(TESTDATAPATH, "test_single_constraint.bed")
    default_args.window = 100
    default_args.procs = 1
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstraint = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "single_constraint_test")
    default_args.save = 1
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_SINGLE")
    return default_args


@pytest.fixture()
def paired_constraint_args(default_args):
    default_args.sequence = os.path.join(TESTDATAPATH, "test.fa.gz")
    default_args.constrain = os.path.join(TESTDATAPATH, "paired_constraints.bed")
    default_args.conslength = 7
    default_args.window = 60
    default_args.span = 60
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstraint = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "paired_constraint_test")
    default_args.save = 1
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_PAIRED")
    return default_args


@pytest.fixture()
def multi_constraint_args(default_args):
    default_args.sequence = os.path.join(TESTDATAPATH, "test.fa.gz")
    default_args.constrain = os.path.join(TESTDATAPATH, "test_constraints.bed")
    default_args.window = 100
    default_args.procs = os.cpu_count() - 1 or 1
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstraint = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "multi_constraint_test")
    default_args.save = 1
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_MULTI")
    return default_args


@pytest.fixture()
def sliding_args(default_args):
    default_args.window = 100
    default_args.procs = os.cpu_count() - 1 or 1
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstraint = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "sliding_test")
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_SLIDING")

    return default_args


@pytest.fixture()
def parafold_args(default_args):
    default_args.sequence = os.path.join(TESTDATAPATH, "parafold_test.fa")
    constrain = os.path.join(TESTDATAPATH, "parafold_test.bed")
    default_args.window = 100
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstraint = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "parafold_test")
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_PARAFOLD")
    default_args.constrain = f"ono,{constrain}"
    default_args.save = 1
    return default_args


def test_data_available():
    assert os.path.exists(os.path.join(TESTDATAPATH, "test.fa.gz")), "Test data not available"
    assert os.path.exists(os.path.join(TESTDATAPATH, "test_single.fa")), "Test data not available"
    assert os.path.isfile(os.path.join(TESTDATAPATH, "test_single_constraint.bed"))
    assert os.path.isfile(os.path.join(TESTDATAPATH, "test_constraints.bed"))
    assert os.path.isfile(os.path.join(TESTDATAPATH, "parafold_test.bed"))
    assert os.path.isfile(os.path.join(TESTDATAPATH, "parafold_test.fa"))


def test_paired_constraint(paired_constraint_args):
    pl_main(paired_constraint_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "paired_constraint_result")
    test_path = paired_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    # TODO: Not done yet but settings for paired constraints work now


def test_parafold(parafold_args):
    pl_main(parafold_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "parafold_result")
    test_path = parafold_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    # Here only files are compared as the parafold function seems to be a bit buggy


def test_single_constraint(single_constraint_args):
    pl_main(single_constraint_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "single_constraint_result")
    test_path = single_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(single_constraint_args.logdir, os.listdir(single_constraint_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Single_Constraint.log")


def test_sliding_window(sliding_args):
    pl_main(sliding_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "sliding_result")
    test_path = sliding_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(sliding_args.logdir, os.listdir(sliding_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Sliding_Constraint.log")


def test_multi_constraint(multi_constraint_args):
    pl_main(multi_constraint_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "multi_constraint_result")
    test_path = multi_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(multi_constraint_args.logdir, os.listdir(multi_constraint_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Multi_Constraint.log")


@pytest.mark.parametrize(
    "seq_id,region,window,span,unconstraint,save,outdir,seq",
    [("onlyA", 7, 100, 60, "raw", 1, "onlyA", "A" * 500),
     ("testseq2", 7, 100, 60, "raw", 2, "testseq2", os.path.join(TESTDATAPATH, "test_single.fa"))]
)
def test_fold_unconstraint(seq_id, region, window, span, unconstraint, save, outdir, seq):
    if os.path.isfile(seq):
        seq = str(SeqIO.read(seq, format="fasta").seq)
    seq = seq.upper().replace("T", "U")
    outdir = os.path.join(TMP_TEST_DIR, outdir)
    # get the resulting np. array via the command line of RNAplfold
    cmd_result = cmd_rnaplfold(seq, window, span, region=region)
    cmd_array = cmd_result.numpy_array

    # runs RIssmed to produce output files using the same input as the command line
    fold_unconstraint(seq, seq_id, region, window, span, unconstraint, save, outdir)
    test_file_path = os.path.join(outdir, seq_id)
    test_files = os.listdir(test_file_path)
    assert len(test_files) != 0, "fold_unconstraint  went wrong"
    for test_file in test_files:
        test_file = os.path.join(test_file_path, test_file)
        if ".gz" in test_file:
            test_result = PLFoldOutput.from_file(test_file)
            test_array = test_result.numpy_array
            compare_arrays(test_array, cmd_array)
        else:
            test_result = PLFoldOutput.from_rissmed_numpy_output(test_file)
            test_array = test_result.numpy_array
            compare_arrays(test_array, cmd_array)


@pytest.mark.parametrize(
    "seq_id,start,end,window,span,region,multi,paired,unpaired,save,outdir,unconstraint,seq",
    [("onlyA", 200, 207, 100, 60, 7, 1, "paired", "unpaired", 1, "onlyA", "raw", "A" * 500),
     ("testseq2", 200, 207, 100, 60, 7, 1, "paired", "unpaired", 1, "testseq2", "raw",
      os.path.join(TESTDATAPATH, "test_single.fa")),
     ("testseq3", 200, 207, 100, 60, 7, 2, "paired", "unpaired", 1, "testseq3", "raw",
      os.path.join(TESTDATAPATH, "test_single.fa"))

     ]
)
def test_fold_constraint(seq_id, start, end, window, span, region, multi, paired, unpaired, save, outdir,
                         unconstraint, seq):
    if os.path.isfile(seq):
        seq = str(SeqIO.read(seq, format="fasta").seq)
    seq = seq.upper().replace("T", "U")
    outdir = os.path.join(TMP_TEST_DIR, outdir)
    # TODO: Maybe rewrite this part in the ConstraintPLfold and just import the function
    tostart, toend = expand_window(start, end, window, multi, len(seq))
    seqtofold = str(seq[tostart - 1:toend])
    locws, locwe = localize_window(start, end, window, len(seq))
    locws = locws - tostart
    locwe = locwe - tostart
    locstart = start - tostart
    locend = end - tostart

    cmd_unpaired = cmd_rnaplfold(seqtofold, window, span, region=region, constraint=[("unpaired", locstart, locend + 1)])
    cmd_unpaired.localize(locws, locwe+1)
    cmd_unpaired_array = cmd_unpaired.numpy_array
    cmd_paired = cmd_rnaplfold(seqtofold, window, span, region=region, constraint=[("paired", locstart, locend + 1)])
    cmd_paired.localize(locws, locwe+1)
    cmd_paired_array = cmd_paired.numpy_array

    constrain_seq(seq_id, seq, start, end, window, span, region, multi, paired, unpaired, save, outdir,
                  unconstraint=unconstraint)
    test_file_path = os.path.join(outdir, seq_id)
    test_files = os.listdir(test_file_path)
    assert len(test_files) >= 3, "constrain seq went wrong"
    for test_file in test_files:
        test_file = os.path.join(test_file_path, test_file)
        if ".gz" in test_file:
            if unpaired in test_file:
                test_result = PLFoldOutput.from_file(test_file)
                test_array = test_result.numpy_array
                compare_arrays(test_array, cmd_unpaired_array)
            if paired in test_file and unpaired not in test_file:
                test_result = PLFoldOutput.from_file(test_file)
                test_array = test_result.numpy_array
                # TODO: seems like the API produces nan at paired constraint positions instead of 0.
                compare_arrays(test_array, cmd_paired_array)





