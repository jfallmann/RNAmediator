from __future__ import annotations

import gzip
import os
import subprocess
import sys
from argparse import Namespace
from tempfile import TemporaryDirectory, NamedTemporaryFile
import numpy as np
import pytest

TESTPATH = os.path.dirname(os.path.abspath(__file__))
PARPATH = os.path.dirname(TESTPATH)
sys.path.append(PARPATH)
from RIssmed.ConstraintPLFold import main as pl_main
from RIssmed.ConstraintPLFold import fold_unconstraint
from Bio import SeqIO

EXPECTED_LOGS = os.path.join(TESTPATH, "Expected_Logs")
EXPECTED_RESULTS = os.path.join(TESTPATH, "Expected_Results")

tmp_dir = TemporaryDirectory()
TMP_TEST_DIR = tmp_dir.name


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
    test_lines = set()
    with open(test_log) as test_file:
        for line in test_file:
            test_line = " ".join(line.split(" ")[4:])
            test_lines.add(test_line)
    with open(expected_log) as expected_file:
        for line in expected_file:
            expected_line = " ".join(line.split(" ")[4:])
        if "CLI:" and "JetBrains" and "Running ConstraintPLFold on" and "DONE: output in" not in expected_line:
            assert expected_line in test_lines


@pytest.fixture()
def single_constraint_args(default_args):
    default_args.constrain = os.path.join(TESTPATH, "test_single_constraint.bed")
    default_args.window = 100
    default_args.procs = 1
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "single_constraint_test")
    default_args.save = 1
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_SINGLE")
    return default_args


@pytest.fixture()
def multi_constraint_args(default_args):
    default_args.sequence = os.path.join(TESTPATH, "test.fa.gz")
    default_args.constrain = os.path.join(TESTPATH, "test_constraints.bed")
    default_args.window = 100
    default_args.procs = os.cpu_count() - 1 or 1
    default_args.conslength = 7
    default_args.region = 7
    default_args.unpaired = "unpaired"
    default_args.paired = "paired"
    default_args.unconstrained = "raw"
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
    default_args.unconstrained = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "sliding_test")
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_SLIDING")

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
    pl_main(sliding_args)
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


@pytest.mark.parametrize(
    "seq_id,region,window,span,unconstraint,save,outdir,seq",
    [("onlyA", 7, 100, 60, "raw", 1, "onlyA", "A" * 500),
     ("testseq2", 7, 100, 60, "raw", 1, "testseq2", os.path.join(TESTPATH, "test_single.fa"))]
)
def test_fold_unconstraint(seq_id, region, window, span, unconstraint, save, outdir, seq):
    if os.path.isfile(seq):
        seq = str(SeqIO.read(seq, format="fasta").seq)
        p = 0
    outdir = os.path.join(TMP_TEST_DIR, outdir)
    # get the resulting np. array via the command line of RNAplfold
    cmd_result = run_pl_fold(seq, 100, 60, u=region)
    cmd_array = cmd_result.get_numpy_array()

    # runs RIssmed to produce output files using the same input as the command line
    fold_unconstraint(seq, seq_id, region, window, span, unconstraint, save, outdir)
    test_file_path = os.path.join(outdir, seq_id)
    test_files = os.listdir(test_file_path)
    for test_file in test_files:
        test_file = os.path.join(test_file_path, test_file)
        if ".gz" in test_file:
            test_result = PLFoldOutput.from_file(test_file)
            test_array = test_result.get_numpy_array()  # TODO: Duplication can get replaced via comparison of outputs
            array_difference = test_array - cmd_array
            max_difference = np.max(np.abs(array_difference), where=~np.isnan(array_difference), initial=-1)
            max_diff_idx = np.unravel_index(np.nanargmax(array_difference), array_difference.shape)
            assert np.array_equal(test_array, cmd_array, equal_nan=True), \
                f"detected high difference between RIssmed and command line result " \
                f"with a max of {max_difference} at index: {max_diff_idx}"
        else:
            test_result = PLFoldOutput.from_rissmed_numpy_output(test_file)
            test_array = test_result.get_numpy_array()
            array_difference = test_array - cmd_array
            max_difference = np.max(np.abs(array_difference), where=~np.isnan(array_difference), initial=-1)
            max_diff_idx = np.unravel_index(np.nanargmax(array_difference), array_difference.shape)
            # TODO: The numpy version of RIssmed looks very strange. Need to talk to Jörg
            # TODO: differences of 0.012628 ... for index 299,0 of second test array seems to be very high
            assert np.array_equal(test_array, cmd_array, equal_nan=True), \
                f"detected high difference between RIssmed and command line result " \
                f"with a max of {max_difference} at index: {max_diff_idx}"


def run_pl_fold(sequence, window, span, start=None, end=None, mode="unpaired", u=30):
    with TemporaryDirectory() as tmp_dir, NamedTemporaryFile(mode="r+") as constraint_file:
        if mode == "paired":
            const = "F"
        elif mode == "unpaired":
            const = "P"
        else:
            raise ValueError("mode must be either paired or unpaired")
        if start is not None and end is not None:
            constraint_string = f"{const} {start} {0} {end - start}"
        else:
            constraint_string = ""
        constraint_file.write(constraint_string)
        constraint_file.seek(0)
        rnaplfold = subprocess.Popen(["RNAplfold", "-W", str(window), "-L", str(span),
                                      "--commands", constraint_file.name, "--auto-id", "-u", str(u)],
                                     stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                     cwd=tmp_dir)
        stdout, stderr = rnaplfold.communicate(sequence.encode("utf-8"))
        assert stderr == b"", f"call to RNApfold went wrong: \n {stderr.decode()}"
        file = os.path.join(tmp_dir, "sequence_0001_lunp")
        rnaplfold_output = PLFoldOutput.from_file(file)
        return rnaplfold_output


class PLFoldOutput:
    def __init__(self, text: str):
        self.text = self._sanitize(text)

        self._array = None

    def __str__(self):
        return self.text

    def __eq__(self, other: PLFoldOutput):
        if type(other) != PLFoldOutput:
            return False
        return np.array_equal(self.get_numpy_array(), other.get_numpy_array(), equal_nan=True)


    @staticmethod
    def _sanitize(text: str):
        list_text = text.split("\n")
        if list_text[0].startswith("1\t"):
            length = len(list_text[0].split("\t")) - 1
            part = [str(x + 1) for x in range(length)]
            list_text = ["#unpaired probabilities", " #$i\tl=" + "\t".join(part)] + list_text
            text = "\n".join(list_text)
        elif list_text[0].startswith("#") and list_text[1].startswith(" #"):
            pass
        else:
            raise ValueError("text seems not to be a valid plfold string")
        text = text.replace("nan", "NA")
        return text

    def get_numpy_array(self):
        if self._array is None:
            array = []
            for line in self.text.split("\n")[:-1]:
                if not line.startswith("#") and not line.startswith(" #"):
                    data = line.split("\t")[1:]
                    data = [float(x) if x != "NA" else np.nan for x in data]
                    array.append(data)
            array = np.array(array)
            self._array = array
        return self._array

    def set_array(self, array: np.ndarray):
        self._array = array

    @classmethod
    def from_file(cls, file_path: str):
        assert os.path.isfile(file_path), f"{file_path} is not a valid file path"
        if ".gz" in file_path:
            with gzip.open(file_path, "rt") as handle:
                file_data = handle.read()
        else:
            with open(file_path, "r") as handle:
                file_data = handle.read()
        return PLFoldOutput(file_data)

    @classmethod
    def from_rissmed_numpy_output(cls, file_path: str):
        assert os.path.isfile(file_path), f"{file_path} is not a valid file path"
        ris_array = np.load(file_path)
        array = np.squeeze(ris_array)
        array_string = "\n".join(
            ['\t'.join([str(x + 1)] + ['%.5f' % num for num in array[x]]) for x in range(len(array))])
        output = PLFoldOutput(array_string)
        output.set_array(array)
        return output
