from __future__ import annotations

import gzip
import os
import subprocess
import sys
from argparse import Namespace
from tempfile import TemporaryDirectory, NamedTemporaryFile

import numpy as np
import pytest

TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
PARPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PARPATH)
from RIssmed.ConstraintPLFold import main as pl_main
from RIssmed.ConstraintPLFold import fold_unconstraint, constrain_seq
from Bio import SeqIO
from RIssmed.lib.Collection import expand_window, localize_window

EXPECTED_LOGS = os.path.join(TESTFOLDER, "Expected_Logs")
EXPECTED_RESULTS = os.path.join(TESTFOLDER, "Expected_Results")
TESTDATAPATH = os.path.join(TESTFOLDER, "testdata")

TMP_DIR = TemporaryDirectory()
TMP_TEST_DIR = TMP_DIR.name


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
                with gzip.open(test_file) as test, gzip.open(expected_file) as expected:
                    for test_line in test:
                        expected_line = expected.readline()
                        assert expected_line == test_line, "lines are different"
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
    # Needs to be converted as api output might include nans in contrast to the zeros of cmd line
    test_array = np.nan_to_num(test_array)
    expected_array = np.nan_to_num(expected_array)
    array_difference = test_array - expected_array
    max_difference = np.max(np.abs(array_difference), where=~np.isnan(array_difference), initial=-1)
    max_diff_idx = np.unravel_index(np.nanargmax(array_difference), array_difference.shape)
    assert np.allclose(test_array, expected_array, equal_nan=True, atol=0.0000001), \
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
    default_args.unconstrained = "raw"
    default_args.outdir = os.path.join(TMP_TEST_DIR, "single_constraint_test")
    default_args.save = 1
    default_args.logdir = os.path.join(TMP_TEST_DIR, "LOG_SINGLE")
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
    compare_logs(test_log=test_log, expected_log=expected_log)


def test_sliding_window(sliding_args):
    pl_main(sliding_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "sliding_result")
    test_path = sliding_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(sliding_args.logdir, os.listdir(sliding_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Sliding_Constraint.log")
    compare_logs(test_log=test_log, expected_log=expected_log)


def test_multi_constraint(multi_constraint_args):
    pl_main(multi_constraint_args)
    expected_path = os.path.join(EXPECTED_RESULTS, "multi_constraint_result")
    test_path = multi_constraint_args.outdir
    compare_output_folders(test_path=test_path, expected_path=expected_path)
    test_log = os.path.join(multi_constraint_args.logdir, os.listdir(multi_constraint_args.logdir)[0])
    expected_log = os.path.join(EXPECTED_LOGS, "Multi_Constraint.log")
    compare_logs(test_log=test_log, expected_log=expected_log)


@pytest.mark.parametrize(
    "seq_id,region,window,span,unconstraint,save,outdir,seq",
    [("onlyA", 7, 100, 60, "raw", 1, "onlyA", "A" * 500),
     ("testseq2", 7, 100, 60, "raw", 1, "testseq2", os.path.join(TESTDATAPATH, "test_single.fa"))]
)
def test_fold_unconstraint(seq_id, region, window, span, unconstraint, save, outdir, seq):
    if os.path.isfile(seq):
        seq = str(SeqIO.read(seq, format="fasta").seq)
    seq = seq.upper().replace("T", "U")
    outdir = os.path.join(TMP_TEST_DIR, outdir)
    # get the resulting np. array via the command line of RNAplfold
    cmd_result = run_pl_fold(seq, window, span, u=region)
    cmd_array = cmd_result.get_numpy_array()

    # runs RIssmed to produce output files using the same input as the command line
    fold_unconstraint(seq, seq_id, region, window, span, unconstraint, save, outdir)
    test_file_path = os.path.join(outdir, seq_id)
    test_files = os.listdir(test_file_path)
    assert len(test_files) != 0, "fold_unconstraint  went wrong"
    for test_file in test_files:
        test_file = os.path.join(test_file_path, test_file)
        if ".gz" in test_file:
            test_result = PLFoldOutput.from_file(test_file)
            test_array = test_result.get_numpy_array()
            compare_arrays(test_array, cmd_array)
        else:
            test_result = PLFoldOutput.from_rissmed_numpy_output(test_file)
            test_array = test_result.get_numpy_array()
            compare_arrays(test_array, cmd_array)


@pytest.mark.parametrize(
    "seq_id,start,end,window,span,region,multi,paired,unpaired,save,outdir,pl_data,unconstraint,seq",
    [("onlyA", 200, 207, 100, 60, 7, 1, "paired", "unpaired", 1, "onlyA", {'up': []}, "raw", "A" * 500),
     ("testseq2", 200, 207, 100, 60, 7, 1, "paired", "unpaired", 1, "testseq2", {'up': []}, "raw",
      os.path.join(TESTDATAPATH, "test_single.fa")),
     ("testseq3", 200, 207, 100, 60, 7, 2, "paired", "unpaired", 1, "testseq3", {'up': []}, "raw",
      os.path.join(TESTDATAPATH, "test_single.fa"))

     ]
)
def test_fold_constraint(seq_id, start, end, window, span, region, multi, paired, unpaired, save, outdir, pl_data,
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

    cmd_unpaired = run_pl_fold(seqtofold, window, span, locstart+1, locend+2, u=region)
    cmd_unpaired.localize(locws, locwe+1)
    cmd_unpaired_array = cmd_unpaired.get_numpy_array()
    cmd_paired = run_pl_fold(seqtofold, window, span, locstart+1, locend+2, mode="paired", u=region)
    cmd_paired.localize(locws, locwe+1)
    cmd_paired_array = cmd_paired.get_numpy_array()

    constrain_seq(seq_id, seq, start, end, window, span, region, multi, paired, unpaired, save, outdir, pl_data,
                  unconstraint=unconstraint)
    test_file_path = os.path.join(outdir, seq_id)
    test_files = os.listdir(test_file_path)
    assert len(test_files) >= 3, "constrain seq went wrong"
    for test_file in test_files:
        test_file = os.path.join(test_file_path, test_file)
        if ".gz" in test_file:
            if unpaired in test_file:
                test_result = PLFoldOutput.from_file(test_file)
                test_array = test_result.get_numpy_array()
                compare_arrays(test_array, cmd_unpaired_array)
            if paired in test_file and unpaired not in test_file:
                test_result = PLFoldOutput.from_file(test_file)
                test_array = test_result.get_numpy_array()
                # TODO: seems like the API produces nan at paired constraint positions instead of 0.
                compare_arrays(test_array, cmd_paired_array)


def run_pl_fold(sequence, window, span, start=None, end=None, mode="unpaired", u=30) -> PLFoldOutput:
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
            for line in self.text.split("\n"):
                if not line.startswith("#") and not line.startswith(" #") and not line == "":
                    data = line.split("\t")[1:]
                    data = [float(x) if x != "NA" else np.nan for x in data]
                    array.append(data)
            array = np.array(array)
            self._array = array
        return self._array

    def set_array(self, array: np.ndarray):
        self._array = array

    def localize(self, start: int, end: int):
        if self._array is not None:
            self._array = self._array[start:end]
        text_list = self.text.split("\n")[start+2:end+2]
        for x, item in enumerate(text_list):
            text_list[x] = "\t".join([str(x+1)] + item.split("\t")[1:])
        self.text = self._sanitize("\n".join(text_list))

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
            ['\t'.join([str(x + 1)] + ['%.8f' % num for num in array[x]]) for x in range(len(array))])
        output = PLFoldOutput(array_string)
        output.set_array(array)
        return output
