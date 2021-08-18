import os
from io import StringIO, TextIOWrapper
from collections import defaultdict
import pytest
from tempfile import TemporaryDirectory
from RIssmed.Tweaks.FileProcessor import parseseq, idfromfa, \
    read_constraints_from_bed, make_outdir, read_paired_constraints_from_bed

TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
TESTDATAPATH = os.path.join(TESTFOLDER, "testdata")
TMP_DIR = TemporaryDirectory()
TMP_TEST_DIR = TMP_DIR.name


@pytest.mark.parametrize(
    "sequence",
    [
        os.path.join(TESTDATAPATH, "test.fa"),
        "AAAAAATTTTTTGGGGG",
        os.path.join(TESTDATAPATH, "test.fa.gz"),
    ]
)
def test_parseseq(sequence):
    seq = parseseq(sequence)
    assert type(seq) in [StringIO, TextIOWrapper]


@pytest.mark.parametrize(
    "fa_id",
    [
        "Seq1-40:random:nochrom:(.)",
        "foo",
        "ENSG00000270742::chr1:61124732-61125202(+)"

    ]
)
def test_idfromfa(fa_id):
    goi, chrom, strand = idfromfa(fa_id)
    if goi == fa_id:
        assert chrom == strand == "na"


@pytest.mark.parametrize("linewise", [True, False])
@pytest.mark.parametrize(
    "bedfile,expected",
    [
        (open(os.path.join(TESTDATAPATH, "test_constraints.bed")), True),
        ("foo", False),
        (StringIO("chr1	8	17	ENSG00000273544	u	-\n"
                  "chr1	9	18	ENSG00000201457	u	-"), True),
        ("chr1	8	17	ENSG00000273544	u	-\n", False)
    ]
)
def test_constraints_from_bed(bedfile, linewise, expected, caplog):
    if isinstance(bedfile, StringIO) or isinstance(bedfile, TextIOWrapper):
        bedfile.seek(0)
    constraints = read_constraints_from_bed(bedfile, linewise)

    assert (constraints is not None) is expected

    if constraints is None:
        assert caplog.text != ""
    else:
        assert type(constraints) == defaultdict
        for entry in constraints:
            for cons in constraints[entry]:
                assert "|" in cons
                assert "-" in cons
        assert caplog.text == ""


@pytest.mark.parametrize("linewise", [True, False])
@pytest.mark.parametrize(
    "bedfile,expected",
    [
        (open(os.path.join(TESTDATAPATH, "paired_constraints.bed")), True),
        ("foo", False),
        (StringIO("chr1	8	17	ENSG00000273544	u	-\t"
                  "chr1	9	18	ENSG00000201457	u	-"), True),
        ("chr1	8	17	ENSG00000273544	u	-\n", False)
    ]
)
def test_paired_constraints_from_bed(bedfile, linewise, expected, caplog):
    if isinstance(bedfile, StringIO) or isinstance(bedfile, TextIOWrapper):
        bedfile.seek(0)
    constraints = read_paired_constraints_from_bed(bedfile, linewise)

    assert (constraints is not None) is expected

    if constraints is None:
        assert caplog.text != ""
        print(caplog.text)
    else:
        assert type(constraints) == defaultdict
        assert len(constraints) > 0
        for entry in constraints:
            assert len(constraints[entry]) > 0
            for cons in constraints[entry]:
                assert "|" in cons
                assert "-" in cons
                assert ":" in cons
        assert caplog.text == ""


@pytest.mark.parametrize(
    "dir_path",
    [
        os.path.join(TMP_TEST_DIR, "Foo"),
        os.path.join(TMP_TEST_DIR, "bla/Foo"),
    ]
)
def test_make_outdir(dir_path, caplog):
    outdir = make_outdir(dir_path)
    assert os.path.exists(outdir)
    assert caplog.text == ""
