import os
from io import StringIO, TextIOWrapper
from collections import defaultdict
import pytest
import gzip
from tempfile import TemporaryDirectory
from RNAmediator.Tweaks.FileProcessor import (
    parseseq,
    idfromfa,
    read_constraints_from_bed,
    make_outdir,
    read_paired_constraints_from_bed,
    parse_annotation_bed,
)

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
        StringIO("ATATATGATATATGAATAGATATATATAGATA"),
    ],
)
def test_parseseq(sequence):
    seq = parseseq(sequence)
    assert type(seq) in [StringIO, TextIOWrapper]


@pytest.mark.parametrize("fa_id", ["Seq1-40:random:nochrom:(.)", "foo", "ENSG00000270742::chr1:61124732-61125202(+)"])
def test_idfromfa(fa_id):
    goi, chrom, strand = idfromfa(fa_id)
    if goi == fa_id:
        assert chrom == strand == "na"


@pytest.mark.parametrize("linewise", [True, False])
@pytest.mark.parametrize(
    "bedfile,expected",
    [
        (open(os.path.join(TESTDATAPATH, "test_constraints.bed")), True),
        (StringIO("chr1\t8\t17\tENSG00000273544\tu\t-\n" "chr1\t9\t18\tENSG00000201457\tu\t-"), True),
    ],
)
def test_constraints_from_bed(bedfile, linewise, expected, caplog):
    if isinstance(bedfile, StringIO) or isinstance(bedfile, TextIOWrapper):
        bedfile.seek(0)
    constraints = read_constraints_from_bed(bedfile, linewise)

    assert (constraints is not None) is expected

    if constraints is None:
        assert caplog.text != ""
    else:
        if not isinstance(constraints, Exception):
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
        ("foo", True),
        (StringIO("chr1	8	17	ENSG00000273544	u	-\t" "chr1	9	18	ENSG00000201457	u	-"), True),
        ("chr1	8	17	ENSG00000273544	u	-\n", True),
    ],
)
def test_paired_constraints_from_bed(bedfile, linewise, expected, caplog):
    if isinstance(bedfile, StringIO) or isinstance(bedfile, TextIOWrapper):
        bedfile.seek(0)

    try:
        constraints = read_paired_constraints_from_bed(bedfile, linewise)

        assert (constraints is not None) is expected

        if constraints is None:
            assert caplog.text != ""
        else:
            if not isinstance(constraints, Exception):
                assert type(constraints) == defaultdict
                for entry in constraints:
                    for cons in constraints[entry]:
                        assert "|" in cons
                        assert "-" in cons
                assert caplog.text == ""
    except ValueError:
        pass


@pytest.mark.parametrize(
    "dir_path",
    [
        "Foo",
        "bla/Foo",
    ],
)
def test_make_outdir(dir_path, caplog, tmp_path):
    dir_path = os.path.join(tmp_path, dir_path)
    outdir = make_outdir(dir_path)
    assert os.path.exists(outdir)
    assert caplog.text == ""


@pytest.fixture()
def gene_coords_string():
    data = "chr1\t50\t100\tgene1\t.\t+\n" "chr3\t150\t200\tgene1\t.\t-"
    return data


@pytest.fixture()
def gene_coords_bed(gene_coords_string):
    file = os.path.join(TMP_TEST_DIR, "gene_coords.bed")
    with open(file, "w") as handle:
        handle.write(gene_coords_string)
    yield file
    os.remove(file)


@pytest.fixture()
def gene_coords_bed(gene_coords_string):
    file = os.path.join(TMP_TEST_DIR, "gene_coords.bed")
    with open(file, "w") as handle:
        handle.write(gene_coords_string)
    yield file
    os.remove(file)


@pytest.fixture()
def gene_coords_bed_gz(gene_coords_string):
    file = os.path.join(TMP_TEST_DIR, "gene_coords.bed.gz")
    with gzip.open(file, "w") as handle:
        handle.write(bytes(gene_coords_string, encoding="utf-8"))
    yield file
    os.remove(file)


@pytest.mark.parametrize("gene_coords", ["gene_coords_bed", "gene_coords_bed_gz"])
def test_parse_annotation_bed(gene_coords, request, caplog):
    gene_coords = request.getfixturevalue(gene_coords)
    annotations = parse_annotation_bed(gene_coords)
    assert type(annotations) == defaultdict
    for entry in annotations:
        for annotation in annotations[entry]:
            assert len(annotation.split("|")) == 2
            coords, strand = annotation.split("|")
            assert len(coords.split("-")) == 2


def test_parse_annotation_bed_file_not_found():
    with pytest.raises(FileNotFoundError):
        annotations = parse_annotation_bed("nonsense")
