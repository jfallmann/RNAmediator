import os
from io import StringIO
from tempfile import TemporaryDirectory, NamedTemporaryFile

import pytest
from Bio import SeqIO
from collections import defaultdict
from RIssmed.Tweaks.RIssmed import (
    SequenceSettings,
    Constraint,
    get_gene_coords,
    set_run_settings_dict,
    read_constraints,
    add_rissmed_constraint,
    preprocess,
    expand_pl_window,
    localize_pl_window,
)

TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
TESTDATAPATH = os.path.join(TESTFOLDER, "testdata")
TMP_DIR = TemporaryDirectory()
TMP_TEST_DIR = TMP_DIR.name


@pytest.fixture
def single_sequence_fasta():
    path = os.path.join(TESTDATAPATH, "test.fa")
    return path


@pytest.fixture
def sequence_string():
    return "AAATTTTTGGGGUUUUUUUCCCCTTTTTttt"


@pytest.fixture
def stringio_fasta():
    fasta = (
        ">ENSG00000270742::chr1:61124732-61125202(+)\n"
        "AAAATTTTTTTTTTTTTTTGGGGGGGGGGGGGGCCCCCCCCCCGGGGGGGGGGGGGGGGGGAA"
    )
    file = StringIO(fasta)
    yield file
    file.close()


@pytest.fixture
def command_line_constraint():
    return "5-20"


@pytest.fixture
def command_line_mutate_constraint():
    return "50-51:A"


@pytest.fixture
def paired_constraints():
    return os.path.join(TESTDATAPATH, "paired_constraints.bed")


@pytest.fixture()
def sliding_constraint():
    return "sliding"


@pytest.fixture()
def scanning_constraint():
    return "scanning"


@pytest.fixture
def genomic_coords_file():
    genomic_coords = "chr1\t200\t1200\n"
    tmp_coords = NamedTemporaryFile("wt", suffix=".bed", delete=False)
    tmp_coords.write(genomic_coords)
    tmp_coords.close()
    yield tmp_coords.name
    os.remove(tmp_coords.name)


@pytest.fixture()
def single_bedfile():
    bed = """chr1	26	35	ENSG00000270742	u	+"""
    tmp_bedfile = NamedTemporaryFile("wt", suffix=".bed", delete=False)
    tmp_bedfile.write(bed)
    tmp_bedfile.close()
    yield tmp_bedfile.name
    os.remove(tmp_bedfile.name)


@pytest.fixture()
def single_longmutate_bedfile():
    bed = """chr1	26	35	ENSG00000270742	ACTCAGACT	+"""
    tmp_bedfile = NamedTemporaryFile("wt", suffix=".bed", delete=False)
    tmp_bedfile.write(bed)
    tmp_bedfile.close()
    yield tmp_bedfile.name
    os.remove(tmp_bedfile.name)


@pytest.fixture()
def single_snp_bedfile():
    bed = """chr1	26	27	ENSG00000270742	A	+"""
    tmp_bedfile = NamedTemporaryFile("wt", suffix=".bed", delete=False)
    tmp_bedfile.write(bed)
    tmp_bedfile.close()
    yield tmp_bedfile.name
    os.remove(tmp_bedfile.name)


@pytest.fixture
def multi_bedfile():
    return os.path.join(TESTDATAPATH, "test_constraints.bed")


@pytest.fixture
def multi_bedfile_gz():
    return os.path.join(TESTDATAPATH, "test_constraints.bed.gz")


@pytest.fixture
def basic_sequence_settings(sequence_string):
    seq: SeqIO.SeqRecord = SeqIO.SeqRecord(sequence_string)
    return SequenceSettings(seq)


@pytest.fixture
def na_strand_sequence_settings(sequence_string):
    seq: SeqIO.SeqRecord = SeqIO.SeqRecord(sequence_string)
    return SequenceSettings(seq, strand="na")


@pytest.fixture
def basic_constraint():
    return Constraint(10, 50, "+")


@pytest.fixture
def minus_constraint():
    return Constraint(5, 20, "-")


@pytest.fixture
def mutate_constraint():
    return Constraint(10, 1, "+", "G")


@pytest.fixture
def mutate_minus_constraint():
    return Constraint(5, 6, "-", "C")


@pytest.fixture
def basic_constraintlist():
    cons_a = (Constraint(10, 20, "+"),)
    cons_b = (Constraint(50, 70, "+"), Constraint(10, 20, "+"))
    return [cons_a, cons_b]


@pytest.fixture
def basic_minus_constraintlist():
    cons_a = (Constraint(10, 20, "-"),)
    cons_b = (Constraint(50, 70, "-"), Constraint(10, 20, "-"))
    return [cons_a, cons_b]


@pytest.fixture
def mutate_constraintlist():
    cons_a = (Constraint(10, 20, "+", "A"),)
    cons_b = (Constraint(50, 70, "+", "U"), Constraint(10, 20, "+", "T"))
    return [cons_a, cons_b]


@pytest.fixture
def mutate_minus_constraintlist():
    cons_a = (Constraint(10, 20, "-", "A"),)
    cons_b = (Constraint(50, 70, "-", "U"), Constraint(10, 20, "-", "T"))
    return [cons_a, cons_b]


@pytest.mark.parametrize(
    "seq, strand, constraintlist",
    [
        (SeqIO.SeqRecord("AAAACCCCGGGGGUUUU"), "+", "basic_constraintlist"),
        (SeqIO.SeqRecord("AAAATTTTTTTT"), "-", "basic_minus_constraintlist"),
    ],
)
def test_sequence_settings_init_with_strand(seq, strand, constraintlist, request):
    constraintlist = request.getfixturevalue(constraintlist)
    seq_settings = SequenceSettings(seq, constrainlist=constraintlist, strand=strand)
    assert "T" not in seq_settings.sequence_record
    for entry in seq_settings.constrainlist:
        for cons in entry:
            assert cons.strand == strand


def test_sequence_settings_auto_strand(sequence_string, caplog):
    seq = SeqIO.SeqRecord(sequence_string)
    seq_settings = SequenceSettings(seq, strand="foo")
    assert seq_settings.strand == "na"
    assert "strand value automatically set to na" in caplog.text


@pytest.mark.parametrize(
    "seqsettings, constraint",
    [
        ("basic_sequence_settings", "basic_constraint"),
        ("na_strand_sequence_settings", "minus_constraint"),
    ],
)
def test_sequence_settings_add_constraint(constraint, seqsettings, request):
    constraint = request.getfixturevalue(constraint)
    seqsettings = request.getfixturevalue(seqsettings)
    seqsettings.add_constraints([constraint])


@pytest.mark.parametrize(
    "seqsettings, constraint",
    [
        ("basic_sequence_settings", "mutate_constraint"),
        ("na_strand_sequence_settings", "mutate_minus_constraint"),
    ],
)
def test_sequence_settings_add_mutate_constraint(constraint, seqsettings, request):
    constraint = request.getfixturevalue(constraint)
    seqsettings = request.getfixturevalue(seqsettings)
    seqsettings.add_constraints([constraint])


@pytest.mark.parametrize(
    "seqsettings, constraint",
    [
        ("basic_sequence_settings", "minus_constraint"),
    ],
)
def test_sequence_settings_add_constraints_assertions(seqsettings, constraint, request, capsys):
    constraint = request.getfixturevalue(constraint)
    seqsettings = request.getfixturevalue(seqsettings)
    seqsettings.add_constraints([constraint])
    out, err = capsys.readouterr()
    assert "AssertionError" in out


def test_constraint(basic_constraint):
    string_repr = str(basic_constraint)
    coords, strand = string_repr.split("|")
    start, end = [int(x) for x in coords.split("-")]
    assert start == basic_constraint.start
    assert end == basic_constraint.end


@pytest.fixture
def gene_coords_dict():
    gcd = {
        "gene1": ["8000-19000|+"],
        "gene2": ["5000-19000|+"],
        "gene3": ["200-700|-"],
        "gene4": ["1050-3080|+"],
    }
    return gcd


@pytest.fixture
def none_fixture():
    return None


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("goi", ["gene1", "gene2", "gene3", "gene4"])
def test_get_gene_coords_with_annotation(gene_coords_dict, goi, strand, caplog):
    expected_strand = gene_coords_dict[goi][0].split("|")[1]
    genomic_start, genomic_end, genomic_strand, value = get_gene_coords(gene_coords_dict, goi, strand)
    assert expected_strand == genomic_strand
    if expected_strand != strand:
        assert "Strand values differ" in caplog.text
    else:
        assert "" == caplog.text


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("goi", ["gene1", "gene2", "gene3", "gene4"])
def test_get_gene_coords_default(none_fixture, goi, strand, caplog):
    genomic_start, genomic_end, genomic_strand, value = get_gene_coords(none_fixture, goi, "+")
    assert "No coords found for gene" in caplog.text
    assert genomic_start == 0
    assert genomic_end == 0
    assert genomic_strand == "."


def test_get_gene_coords_missing_entry(gene_coords_dict, caplog):
    genomic_start, genomic_end, genomic_strand, value = get_gene_coords(gene_coords_dict, "foo", "+")
    assert "No coords found for gene" in caplog.text
    assert genomic_start == 0
    assert genomic_end == 0
    assert genomic_strand == "."


@pytest.mark.parametrize(
    "sequence, constrain, genomic_coords, conslength, constraintype",
    [
        ("sequence_string", "3-5", None, 1, "hard"),
        ("sequence_string", "3-4|A", None, 1, "mutate"),
        ("single_sequence_fasta", "ono,10-15", "genomic_coords_file", 1, "hard"),
        ("stringio_fasta", "single_bedfile", "genomic_coords_file", 1, "hard"),
        ("stringio_fasta", "sliding", "genomic_coords_file", 3, "hard"),
    ],
)
def test_set_run_settings_dict(sequence, constrain, conslength, constraintype, genomic_coords, request, caplog):
    sequence = request.getfixturevalue(sequence)
    try:
        constrain = request.getfixturevalue(constrain)
    except pytest.FixtureLookupError:
        pass
    genomic_coords = request.getfixturevalue(genomic_coords) if genomic_coords is not None else ""
    constraintype = constraintype if constraintype is not None else "hard"
    run_settings = set_run_settings_dict(sequence, constrain, conslength, genomic_coords, constraintype)
    for entry in run_settings:

        assert isinstance(run_settings[entry], SequenceSettings)


@pytest.mark.parametrize("linewise", [True, False])
@pytest.mark.parametrize(
    "constraints",
    [
        "multi_bedfile",
        "multi_bedfile_gz",
        "command_line_constraint",
        "paired_constraints",
        "sliding_constraint",
        "scanning_constraint",
    ],
)
def test_read_constraints(constraints, request, linewise):
    constraints = request.getfixturevalue(constraints)
    constraintlist = read_constraints(constraints, linewise)
    if linewise and constraints != "sliding" and constraints != "scanning":
        assert type(constraintlist) == defaultdict
        assert "lw" in constraintlist
    else:
        if os.path.isfile(constraints):
            assert type(constraintlist) == defaultdict
        else:
            assert type(constraintlist) == list


@pytest.mark.parametrize("cons1", ["5-10|+"])
@pytest.mark.parametrize(
    "cons2, expected_length",
    [
        ("5-10|+:10-20|+", 2),
        ("7-8:12-13", 2),
        ("7-8|+", 1),
    ],
)
def test_add_rissmed_constraint(cons1, cons2, expected_length, sequence_string):
    seq = SeqIO.SeqRecord(sequence_string)
    run_settings_dict = dict()
    add_rissmed_constraint(run_settings_dict, cons1, seq)
    add_rissmed_constraint(run_settings_dict, cons2, seq)
    assert len(run_settings_dict) == 1
    for entry in run_settings_dict:
        conslist = run_settings_dict[entry].constrainlist
        assert len(conslist) == 2
        assert len(conslist[1]) == expected_length


@pytest.mark.parametrize("outdir", ["", os.path.join(TMP_TEST_DIR, "preprocess")])
def test_preprocess(single_sequence_fasta, single_bedfile, outdir):
    run_settings, outdir = preprocess(single_sequence_fasta, single_bedfile, 1, "hard", outdir, "")
    assert os.path.isdir(outdir)
    assert isinstance(run_settings, dict)


@pytest.mark.parametrize(
    "start,end,window,multiplyier,seqlen,expected",
    [
        (20, 25, 50, 2, 150, [1, 125]),
        (100, 125, 70, 1, 300, [30, 195]),
        (100, 125, 70, 2, 150, [1, 150]),
    ],
)
def test_expand_pl_window(start, end, window, multiplyier, seqlen, expected):
    new_start, new_end = expand_pl_window(start, end, window, multiplyier, seqlen)
    expected_start = expected[0]
    expected_end = expected[1]
    assert new_start == expected_start
    assert new_end == expected_end


@pytest.mark.parametrize(
    "start,end,window,multiplyier,seqlen,expected",
    [
        (20, 25, 50, 2, 150, [1, 125]),
        (100, 125, 70, 1, 300, [30, 195]),
        (100, 125, 70, 2, 150, [1, 150]),
    ],
)
def test_localize_pl_window(start, end, window, multiplyier, seqlen, expected):
    local_start, local_end = localize_pl_window(start, end, window, seqlen, multiplyier)
