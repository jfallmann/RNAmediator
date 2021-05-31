import os

TESTFOLDER = os.path.dirname(os.path.abspath(__file__))
PARPATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# sys.path.append(PARPATH)
from RIssmed.RNAtweaks.RNAtweaks import api_rnaplfold, cmd_rnaplfold
import pytest
import random


def randomsequence():
    random.seed(1)
    return "".join(random.choices(["A", "T", "U", "G"], k=350))

@pytest.mark.parametrize(
    "seq,window,span,constraint",
    [
        ("A" * 500, 70, 70, []),
        ("A" * 500, 70, 70, [("p", 20, 30)]),
        (randomsequence(), 70, 70, [("p", 20, 30)]),
    ]
)
def test_api_and_cmd_plfold(seq, window, span, constraint):
    cmd_result = cmd_rnaplfold(seq, window, span, constraint=constraint)
    api_result = api_rnaplfold(seq, window, span, constraint=constraint)
    assert api_result == cmd_result
