import pytest
import os

@pytest.fixture()
def sliding_window():

    return 1


def test_data_available():
    assert os.path.exists("test.fa.gz"), "Test data not available"