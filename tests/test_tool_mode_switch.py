"""The one contract that matters most tonight: swapping WILDTYPE_MODE
must not change which methods exist, only what they return. If this test
ever fails, the agent loop and mocked-mode work done tonight won't carry
over cleanly to real Proto calls tomorrow.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from wildtype.tools.base import get_toolset


def test_mock_toolset_shape():
    ts = get_toolset("mock")
    assert ts.mode == "mock"
    score = ts.scoring.score_missense("GPT", "MKV" * 20, 5, "A")
    assert score.gene == "GPT"
    struct = ts.structure.fold("GPT", "MKV" * 20)
    assert struct.plddt_mean > 0
    evo1 = ts.dna.score_structural("FGF4", "ACGT" * 20, "retrogene insertion")
    assert evo1.locus == "FGF4"
    hits = ts.literature.search("GPT dog")
    assert len(hits) > 0


def test_local_toolset_has_same_shape_as_mock():
    mock_ts = get_toolset("mock")
    local_ts = get_toolset("local")
    assert set(dir(mock_ts.scoring)) & {"score_missense", "functional_embedding"} == {
        "score_missense",
        "functional_embedding",
    }
    assert set(dir(local_ts.scoring)) & {"score_missense", "functional_embedding"} == {
        "score_missense",
        "functional_embedding",
    }


def test_unknown_mode_raises():
    import pytest

    with pytest.raises(ValueError):
        get_toolset("not-a-real-mode")
