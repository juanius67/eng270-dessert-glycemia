from __future__ import annotations
import sys
import ctypes
import pytest
from pathlib import Path
from unittest.mock import MagicMock, patch
from src import bindings
from src.bindings import (
    _shared_name,
    _candidate_commands,
    _hash_file,
    model_source_hash,
    set_params_from_dict,
    get_build_metadata,
    BergmanParams,
)

def test_shared_name():
    with patch("sys.platform", "linux"):
        assert _shared_name() == "libmodel.so"
    with patch("sys.platform", "darwin"):
        assert _shared_name() == "libmodel.dylib"
    with patch("sys.platform", "win32"):
        assert _shared_name() == "model.dll"

def test_candidate_commands():
    output = Path("output.so")

    with patch("sys.platform", "linux"):
        cmds = _candidate_commands(output)
        assert len(cmds) == 2
        assert "gcc" in cmds[0]
        assert "clang" in cmds[1]

    with patch("sys.platform", "darwin"):
        cmds = _candidate_commands(output)
        assert len(cmds) == 2
        assert "clang" in cmds[0]
        assert "gcc" in cmds[1]

    with patch("sys.platform", "win32"):
        cmds = _candidate_commands(output)
        assert len(cmds) == 2
        assert "cl" in cmds[0]
        assert "gcc" in cmds[1]

def test_hash_file(tmp_path):
    test_file = tmp_path / "test_hash.txt"
    content = b"test content"
    test_file.write_bytes(content)

    import hashlib
    expected_hash = hashlib.sha256(content).hexdigest()
    assert _hash_file(test_file) == expected_hash

@patch("src.bindings._hash_file")
def test_model_source_hash(mock_hash):
    mock_hash.return_value = "dummy_hash"
    assert model_source_hash() == "dummy_hash"
    mock_hash.assert_called_once()

def test_set_params_from_dict_error():
    # Save original lib
    original_lib = bindings.lib
    bindings.lib = None
    try:
        with pytest.raises(RuntimeError, match="C library not loaded"):
            set_params_from_dict({})
    finally:
        bindings.lib = original_lib

def test_set_params_from_dict_success():
    original_lib = bindings.lib
    mock_lib = MagicMock()
    bindings.lib = mock_lib

    params = {
        "S_G_min1": 1.0,
        "p2_min1": 2.0,
        "p3_min1_per_uU_mL": 3.0,
        "n_min1": 4.0,
        "Gb_mg_dL": 5.0,
        "Ib_uU_mL": 6.0,
    }

    try:
        set_params_from_dict(params)
        mock_lib.set_params.assert_called_once()
        args = mock_lib.set_params.call_args[0][0]
        assert isinstance(args, BergmanParams)
        assert args.S_G == 1.0
        assert args.p2 == 2.0
        assert args.p3 == 3.0
        assert args.n == 4.0
        assert args.Gb == 5.0
        assert args.Ib == 6.0
    finally:
        bindings.lib = original_lib

@patch("src.bindings.BUILD_METADATA", {"key": "val"})
@patch("src.bindings.model_source_hash", return_value="hash123")
@patch("sys.version", "3.10")
@patch("platform.platform", return_value="Linux-Test")
def test_get_build_metadata(mock_platform, mock_source_hash):
    meta = get_build_metadata()
    assert meta["key"] == "val"
    assert meta["source_hash"] == "hash123"
    assert meta["platform"] == "Linux-Test"
    assert meta["python"] == "3.10"
    assert "library" in meta
    assert "compiler" in meta
    assert "flags" in meta
    assert "command" in meta
