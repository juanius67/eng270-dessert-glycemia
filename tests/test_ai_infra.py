from __future__ import annotations
import json
import platform
import numpy as np
from src.ai_infra import collect_env_metadata, write_json

def test_collect_env_metadata():
    metadata = collect_env_metadata()
    assert isinstance(metadata, dict)
    assert metadata["python"] == platform.python_version()
    assert metadata["numpy"] == np.__version__
    assert metadata["os"] == platform.platform()
    assert "python_build" in metadata

def test_write_json(tmp_path):
    data = {"key": "value", "number": 123}
    test_file = tmp_path / "test.json"
    write_json(test_file, data)

    assert test_file.exists()
    with test_file.open("r", encoding="utf-8") as f:
        loaded_data = json.load(f)

    assert loaded_data == data

    # Check nested directory creation
    nested_file = tmp_path / "subdir" / "test.json"
    write_json(nested_file, data)
    assert nested_file.exists()
