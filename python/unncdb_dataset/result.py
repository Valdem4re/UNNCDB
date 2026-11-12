"""Save and reload dataset DataFrames with a JSON manifest."""
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def save_csv(frame: pd.DataFrame, path, *, overwrite=False):
    """Save a dataset frame and the manifest in frame.attrs."""
    manifest = dict(frame.attrs["manifest"])
    path = Path(path).expanduser().resolve()
    sidecar = path.with_suffix(".manifest.json")
    if path == sidecar:
        raise ValueError("CSV and manifest paths must differ")
    if not overwrite and (path.exists() or sidecar.exists()):
        raise FileExistsError(f"Dataset or manifest already exists: {path}")
    missing = set(manifest["feature_columns"] + [manifest["target_column"], "family_id"]) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing dataset columns: {sorted(missing)}")
    if frame[manifest["target_column"]].isna().any() or not frame[manifest["target_column"]].isin([0, 1]).all():
        raise ValueError("Target labels must be 0 or 1")
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".dataset-", dir=path.parent) as directory:
        tmp_csv = Path(directory) / "data.csv"
        tmp_manifest = Path(directory) / "manifest.json"
        frame.to_csv(tmp_csv, index=False)
        manifest["csv"] = {"filename": path.name, "sha256": _hash(tmp_csv), "rows": len(frame),
                           "columns": list(frame.columns), "missing_value": "empty field"}
        tmp_manifest.write_text(json.dumps(manifest, ensure_ascii=False, indent=2,
                                           default=_json_default) + "\n", encoding="utf-8")
        # Publish manifest last. A partial/mismatched pair is detected on reload.
        os.replace(tmp_csv, path)
        os.replace(tmp_manifest, sidecar)
    return path

def load_csv(path, *, verify=True) -> pd.DataFrame:
    """Load a dataset frame and restore its manifest."""
    path = Path(path).expanduser().resolve()
    manifest = json.loads(path.with_suffix(".manifest.json").read_text(encoding="utf-8"))
    if manifest.get("format_version") != 1:
        raise ValueError("Unsupported dataset manifest version")
    if verify and _hash(path) != manifest["csv"]["sha256"]:
        raise ValueError("CSV checksum does not match the manifest")
    frame = pd.read_csv(path, float_precision="round_trip", keep_default_na=False,
                        na_values=[""], dtype={"refcode": str, "family_id": str,
                                              "source_sha256": str, "cif_sha256": str,
                                              "structure_hash": str})
    if len(frame) != manifest["csv"]["rows"] or list(frame.columns) != manifest["csv"]["columns"]:
        raise ValueError("CSV shape/columns do not match the manifest")
    if not frame[manifest["target_column"]].isin([0, 1]).all():
        raise ValueError("Invalid target labels")
    frame.attrs["manifest"] = manifest
    frame.attrs["rejections"] = []
    return frame
