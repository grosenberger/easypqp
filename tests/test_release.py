"""Release validation tests do not require the scientific runtime dependencies."""

import io
from pathlib import Path
import runpy
import tarfile
import zipfile

import pytest


check_release = runpy.run_path(
    str(Path(__file__).resolve().parents[1] / "scripts" / "check_release.py")
)["check_release"]


def make_distributions(tmp_path, version="0.1.62", *, database=True, sdist_version=None):
    metadata = f"Name: easypqp\nVersion: {version}\n".encode()
    with zipfile.ZipFile(tmp_path / "easypqp.whl", "w") as wheel:
        wheel.writestr("easypqp.dist-info/METADATA", metadata)
        if database:
            wheel.writestr("easypqp/data/unimod.xml", "<unimod/>")
    with tarfile.open(tmp_path / "easypqp.tar.gz", "w:gz") as sdist:
        metadata = f"Name: easypqp\nVersion: {sdist_version or version}\n".encode()
        info = tarfile.TarInfo("easypqp/PKG-INFO")
        info.size = len(metadata)
        sdist.addfile(info, io.BytesIO(metadata))


@pytest.mark.parametrize("tag", ["0.1.62", "v0.1.62", "0.1.62rc1"])
def test_matching_release(tmp_path, tag):
    make_distributions(tmp_path, tag.removeprefix("v"))
    check_release(tag, tmp_path)


@pytest.mark.parametrize("tag", [
    "", "master", "refs/heads/master", "0.1.062", "0.1.62+local", "0.1.62.dev1",
    "$(touch /tmp/should-not-execute)",
])
def test_invalid_release_tag(tmp_path, tag):
    with pytest.raises(ValueError, match="release tag|Release tag"):
        check_release(tag, tmp_path)


@pytest.mark.parametrize("wheel_version,sdist_version", [
    ("0.1.60", "0.1.60"),
    ("0.1.62", "0.1.60"),
    ("0.1.62.dev1+g123", "0.1.62.dev1+g123"),
])
def test_mismatched_release(tmp_path, wheel_version, sdist_version):
    make_distributions(tmp_path, wheel_version, sdist_version=sdist_version)
    with pytest.raises(ValueError, match="does not match"):
        check_release("0.1.62", tmp_path)


def test_missing_database(tmp_path):
    make_distributions(tmp_path, database=False)
    with pytest.raises(ValueError, match="UniMod"):
        check_release("0.1.62", tmp_path)


def test_extra_artifacts(tmp_path):
    make_distributions(tmp_path)
    (tmp_path / "old-release.whl").touch()
    with pytest.raises(ValueError, match="exactly one"):
        check_release("0.1.62", tmp_path)


def test_missing_artifacts(tmp_path):
    with pytest.raises(ValueError, match="exactly one"):
        check_release("0.1.62", tmp_path)
