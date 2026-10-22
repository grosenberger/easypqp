"""Reject untagged, mismatched, or incomplete distributions before publishing."""

import argparse
from email.parser import BytesParser
from pathlib import Path
import tarfile
import zipfile

from packaging.version import InvalidVersion, Version


def check_release(tag, dist_dir):
    raw_version = tag.removeprefix("v")
    try:
        version = Version(raw_version)
    except InvalidVersion as exc:
        raise ValueError(f"Invalid release tag: {tag!r}") from exc
    if str(version) != raw_version or version.local is not None or version.is_devrelease:
        raise ValueError(f"Release tag must be a canonical public version: {tag!r}")

    paths = sorted(Path(dist_dir).iterdir())
    wheels = [p for p in paths if p.name.endswith(".whl")]
    sdists = [p for p in paths if p.name.endswith(".tar.gz")]
    if len(paths) != 2 or len(wheels) != 1 or len(sdists) != 1:
        raise ValueError("Expected exactly one wheel and one source distribution")

    with zipfile.ZipFile(wheels[0]) as wheel:
        metadata = [n for n in wheel.namelist() if n.endswith(".dist-info/METADATA")]
        if len(metadata) != 1 or "easypqp/data/unimod.xml" not in wheel.namelist():
            raise ValueError("Wheel is missing metadata or the packaged UniMod database")
        wheel_metadata = BytesParser().parsebytes(wheel.read(metadata[0]))

    with tarfile.open(sdists[0], "r:gz") as sdist:
        metadata = [
            m for m in sdist.getmembers()
            if m.name.count("/") == 1 and m.name.endswith("/PKG-INFO")
        ]
        if len(metadata) != 1:
            raise ValueError("Source distribution is missing package metadata")
        sdist_metadata = BytesParser().parsebytes(sdist.extractfile(metadata[0]).read())

    for metadata in (wheel_metadata, sdist_metadata):
        if metadata["Name"] != "easypqp" or metadata["Version"] != raw_version:
            raise ValueError(
                f"Tag {tag!r} does not match distribution "
                f"{metadata['Name']} {metadata['Version']}"
            )
    print(f"Validated wheel and source distribution for easypqp {version}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag")
    parser.add_argument("dist_dir", type=Path)
    args = parser.parse_args()
    try:
        check_release(args.tag, args.dist_dir)
    except ValueError as exc:
        parser.exit(1, f"{exc}\n")
