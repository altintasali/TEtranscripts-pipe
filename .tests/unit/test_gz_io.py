"""gz_io picks its opener from the filename extension alone."""
import gzip

import gz_io


def test_plain_roundtrip(tmp_path):
    p = tmp_path / "x.tsv"
    with gz_io.open_write(p) as fh:
        fh.write("a\tb\n")
    assert p.read_text() == "a\tb\n"
    with gz_io.open_read(p) as fh:
        assert fh.read() == "a\tb\n"


def test_gz_roundtrip(tmp_path):
    p = tmp_path / "x.tsv.gz"
    with gz_io.open_write(p) as fh:
        fh.write("a\tb\n")
    # actually compressed, not just named .gz
    assert p.read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(p, "rt") as fh:
        assert fh.read() == "a\tb\n"
    with gz_io.open_read(p) as fh:
        assert fh.read() == "a\tb\n"


def test_accepts_path_objects_and_strings(tmp_path):
    p = tmp_path / "x.tsv.gz"
    with gz_io.open_write(str(p)) as fh:
        fh.write("s\n")
    with gz_io.open_read(p) as fh:
        assert fh.read() == "s\n"
