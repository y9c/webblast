#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for record reading (FASTA + bare sequence)."""

import io
from unittest import mock

from webblast.utils import _bare_sequence, _parse_fasta_text, read_records


def test_bare_sequence_collapses():
    assert _bare_sequence("ACGT\nACGT  acgt") == "ACGTACGTacgt"


def test_parse_fasta_text():
    recs = _parse_fasta_text(">a\nACGT\n>b\nTTTT\n", None)
    assert recs == [("a", "ACGT"), ("b", "TTTT")]


def _stdin(text):
    return mock.patch("sys.stdin", io.StringIO(text))


def test_read_bare_sequence_from_stdin():
    with _stdin("AGTCAAAACCACAA"):
        recs = read_records([])
    assert recs == [("seq", "AGTCAAAACCACAA")]


def test_read_fasta_from_stdin():
    with _stdin(">a\nACGT\n>b\nTTTT\n"):
        recs = read_records([])
    assert recs == [("a", "ACGT"), ("b", "TTTT")]


def test_read_empty_raises():
    import pytest

    with _stdin("   \n\n"):
        with pytest.raises(ValueError):
            read_records([])
