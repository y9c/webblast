#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the CLI renderers (text/tsv/fasta/json) using a real fixture."""

import json
from pathlib import Path

from webblast.cli import _render_fasta, _render_tsv
from webblast.client import PROGRAMS
from webblast.parse import parse_blast_json

FIXTURE = Path(__file__).parent / "fixtures" / "blastn_json2s.json"


def _report():
    return parse_blast_json(FIXTURE.read_bytes())


def test_tsv_render():
    out = _render_tsv(_report())
    assert out
    lines = [l for l in out.splitlines() if l.strip()]
    # headerless tab-delimited rows
    assert "\t" in lines[0]
    assert len(lines) == _report().total_hits


def test_fasta_render():
    out = _render_fasta(_report())
    assert out.startswith(">")
    lines = out.splitlines()
    # fasta: even number of header/seq lines
    assert len(lines) % 2 == 0


def test_programs_and_databases_exist():
    from webblast.client import DATABASES, PROGRAMS

    assert "megablast" in PROGRAMS
    assert "nt" in DATABASES


def test_megablast_emits_megablast_param():
    from webblast.client import BlastClient

    # ensure megablast actually selects the megablast algorithm (MEGABLAST=on),
    # otherwise the server falls back to a slow sensitive blastn
    assert PROGRAMS["megablast"]["params"]["MEGABLAST"] == "on"
    client = BlastClient(program="megablast", database="nt", email="x@example.org")
    params = client._build_params(">x\nACGT", "JSON2_S")
    assert params["MEGABLAST"] == "on"
    assert params["BLAST_PROGRAMS"] == "megaBlast"
