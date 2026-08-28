#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for parsing real NCBI BLAST JSON2_S output."""

import json
from pathlib import Path

import pytest

from webblast.parse import BlastReport, parse_blast_json

FIXTURE = Path(__file__).parent / "fixtures" / "blastn_json2s.json"


def test_parse_top_level():
    data = FIXTURE.read_bytes()
    rep = parse_blast_json(data)
    assert isinstance(rep, BlastReport)
    assert rep.program == "blastn"
    assert "BLASTN" in (rep.version or "")
    assert rep.db
    assert len(rep.queries) == 1


def test_parse_hit_hsp_fields():
    rep = parse_blast_json(FIXTURE.read_bytes())
    q = rep.queries[0]
    assert q.query_title == "probe"
    assert q.query_len == 57
    assert len(q.hits) > 0

    hit = q.hits[0]
    assert hit.accession
    assert hit.title
    assert hit.len > 0
    assert hit.hsps

    hsp = hit.hsps[0]
    assert hsp.evalue is not None
    assert hsp.pident is not None
    assert 0 <= hsp.pident <= 100
    assert hsp.align_len == 57
    assert hsp.qseq
    assert len(hsp.qseq) == len(hsp.hseq)


def test_parse_rejects_non_blast_payload():
    with pytest.raises(ValueError):
        parse_blast_json(json.dumps({"hello": "world"}))
