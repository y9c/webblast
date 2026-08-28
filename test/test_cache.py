#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the content-addressed result cache."""

from webblast.cache import cache_key, get, put


def test_cache_roundtrip(tmp_path):
    key = cache_key("blastn", "nt", ">x\nACGT", {"max_num_seq": 5})
    assert key == cache_key("blastn", "nt", ">x\nACGT", {"max_num_seq": 5})  # stable
    assert key != cache_key("blastn", "nt", ">x\nACGT", {"max_num_seq": 10})

    assert get(key, tmp_path) is None  # miss
    path = put(key, b"hello", tmp_path)
    assert path.endswith(".json")
    assert get(key, tmp_path) == b"hello"  # hit


def test_cache_key_differs_by_query():
    assert cache_key("blastn", "nt", "ACGT") != cache_key("blastn", "nt", "ACGA")
