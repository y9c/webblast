#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple content-addressed cache for BLAST reports.

Results are stored as the raw JSON report bytes, keyed by a hash of the full
request (program + database + search options + query), so identical searches
are served instantly from disk.
"""

from __future__ import annotations

import hashlib
import json
import os
from typing import Optional, Union

DEFAULT_CACHE_DIR = "blast_cache"


def _serialise(value) -> str:
    if value is None:
        return ""
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True, default=str)
    return str(value)


def cache_key(program: str, database: str, query: str, options: Optional[dict] = None) -> str:
    """Build a stable md5 key for a search request."""
    parts = [_serialise(program), _serialise(database), _serialise(options or {}), query]
    return hashlib.md5("\x00".join(parts).encode("utf-8")).hexdigest()


def cache_path(key: str, cache_dir: str = DEFAULT_CACHE_DIR) -> str:
    return os.path.join(cache_dir, f"{key}.json")


def get(key: str, cache_dir: str = DEFAULT_CACHE_DIR) -> Optional[bytes]:
    """Return cached report bytes, or None."""
    path = cache_path(key, cache_dir)
    if os.path.exists(path):
        try:
            with open(path, "rb") as fh:
                return fh.read()
        except OSError:
            return None
    return None


def put(key: str, data: Union[bytes, bytearray, str], cache_dir: str = DEFAULT_CACHE_DIR) -> str:
    """Store report bytes and return the written path."""
    if isinstance(data, str):
        data = data.encode("utf-8")
    os.makedirs(cache_dir, exist_ok=True)
    path = cache_path(key, cache_dir)
    with open(path, "wb") as fh:
        fh.write(data)
    return path
