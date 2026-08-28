#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Batch species-name translation via an OpenAI-compatible local LLM gateway.

Used only as a *supplement* for species that have no entry in the bundled
dictionary. To keep it fast and cheap, names are sent in large batches (a few
thousand per request) rather than one-at-a-time, and every returned name is
flagged as an **AI guess (unverified)**.

Works with any OpenAI-compatible endpoint (vLLM, LiteLLM, LM Studio, a local
DeepSeek gateway, ...) via ``POST {base_url}/chat/completions``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import requests

DEFAULT_LLM_BASE = "http://localhost:8000/v1"
DEFAULT_LLM_MODEL = "deepseek-chat"


@dataclass
class LLMConfig:
    base_url: str = DEFAULT_LLM_BASE
    model: str = DEFAULT_LLM_MODEL
    api_key: str = ""
    batch_size: int = 500          # names per request (you can raise to ~2000)
    timeout: float = 300.0


def _chat_completions_url(base_url: str) -> str:
    return base_url.rstrip("/") + "/chat/completions"


def _headers(config: LLMConfig) -> dict:
    h = {"Content-Type": "application/json"}
    if config.api_key:
        h["Authorization"] = f"Bearer {config.api_key}"
    return h


def _build_prompt(names: Sequence[str]) -> str:
    # Tab-delimited output is the most reliable across local OpenAI-compatible
    # gateways (DeepSeek/vLLM can garble `json_object` responses).
    lines = "\n".join(names)
    return (
        "You translate scientific binomial species names into their standard "
        "Chinese common names (中文名). Reply with exactly ONE line per input "
        "name, in this format: <input name><TAB><Chinese name>. If a species has "
        "no standard Chinese name, use: <input name><TAB>N/A. Do not add any other "
        "text.\n\n"
        f"{lines}"
    )


def _parse_response(content: str, names: Sequence[str]) -> Dict[str, Optional[str]]:
    """Parse the model reply into {name: chinese|None}, tolerant of formats."""
    # 1) JSON object
    try:
        data = json.loads(content)
        if isinstance(data, dict):
            out = {}
            for n in names:
                v = data.get(n)
                out[n] = v.strip() if isinstance(v, str) and v.strip() and v.strip().upper() not in (
                    "N/A", "NA", "NONE", "NULL", "EMPTY") else None
            return out
    except (ValueError, TypeError):
        pass

    # 2) tab-delimited lines: <name>\t<chinese>, or <name>\tN/A
    out = {}
    for n in names:
        out[n] = None
    for line in content.splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            v = v.strip()
            if k.strip() in out and v and v.upper() not in ("N/A", "NA", "NONE", "NULL", "EMPTY"):
                out[k.strip()] = v
        elif ":" in line and line.count(":") == 1:
            k, v = line.split(":", 1)
            v = v.strip()
            if k.strip() in out and v and v.upper() not in ("N/A", "NA", "NONE", "NULL", "EMPTY"):
                out[k.strip()] = v
    return out


def batch_translate(names: Sequence[str], config: LLMConfig = LLMConfig()) -> Dict[str, Optional[str]]:
    """Translate many names in batches. Returns {name: chinese|None} for all input.

    Raises the last HTTP error if the endpoint is unreachable, so callers can
    decide whether to fail or fall back silently.
    """
    uniq = list(dict.fromkeys(n for n in names if n))  # dedupe, preserve order
    out: Dict[str, Optional[str]] = {n: None for n in uniq}
    last_err: Optional[Exception] = None

    for i in range(0, len(uniq), config.batch_size):
        chunk = uniq[i : i + config.batch_size]
        try:
            resp = requests.post(
                _chat_completions_url(config.base_url),
                headers=_headers(config),
                json={
                    "model": config.model,
                    "messages": [
                        {"role": "system", "content": "You are a precise species-name translator."},
                        {"role": "user", "content": _build_prompt(chunk)},
                    ],
                    "temperature": 0.0,
                    "max_tokens": 32768,
                },
                timeout=config.timeout,
            )
            resp.raise_for_status()
            content = resp.json()["choices"][0]["message"]["content"]
            out.update(_parse_response(content, chunk))
            last_err = None
        except Exception as exc:  # pragma: no cover - network dependent
            last_err = exc
            # leave this chunk as None; record and continue
    if last_err and not any(v for v in out.values()):
        raise last_err
    return out


def available(config: LLMConfig = LLMConfig()) -> bool:
    """True if the gateway responds."""
    try:
        r = requests.get(config.base_url.rstrip("/") + "/models", timeout=3)
        return r.status_code == 200
    except Exception:
        return False
