#!/usr/bin/env python3
"""One-time data preparation: use a local/private LLM gateway to translate
species missing from the dictionary, then bake the results into the bundled
data. This is an OFFLINE maintainer step — the installed package never calls AI.

The API key is taken ONLY from environment variables; it is never stored in the
repository.

Usage (set the secret in your environment, not in the command line):
    export LLM_BASE="https://your-gateway/v1"
    export LLM_MODEL="deepseek-chat"
    export LLM_KEY="sk-..."
    python scripts/prepare_species.py --input names.txt

    # dry run: only show how many names would be translated
    python scripts/prepare_species.py --input names.txt --report-only
"""

from __future__ import annotations

import argparse
import lzma
import os
import pickle
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from webblast.llm import LLMConfig, batch_translate  # noqa: E402

DICT_FILE = ROOT / "webblast" / "data" / "species.zh.pkl.xz"


def load_dict() -> dict:
    with lzma.open(DICT_FILE, "rb") as fh:
        return pickle.load(fh)


def save_dict(d: dict) -> None:
    with lzma.open(DICT_FILE, "wb") as fh:
        pickle.dump(d, fh, protocol=5)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="file of scientific names, one per line")
    ap.add_argument("--batch", type=int, default=200)
    ap.add_argument("--report-only", action="store_true",
                    help="only report coverage; don't call the model")
    ap.add_argument("--out-tsv", default=None, help="also write a TSV of new names")
    args = ap.parse_args()

    d = load_dict()
    names = [n.strip() for n in open(args.input, encoding="utf-8") if n.strip()]
    missing = [n for n in dict.fromkeys(names) if n not in d]
    print(f"bundled dict: {len(d):,} entries ; requested {len(names):,} ; missing {len(missing):,}")

    if not missing:
        print("nothing missing — nothing to translate.")
        return
    if args.report_only:
        print("report-only: would translate", len(missing), "names")
        return

    # Gateway secret is supplied via the environment for this one-off run; it is
    # never stored in the repository.
    config = LLMConfig(
        base_url=os.environ.get("LLM_BASE", ""),
        model=os.environ.get("LLM_MODEL", ""),
        api_key=os.environ.get("LLM_KEY", ""),
        batch_size=args.batch,
    )
    if not (config.base_url and config.model):
        print("set LLM_BASE and LLM_MODEL (and LLM_KEY) in the environment.")
        return

    print(f"translating {len(missing):,} names via {config.model} @ {config.base_url} ...")
    result = batch_translate(missing, config)
    added = {k: v for k, v in result.items() if v}
    print(f"got {len(added):,} translations, {len(result) - len(added):,} none/empty")

    before = len(d)
    d.update(added)
    save_dict(d)
    print(f"bundled dict updated: {before:,} -> {len(d):,}  ({DICT_FILE})")

    if args.out_tsv:
        with open(args.out_tsv, "w", encoding="utf-8") as fh:
            for k, v in added.items():
                fh.write(f"{k}\t{v}\n")
        print(f"wrote new-names TSV: {args.out_tsv}")

    for k in list(added)[:10]:
        print(f"  {k} -> {added[k]}")


if __name__ == "__main__":
    main()
