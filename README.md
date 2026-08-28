# webblast

[![Pypi Releases](https://img.shields.io/pypi/v/webblast.svg)](https://pypi.python.org/pypi/webblast)
[![Downloads](https://pepy.tech/badge/webblast)](https://pepy.tech/project/webblast)

**Fast NCBI BLAST searches from the command line.** `webblast` submits
nucleotide / protein queries to NCBI BLAST through the **same interactive path
the website uses** (a real NCBI session + browser `User-Agent` + the full web
parameter set), then polls **early and adaptively** — so a typical search is
done in **a few seconds**, not the ~30s that NCBI's conservative `RTOE` estimate
suggests.

It parses the compact `JSON2_S` report (instead of multi-megabyte HTML) into
clean, machine-readable records.

## Install

```bash
pip install webblast
```

## Quick start

```bash
# search a FASTA file against nt with megablast (default)
webblast -p megablast -d nt query.fa

# TSV output, limit to 50 target sequences
webblast -d nt -f tsv --max-num-seq 50 query.fa

# protein search
webblast -p blastp -d nr protein.fa

# read a pasted FASTA from stdin
echo '>seq
AGTCAAAACCACAATGAGATACCATCTCATGTCAGTCAGAATGGCTATTACTAAAAA' | webblast --limit 1
```

Runs as a Python library too:

```python
from webblast import BlastClient, parse_blast_json

client = BlastClient(program="megablast", database="nt", email="you@example.org")
report = client.run(">seq\nACGTACGTACGT", format_type="JSON2_S")
parsed = parse_blast_json(report)
print(parsed.program, parsed.version, parsed.total_hits)
```

## Why it's fast

1. **Interactive web queue.** NCBI gives priority to interactive (website) users
   and deprioritizes automated API traffic. `BlastClient` reproduces the browser's
   submission — session cookie, browser User-Agent, and web parameters
   (`PAGE`, `BLAST_PROGRAMS`, `MEGABLAST=on`, ...) — so jobs land on the fast
   queue.
2. **Don't sleep the RTOE.** NCBI returns a very conservative `RTOE` (often ~30s)
   that is a poor estimate for short queries, which usually finish in **3–10s**.
   `webblast` polls early and backs off geometrically, catching the result as
   soon as it's ready.
3. **Compact JSON.** Results are read back as `JSON2_S` (hundreds of KB, not
   multi-megabyte HTML), so parsing is light.

## Output formats

| `--format` | Description |
|------------|-------------|
| `text` (default) | Rich terminal table: accession, description, E-value, %identity, alignment length, bit score, 中文名 |
| `tsv` | Clean tab-delimited rows (query, accession, title, evalue, pident, align_len, bitscore, zh) |
| `json` | Compact JSON with per-query/hit/hsp fields (+ `zh` / `sciname`) |
| `fasta` | Subject segments of the top HSP per hit (pseudo-FASTA) |

## 中文名 (Chinese species names)

Each hit gets its Chinese common name from a bundled dictionary
(`webblast/data/species.zh.pkl.xz`, lzma-compressed) covering **~489,000 species**.
It's loaded lazily (once) and looked up in O(1), so it's offline, deterministic
and fast — no model, no network calls at runtime.

- `--no-translate` disables the 中文名 column (it's on by default).

To grow the dictionary, drop extra `scientific_name<TAB>中文名` rows into a file
and point `webblast` at it via the `WEBLAST_SPECIES_TSV` env var (or
`~/.config/webblast/species.tsv`); those rows override the bundled ones.

## Options

```
-p, --program       megablast (default), discontiguous-megablast, blastn, blastp,
                    quickblastp, blastx, tblastn, tblastx, rpsblast
-d, --database      nt (default), nr, refseq_rna, refseq_protein, swissprot, ...
-o, --out FILE      write output to a file
--max-num-seq N     max target sequences per query (default 100)
--expect E          E-value cutoff (e.g. 1e-5)
--email ADDR        contact email (NCBI usage policy)
--api-key KEY       NCBI API key (raises rate limit)
--cache/--no-cache  on-disk result cache (default on)
--limit N           only use the first N query records
--wait-timeout S    max wait for the search (default 1200s)
```

## Caching

Identical searches (keyed on program + database + options + query) are cached on
disk under `blast_cache/`, so re-running the same query is instant.

## Input formats

FASTA / FASTQ / (optionally) BAM/SAM via `pysam` (`pip install webblast[bam]`),
from files or stdin.

## Notes & fair use

* NCBI is a shared resource. Keep submissions modest, use `--email`, and
  consider off-peak hours for large batches.
* Running tens of thousands of searches is better served by a
  local/cloud BLAST install (the URL API only suits small, occasional jobs).

## Development

```bash
pip install -e .[dev]
pytest
```
