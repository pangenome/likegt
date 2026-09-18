# Genotyping from FASTA/FASTQ: implementation report and validation

**Branch:** `geno-fasta-input`
**Base:** `ekg/likegt@58b91ad` (main, includes PR #1 `pangenome/geno-external-bam`)
**Status:** all tests pass; real-data regression validation byte-identical to historical results.

## Summary

`likegt geno` previously required a BAM/CRAM aligned to an external reference plus a
samtools-style `--region`. It now also accepts plain **FASTA** or **FASTQ** (optionally
gzipped) with raw reads as input. The input format is detected by content sniffing, not
by file extension. For sequence inputs the `--region` argument is not required (and is
ignored with a warning, since raw reads carry no reference coordinates); reads go
straight to realignment against the graph path sequences. The BAM/CRAM path is
unchanged.

## Motivation

A FASTA/FASTQ file contains only raw sequences with names — no reference coordinates.
Coordinates in a region string exist only because a BAM/CRAM records where each read was
aligned; `samtools view <bam> <region>` is a cheap pre-filter to pull out the reads of
interest. With FASTA input the equivalent filter is simply preparing the FASTA with the
relevant sequences; from that point the pipeline is format-agnostic (realign → inject →
coverage → COSIGT ranking).

## Changes

### `src/commands/geno.rs`

- **`InputFormat`** enum: `Bam`, `Cram`, `Fasta`, `Fastq`.
- **`detect_input_format(path)`**: content sniffing.
  - Plain `BAM\x01` / `CRAM` magic → BAM/CRAM.
  - gzip/BGZF magic (`1f 8b`) → decompress, check payload for `BAM\x01` (real BAM files
    are BGZF-compressed, so the BAM magic is *inside* the gzip stream) or `CRAM`.
  - Otherwise plain text: first line `>` → FASTA, `@` → FASTQ; otherwise a clear error
    ("expected BAM, CRAM, FASTA or FASTQ").
  - gzipped FASTA/FASTQ are transparently handled in all branches.
- **`prepare_sequence_reads(input, format, out_fastq)`**: FASTA → FASTQ conversion in
  pure Rust (read name = first token of the `>` header; dummy `I` qualities of matching
  length; multi-line records supported; `.gz` input decompressed via `flate2`); FASTQ is
  passed through. Returns the read count.
- **`run_geno`**: branches on the detected format.
  - BAM/CRAM: unchanged pipeline (region extraction → `samtools sort -n` →
    `samtools fastq`); bails with a clear error if `--region` is missing.
  - FASTA/FASTQ: all reads are converted/passed to the FASTQ work file; no region
    extraction, no `--alignment-reference` needed.
  - `--region`, `--min-mapq`, `--exclude-flags` supplied with sequence input produce
    warnings (they are meaningless without alignment records).
- `GenoConfig.region` is now `Option<String>`; `reads_extracted` is the sequence record
  count for sequence inputs; `derive_sample_id` also strips `.fa/.fasta/.fna/.fq/.fastq`
  (± `.gz`) suffixes.
- **Panplexity parser fix retained**: the positional bare-value format written by
  panplexity (one 0/1 mask value or weight per line, node ID = line ordinal) is parsed
  correctly, alongside the ID+value format. This fix was already present in the working
  tree (never committed upstream); a dedicated regression test was added (see below).

### `src/main.rs`

- `--region` is optional; help text documents that it is required for BAM/CRAM and
  ignored for FASTA/FASTQ.
- `--alignment` help updated: "Input BAM/CRAM aligned to an external reference, or
  FASTA/FASTQ with raw reads (plain or gzipped)".

### `README.md`

- `geno` section retitled "External BAM/CRAM or FASTA/FASTQ Genotyping"; documents both
  paths and includes a FASTA example.

## Tests performed

### Unit tests (`cargo test`)

New tests in `src/commands/geno.rs`:

- `test_detect_and_convert_sequence_inputs` — format detection for: FASTA, FASTQ,
  gzipped FASTA, plain-magic "BAM", **gzip/BGZF-wrapped BAM** (real BAM shape; the
  uncompressed-magic case alone is not representative), and FASTA→FASTQ conversion
  output including multi-line records and dummy qualities, plus FASTQ passthrough.
- `test_parse_panplexity_positional_format` — regression protection for the panplexity
  positional mask/weights fix (bare `0/1` per line → ordinal node IDs; bare weights →
  ordinal IDs with values). The pre-existing ID+value test continues to pass.

Existing suite (29 lib tests incl. the BAM/CRAM end-to-end pipeline integration test,
plus integration test binaries) — **all pass**.

### End-to-end smoke tests (login node, `likegt` conda env)

Using the HLA-F test graph (`tests/data/hla-f.k51.gfa`) and 60 simulated 200 bp reads:

| Input | Result |
|---|---|
| `reads.fa` (FASTA) | 60 reads extracted, 60 realigned, genotype reported |
| `reads.fq` (FASTQ) | identical genotype result |
| `reads.fa.gz` (gzipped FASTA) | works; same read counts |
| FASTA + `--region` supplied | warning printed, region ignored |
| garbage text file | clean error: "expected BAM, CRAM, FASTA or FASTQ" |

### Real-data regression validation (SLURM, 4 sample/locus pairs)

Reference data: 223,128 result files previously produced by
`job/geno_panplexity_168.sh` with the patched binary (panplexity parser fix).
Validation script: `job/geno_panplexity_val_newbin.sh` (SLURM array, 4 tasks, 2 CPUs).
New binary: release build containing **both** the panplexity fix (from the patched
source tree) **and** this FASTA/FASTQ change.

Pairs chosen to include non-trivial true-genotype ranks; all inputs are real CRAMs from
`cosigt_validation/bams`, real SV graphs from `gfa_paf_sv`, and **positional-format**
panplexity mask/weights (e.g. DMBT1: 15,882 bare values, exercising exactly the fixed
parser path):

| Locus | Sample | Mode | true rank | r1 cosine | margin | nodes | Result |
|---|---|---|---|---|---|---|---|
| ABCA13 | HG00097 | mask / mask+weights | 1 | 0.9816344558579119 | 0 | 22691 | identical |
| DMBT1 | HG03704 | mask / mask+weights | 10 | 0.9910295486287686 | 0.000722776 | 14071 | identical |
| RUNX1 | HG03784 | mask / mask+weights | 2 | 0.9863319028428891 | 0 | 12088 | identical |
| RYR1 | NA18970 | mask / mask+weights | 2 | 0.9802654610644870 | 0.000340574 | 11197 | identical |

- All 8 summary rows (true rank, rank-1 cosine, true cosine, margin, nodes_used) are
  **bit-identical** to the stored results.
- All 8 full **top-50 genotype rankings are byte-identical** to the historical
  `.top.tsv` files.
- Conclusion: the new binary reproduces historical BAM/CRAM genotyping exactly; the
  FASTA/FASTQ change causes no behavioural change on the BAM/CRAM path.

Outputs: `/moosefs/raid5/lpignata/old_dups/geno_panplexity_val_newbin/`

### Real-data FASTA-input validation (SLURM)

Question: does the *new FASTA path* produce the same genotyping result as the CRAM
path on the same reads?

- Source: `cosigt_validation/bams/HG03784/HG03784_RUNX1.cram`, region
  `chr21:34624901-37242226` (from `svs_refined.hgsvcv3.bed`).
- Extracted the same region reads as FASTA (`samtools view -b -F 0x904 -T <ref>` →
  `sort -n` → `samtools fasta`): **53,736 real reads**.
- Ran the new binary with `-b reads.fa` (no `--region`, no `--alignment-reference`),
  graph `RUNX1.gfa`, paths `RUNX1.fa.gz`, positional mask+weights.

| | true rank | r1 cosine | true cosine | nodes used |
|---|---|---|---|---|
| From FASTA | 2 | 0.9863319028428891 | 0.9863319028428891 | 12088 |
| From CRAM | 2 | 0.9863319028428891 | 0.9863319028428891 | 12088 |

The full top-50 genotype ranking is **byte-identical** between the FASTA and CRAM runs
(true haplotypes `HG03784#1#JBIRDL010000008.1... + HG03784#2#JBIRDM010000002.1...` at
rank 2 in both).

Outputs: `/moosefs/raid5/lpignata/old_dups/geno_fasta_test/`

## Issues found along the way

1. **BGZF detection bug (caught during validation):** the first detection
   implementation checked for `BAM\x01` only at the start of the file; real BAM files
   are BGZF/gzip-compressed, so they start with the gzip magic. The external-BAM
   integration test (which runs only when samtools/minimap2/gfainject/gafpack are in
   PATH) caught this; detection now decompresses gzip files and checks the payload.
2. **Environment:** `gfainject` is no longer installed in `~/.cargo/bin` or
   `panplexity/target/release`; it currently resolves only via the Guix profile
   (`/gnu/store/3yj1pbygmz3zd31lgxqgbbxadihmbrfj-profile/bin`). The validation scripts
   append that profile at the end of PATH (conda-env tools keep precedence).
   `job/geno_panplexity_168.sh` may need the same addition before its next run.
3. **SLURM note:** `/tmp` on the login node is not visible to workers; test data for
   jobs must live on MooseFS.

## Notes for review

- The panplexity positional-format parser fix travels together with this change because
  it exists only as an uncommitted working-tree patch on the machine that builds the
  production binary. Committing it here protects it with tests and gives it upstream
  visibility; the regression validation above confirms current production results are
  reproduced exactly.
- `--region` is now optional; existing BAM/CRAM users must still supply it (clear error
  message if missing), so no silent behaviour change.
