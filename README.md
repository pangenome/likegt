# LikeGT - Pangenome Graph-based Genotyping Toolkit

A high-performance toolkit for validating and benchmarking pangenome graph-based genotyping using the COSIGT (COSine similarity Genotyping) algorithm.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Important Notes](#important-notes)
- [Performance](#performance)
- [Troubleshooting](#troubleshooting)

## Overview

LikeGT implements comprehensive validation pipelines for pangenome genotyping, including:
- Hold-k-out validation (k=0, 1, 2, ..., n)
- Maximum attainable QV computation
- Graph construction and validation
- Reference bias analysis
- Coverage-based genotype calling using cosine similarity

## Features

### Core Functionality
- **Hold-out Validation**: Complete pipeline from read simulation to genotype calling
- **Max QV Analysis**: Compute theoretical upper bounds on genotyping accuracy
- **Graph Building**: Construct pangenome graphs using allwave + seqwish + odgi, or the current `impg graph` pipeline
- **Graph Annotation**: Optional Panplexity low-complexity annotation, masks, weights, and Bandage coloring
- **Graph Validation**: Check graph suitability for genotyping
- **Batch Processing**: Process multiple samples efficiently
- **Multiple Output Formats**: text, JSON, CSV, TSV, table

### Algorithms
- **COSIGT**: Cosine similarity-based genotyping from coverage vectors
- **Coverage-based QV**: Quality values computed from similarity scores
- **Sequence QV**: Optional alignment-based quality values using biWFA

## Installation

### Prerequisites

LikeGT requires several bioinformatics tools to be installed and available in your PATH.

**Core Dependencies:**
- `cargo` (>= 1.70) - Rust compiler and package manager
- `odgi` (>= 0.8) - Optimized dynamic genome/graph implementation
- `minimap2` (>= 2.24) - Fast sequence alignment program
- `wgsim` (>= 1.0) - Read simulator for next-generation sequencing
- `gfainject` - Tool to project SAM/BAM alignments onto pangenome graphs
- `gafpack` - Coverage calculator for GAF alignments
- `samtools` (>= 1.14) - SAM/BAM file manipulation
- `seqtk` (>= 1.3) - FASTA/FASTQ processing toolkit

**Additional Dependencies:**
- `allwave` - All-vs-all sequence alignment (required for `max-qv` and `build` commands)
- `seqwish` - Graph inducer from alignments (required for `build` command)
- `panplexity` - Low-complexity graph annotation (optional, for `build --panplexity`)
- `impg` - Implicit pangenome graph toolkit (optional, for `build --builder impg`)
- `wfmash` - Approximate sequence alignment (optional, for `max-qv --method wfmash`; installed by the `impg` source build)
- `bwa` (>= 0.7.17) - Alternative aligner (optional, for `--aligner bwa-mem`)
- `zcat`/`gzip` - For handling compressed files

**Python Dependencies (for analysis scripts):**
- Python 3.8+
- Standard library only (no external packages required)

### Quick Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/likegt.git
cd likegt

# Check and install dependencies
./install_dependencies.sh

# Build release version (recommended)
cargo build --release

# Run tests
cargo test

# Install to cargo bin directory
cargo install --path .
```

### Guix / Octopus Setup

On Octopus, use the repository environment wrapper before building or running
tests:

```bash
source ./env.sh
cargo build --locked
cargo test --locked
```

The `geno` command also needs external graph/alignment tools. Check the local
PATH with:

```bash
scripts/check-geno-tools.sh
```

The repo includes `.guix/channels.scm` and `.guix/manifest.scm` for the
Guix-resolvable part of that toolchain (`samtools`, `minimap2`, `bwa`, `odgi`,
and a Guix-packaged `gafpack`):

```bash
guix time-machine -C .guix/channels.scm -- shell -m .guix/manifest.scm
```

`gfainject` and the current `gafpack` can also be installed as pinned
Cargo-built binaries:

```bash
scripts/install-geno-cargo-tools.sh
export PATH="$PWD/target/tools/bin:$PATH"
```

To run the external BAM/graph smoke test locally:

```bash
scripts/run-geno-smoke.sh
```

The smoke test installs the pinned Cargo-built `gfainject` and `gafpack`
binaries into `target/tools/bin`, checks the required tools, and runs the
external BAM genotyping test. The production `geno` command supports both the
current `gafpack --gfa/--gaf` CLI and the older `gafpack --graph/--alignments`
CLI.

### Manual Installation

If you prefer to install dependencies manually:

1. **Install Rust** (if not already installed):
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

2. **Install bioinformatics tools**:

Using conda/mamba (recommended):
```bash
conda create -n likegt python=3.10
conda activate likegt
conda install -c bioconda -c conda-forge \
    minimap2 samtools seqtk bwa odgi seqwish wgsim
```

Using apt (Ubuntu/Debian):
```bash
sudo apt update
sudo apt install minimap2 samtools seqtk bwa
```

3. **Install tools from source** (required for some dependencies):
- odgi: https://github.com/pangenome/odgi
- gfainject: https://github.com/ekg/gfainject  
- gafpack: https://github.com/ekg/gafpack
- allwave: https://github.com/ekg/allwave
- impg: https://github.com/pangenome/impg
- panplexity: https://github.com/AndreaGuarracino/panplexity

4. **Build LikeGT**:
```bash
cargo build --release
```

## Quick Start

### 1. Basic Hold-2-out Validation

```bash
# Single individual
likegt hold-out -f sequences.fa.gz -g graph.gfa -i HG00096 -v

# Multiple individuals
likegt hold-out -f sequences.fa.gz -g graph.gfa -i "HG00096,HG00171,HG00268" -v

# All individuals
likegt hold-out -f sequences.fa.gz -g graph.gfa -i all --format tsv > results.tsv
```

### 2. Genotype an External BAM/CRAM Region or Raw FASTA/FASTQ Reads

```bash
likegt geno \
  --alignment sample.bam \
  --region chr6:29600000-29720000 \
  --graph hla.gfa \
  --index-sequence hla.graph_paths.fa \
  --sample HG00096 \
  --output geno_results

# Raw reads work too: --region is not needed (and ignored) for FASTA/FASTQ input
likegt geno -b reads.fa -g hla.gfa -x hla.graph_paths.fa -s HG00096 -o geno_results

# Optionally exclude Panplexity-masked nodes from COSIGT scoring
likegt geno -b sample.cram -r chr6:29600000-29720000 -g hla.gfa \
  --alignment-reference grch38.fa \
  --index-sequence hla.graph_paths.fa \
  --panplexity-mask hla.panplexity.mask
```

### 3. Compute Maximum Attainable QV

```bash
# Single individual - find best non-self match
likegt max-qv -f sequences.fa.gz -i HG00096 -v

# Multiple individuals
likegt max-qv -f sequences.fa.gz -i "HG00096,HG00171,HG00268" -v

# All individuals in the dataset (shows allwave progress with -v)
likegt max-qv -f sequences.fa.gz -i all -v

# Example verbose output with progress:
# 📊 Found 65 individuals to process
# 🧬 Running allwave for all-vs-all alignment...
# [1.0s] 436/1242 (35.1%) 435.9 alignments/sec ETA: 1.8s
# [2.3s] 609/1242 (49.0%) 263.7 alignments/sec ETA: 2.4s
# ✅ Got 1310 alignment records
# [1/65] Processing HG00096...
#   Best match: HG00268#1 + HG04036#2
#   Max attainable QV: 31.4

# Save results to file
likegt max-qv -f sequences.fa.gz -i all -o max_qv_results.tsv
```

### 4. Build Pangenome Graph

```bash
# Build graph with specific k-mer size
likegt build -f sequences.fa -o output_prefix -k 51 -t 8

# Build multiple k-mer graphs
likegt build -f sequences.fa -o output_prefix -k 25,51,101 -t 8

# Build with the current impg graph pipeline
likegt build -f sequences.fa -o output.gfa --builder impg -t 16

# Add Panplexity low-complexity annotations to generated GFAs
likegt build -f sequences.fa -o output_prefix -k 51 --panplexity
```

## Commands

### `hold-out` - Hold-out Validation Pipeline

The hold-out validation pipeline evaluates genotyping accuracy by simulating a realistic genotyping scenario where the true haplotypes are excluded from the reference panel.

**How the Pipeline Works:**

1. **Sequence Extraction**: Extracts the target individual's haplotype sequences from the input FASTA
2. **Read Simulation**: Simulates paired-end sequencing reads from the extracted sequences using:
   - `wgsim` for short reads (default: 150bp reads, 30x coverage)
   - Configurable error rates and fragment sizes
3. **Read Alignment**: Maps simulated reads to the pangenome graph:
   - First aligns to linear sequences using `minimap2` or `bwa-mem`
   - Produces SAM/BAM alignment file
4. **Graph Projection**: Projects linear alignments onto the graph:
   - Converts SAM to GAF format using `gfainject`
   - Maps linear coordinates to graph node space
5. **Coverage Calculation**: Computes per-node coverage using `gafpack`:
   - Counts read support for each graph node
   - Generates coverage matrix for the sample
6. **Genotype Calling**: Uses COSIGT algorithm to find best matching genotype:
   - Computes cosine similarity between sample and reference coverage vectors
   - Identifies the pair of reference haplotypes with highest similarity
7. **Quality Assessment**: Calculates QV from cosine similarity:
   - QV = -10 * log10(1 - similarity)
   - Reports whether correct genotype was called

**Options:**
- `-f, --fasta`: Input FASTA file with haplotype sequences
- `-g, --graph`: Pangenome graph (.gfa or .og format)
- `-i, --individual`: Individual(s) to test (name, comma-separated, or "all")
- `--hold`: Number of haplotypes to hold out (default: 2)
- `-t, --threads`: Number of threads (default: 4)
- `--coverage-depth`: Simulated read coverage (default: 30)
- `--read-length`: Read length for simulation (default: 150)
- `--aligner`: Aligner to use (minimap2, bwa-mem, GraphAligner)
- `--simulator`: Read simulator (wgsim, mason, pbsim3)
- `--format`: Output format (text, json, csv, tsv, table)
- `-v, --verbose`: Verbose output

### `geno` - External BAM/CRAM or FASTA/FASTQ Genotyping

Genotypes reads from an existing BAM/CRAM alignment or from raw FASTA/FASTQ sequences. For BAM/CRAM, the command extracts reads from a samtools-style region, converts them to FASTQ, realigns them to graph path sequences, projects the alignments into the graph with `gfainject`, calculates node coverage with `gafpack`, and ranks the closest haplotype combinations with COSIGT. For FASTA/FASTQ input, all reads are used directly (no region extraction; FASTA is converted to FASTQ with dummy qualities).

**Options:**
- `-b, --alignment`: Input BAM/CRAM aligned to an external reference, or FASTA/FASTQ with raw reads (plain or gzipped)
- `-r, --region`: Region to extract from a BAM/CRAM, such as `chr6:29600000-29720000`; required for BAM/CRAM, ignored for FASTA/FASTQ
- `-g, --graph`: Pangenome graph to genotype against
- `-x, --index-sequence`: FASTA of graph paths to align against; existing BWA indexes are reused
- `--alignment-reference`: Reference FASTA for CRAM decoding
- `--reference-coverage`: Existing graph path coverage TSV(.gz); otherwise generated with `odgi paths -H`
- `--aligner`: Realigner (`minimap2`, default; or `bwa-mem`)
- `--panplexity`: Auto-load sibling `*.panplexity.mask` or `*.panplexity.weights.txt` files
- `--panplexity-mask`: Exclude masked/low-complexity graph nodes from COSIGT scoring
- `--panplexity-weights`: Use per-node weights in weighted cosine similarity

### `max-qv` - Maximum Attainable QV using Sequence Alignment

Computes the best possible QV achievable when genotyping held-out individuals by finding their most similar non-self haplotypes using actual sequence alignment with allwave.

**Options:**
- `-f, --fasta`: Input FASTA file with haplotype sequences
- `-i, --individual`: Individual to analyze (name, comma-separated list, or "all")
- `-o, --output`: Output file (optional, defaults to stdout)
- `-t, --threads`: Number of threads for allwave (default: 4)
- `-v, --verbose`: Verbose output

### `build` - Build Pangenome Graph

Constructs pangenome graphs using the default allwave + seqwish + odgi pipeline or the current `impg graph` pipeline.

**Options:**
- `-f, --fasta`: Input FASTA file
- `-o, --output`: Output prefix for generated files
- `-k, --kmer-sizes`: K-mer sizes (comma-separated)
- `-t, --threads`: Number of threads
- `--keep-intermediates`: Keep intermediate files
- `--builder`: Graph backend (`allwave-seqwish`, default; or `impg`)
- `--impg-gfa-engine`: impg GFA engine (`pggb`, `seqwish`, `poa`, or partitioned forms like `pggb:10000`)
- `--panplexity`: Run Panplexity on each final GFA
- `--panplexity-window-size`: Window size for complexity calculation
- `--panplexity-threshold`: Numeric threshold or `auto`
- `--panplexity-iqr-multiplier`: IQR multiplier for `auto` thresholding
- `--panplexity-complexity`: Complexity metric (`linguistic` or `entropy`)

When `--panplexity` is enabled, LikeGT writes sibling files for each final GFA:
`*.panplexity.gfa`, `*.low-complexity.bed`, `*.panplexity.mask`,
`*.panplexity.weights.txt`, and `*.panplexity.bandage.csv`.

### `check` - Validate Graph

Checks if a graph is suitable for genotyping.

**Options:**
- `-g, --graph`: Graph file to check (.gfa or .og)
- `-v, --verbose`: Verbose output

### `validate` - Run Validation Tests

Runs comprehensive validation tests on coverage data.

**Options:**
- `-f, --fasta`: Input FASTA file
- `-o, --output`: Output directory
- `-t, --threads`: Number of threads
- `-k, --kmer-size`: K-mer size for analysis

## Important Notes

### How max-qv Works

The `max-qv` command answers the question: **"What's the best possible genotyping accuracy we could achieve for a held-out individual using only other individuals in the dataset?"**

**Process:**
1. Runs allwave for all-vs-all sequence alignment
2. For each individual, finds the best matching non-self haplotypes
3. Computes QV from alignment identity: `QV = -10 * log10(1 - identity)`
4. Reports the best achievable genotype (pair of haplotypes) and their QV scores

**Understanding the Output:**
- **QV 20** = 99% sequence identity
- **QV 30** = 99.9% sequence identity  
- **QV 40** = 99.99% sequence identity
- **Typical max QV**: 25-35 for human genome regions

**Example Output:**
```
individual  target_hap1  target_hap2  qv_hap1  qv_hap2  avg_qv  identity_hap1  identity_hap2
HG00096     HG00268#1    HG04036#2    29.1     33.8     31.4    0.9988         0.9996
```
This shows that HG00096's best non-self match would be HG00268's haplotype 1 + HG04036's haplotype 2, achieving an average QV of 31.4.

### Coverage-based vs Sequence-based QV

**hold-out command (Coverage-based QV):**
- QV computed from cosine similarity between graph node coverage vectors
- Measures how well genotyping via graph coverage matches the truth
- Fast and reflects actual genotyping algorithm performance

**max-qv command (Sequence-based QV):**
- QV computed from actual DNA sequence alignment using allwave
- Measures theoretical best possible accuracy based on sequence similarity
- Shows upper bound on genotyping performance

### Performance Considerations

- Use `--release` builds for production (10-100x faster)
- Coverage matrices can be large; ensure sufficient RAM
- Batch processing is more efficient than individual runs
- Default threads: 4 (adjust based on your system)

## Performance

### Validated Results on HLA-F Dataset

| K-value | Mean Similarity | Mean QV | Samples | Status |
|---------|----------------|---------|---------|---------|
| k=51    | 99.70%        | 26.1    | 46      | **Optimal** |
| k=25    | 93.97%        | 12.7    | 66      | Acceptable |
| k=101   | 89.46%        | 10.1    | 67      | Poor |

### Benchmarks

- **Hold-2-out validation**: ~1.2s per individual
- **Max QV computation**: <1s for 130 haplotypes
- **Graph building**: Varies with input size
- **Memory usage**: ~2GB for typical HLA dataset

## Troubleshooting

### Common Issues

1. **"Command not found" errors**
   - Ensure all required tools are installed and in PATH
   - Check `geno` tools with: `scripts/check-geno-tools.sh`
   - Check core tools manually with: `which minimap2 gfainject gafpack`
   - Check optional build tools with: `which impg panplexity wfmash`

2. **High QV values (60.0)**
   - Indicates perfect match - check if hold-out is working
   - Verify reference doesn't contain test individual

3. **Dimension mismatch errors**
   - Coverage matrices must have consistent dimensions
   - Check for metadata columns that need skipping

4. **Memory errors**
   - Large coverage matrices require significant RAM
   - Consider processing in batches or using a machine with more memory

### Debug Mode

Enable verbose output for detailed pipeline information:
```bash
likegt hold-out -f input.fa -g graph.gfa -i HG00096 --verbose
```

### File Formats

- **FASTA**: Gzipped or uncompressed, with haplotype IDs like "HG00096#1"
- **GFA**: Graph format from odgi/vg tools
- **Coverage Matrix**: TSV with haplotypes as rows, nodes as columns
- **Output**: JSON, CSV, TSV, or human-readable text

## Citation

If you use LikeGT in your research, please cite:

```
LikeGT: A toolkit for pangenome graph-based genotyping validation
[Your publication details here]
```

## License

[Specify your license here]

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Acknowledgments

This toolkit builds upon several excellent tools:
- odgi for graph manipulation
- minimap2 for alignment
- gafpack for coverage calculation
- biWFA for sequence alignment
