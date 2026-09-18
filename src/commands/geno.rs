use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use itertools::Itertools;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use crate::math;
use crate::validation;

#[derive(Debug, Clone)]
pub struct GenoConfig {
    pub alignment: String,
    pub region: String,
    pub graph: String,
    pub index_sequence: Option<String>,
    pub alignment_reference: Option<String>,
    pub reference_coverage: Option<String>,
    pub output_dir: String,
    pub sample_id: Option<String>,
    pub ploidy: usize,
    pub threads: usize,
    pub aligner: String,
    pub preset: String,
    pub no_build_index: bool,
    pub min_mapq: u32,
    pub exclude_flags: String,
    pub panplexity: bool,
    pub panplexity_mask: Option<String>,
    pub panplexity_weights: Option<String>,
    pub exclude_haplotype_patterns: Vec<String>,
    pub top_n: usize,
    pub keep_files: bool,
    pub format: String,
    pub verbose: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct GenoCombination {
    pub rank: usize,
    pub haplotypes: Vec<String>,
    pub cosine_similarity: f64,
    pub graph_qv: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct GenoOutputFiles {
    pub genotype_tsv: String,
    pub sorted_combinations_tsv_gz: String,
    pub summary_json: String,
    pub reference_coverage_tsv_gz: String,
    pub sample_coverage_tsv_gz: String,
    pub index_sequence_fasta: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct GenoRunResult {
    pub sample_id: String,
    pub alignment: String,
    pub region: String,
    pub graph: String,
    pub ploidy: usize,
    pub aligner: String,
    pub reads_extracted: usize,
    pub reads_written_to_fastq: usize,
    pub reads_realigned: usize,
    pub reference_haplotypes: usize,
    pub graph_nodes_common: usize,
    pub graph_nodes_used: usize,
    pub panplexity_nodes_excluded: usize,
    pub top_genotypes: Vec<GenoCombination>,
    pub output_files: GenoOutputFiles,
}

#[derive(Debug)]
struct WorkFiles {
    reference_coverage: PathBuf,
    index_sequence: PathBuf,
    region_bam: PathBuf,
    name_sorted_bam: PathBuf,
    reads_fastq: PathBuf,
    realigned_sam: PathBuf,
    realigned_bam: PathBuf,
    sample_gaf: PathBuf,
    sample_coverage: PathBuf,
    genotype_tsv: PathBuf,
    sorted_combinations: PathBuf,
    summary_json: PathBuf,
}

#[derive(Debug, Clone)]
struct CoverageTable {
    columns: Vec<String>,
    ids: Vec<String>,
    values: Vec<Vec<f64>>,
}

#[derive(Debug, Default, Clone)]
struct NodeFilters {
    excluded: HashSet<String>,
    weights: HashMap<String, f64>,
}

#[derive(Debug)]
struct GenotypeInputs {
    ids: Vec<String>,
    references: Vec<Vec<f64>>,
    sample: Vec<f64>,
    weights: Option<Vec<f64>>,
    common_nodes: usize,
    used_nodes: usize,
    panplexity_nodes_excluded: usize,
}

/// Genotype an externally aligned BAM/CRAM region by realigning reads to graph paths.
pub async fn run_geno(config: GenoConfig) -> Result<GenoRunResult> {
    validate_config(&config)?;

    fs::create_dir_all(&config.output_dir)
        .with_context(|| format!("Failed to create output directory {}", config.output_dir))?;

    let output_dir = PathBuf::from(&config.output_dir);
    let sample_id = config
        .sample_id
        .clone()
        .unwrap_or_else(|| derive_sample_id(&config.alignment));
    let files = make_work_files(&output_dir, &sample_id);

    if config.verbose {
        println!("Starting genotyping for sample {}", sample_id);
        println!("  Alignment: {}", config.alignment);
        println!("  Region: {}", config.region);
        println!("  Graph: {}", config.graph);
    }

    let reference_coverage = prepare_reference_coverage(&config, &files)?;
    let index_sequence = prepare_index_sequence(&config, &files)?;
    let (mask_path, weights_path) = resolve_panplexity_files(&config)?;

    extract_region_bam(&config, &files.region_bam)?;
    let reads_extracted = count_bam_records(&files.region_bam, false, config.threads)?;

    bam_to_fastq(
        &files.region_bam,
        &files.name_sorted_bam,
        &files.reads_fastq,
        config.threads,
    )?;
    let reads_written_to_fastq = count_fastq_reads(&files.reads_fastq)?;
    if reads_written_to_fastq == 0 {
        bail!(
            "No reads were extracted from {} in region {}",
            config.alignment,
            config.region
        );
    }

    realign_reads(
        &config,
        &index_sequence,
        &files.reads_fastq,
        &files.realigned_sam,
    )?;
    sam_to_bam(&files.realigned_sam, &files.realigned_bam, config.threads)?;
    let reads_realigned = count_bam_records(&files.realigned_bam, true, config.threads)?;

    inject_bam_to_graph(&config.graph, &files.realigned_bam, &files.sample_gaf)?;
    generate_sample_coverage(&config.graph, &files.sample_gaf, &files.sample_coverage)?;

    let filters = load_node_filters(mask_path.as_deref(), weights_path.as_deref())?;
    let inputs = build_genotype_inputs(
        &reference_coverage,
        &files.sample_coverage,
        &filters,
        &config.exclude_haplotype_patterns,
    )?;

    let combinations = evaluate_genotypes(&inputs, config.ploidy, config.threads)?;
    if combinations.is_empty() {
        bail!("No genotype combinations were generated");
    }

    let top_count = if config.top_n == 0 {
        combinations.len()
    } else {
        config.top_n.min(combinations.len())
    };
    let top_genotypes = combinations[..top_count].to_vec();

    write_genotype_tsv(
        &files.genotype_tsv,
        &sample_id,
        config.ploidy,
        &combinations[0],
    )?;
    write_sorted_combinations(
        &files.sorted_combinations,
        config.ploidy,
        &combinations,
        top_count,
    )?;

    let result = GenoRunResult {
        sample_id: sample_id.clone(),
        alignment: config.alignment.clone(),
        region: config.region.clone(),
        graph: config.graph.clone(),
        ploidy: config.ploidy,
        aligner: config.aligner.clone(),
        reads_extracted,
        reads_written_to_fastq,
        reads_realigned,
        reference_haplotypes: inputs.ids.len(),
        graph_nodes_common: inputs.common_nodes,
        graph_nodes_used: inputs.used_nodes,
        panplexity_nodes_excluded: inputs.panplexity_nodes_excluded,
        top_genotypes,
        output_files: GenoOutputFiles {
            genotype_tsv: files.genotype_tsv.to_string_lossy().to_string(),
            sorted_combinations_tsv_gz: files.sorted_combinations.to_string_lossy().to_string(),
            summary_json: files.summary_json.to_string_lossy().to_string(),
            reference_coverage_tsv_gz: reference_coverage.to_string_lossy().to_string(),
            sample_coverage_tsv_gz: files.sample_coverage.to_string_lossy().to_string(),
            index_sequence_fasta: index_sequence.to_string_lossy().to_string(),
        },
    };

    fs::write(&files.summary_json, serde_json::to_string_pretty(&result)?)
        .with_context(|| format!("Failed to write {}", files.summary_json.display()))?;

    if !config.keep_files {
        cleanup_intermediates(&files, &config)?;
    }

    print_result(&result, &config.format)?;
    Ok(result)
}

fn validate_config(config: &GenoConfig) -> Result<()> {
    if config.ploidy == 0 {
        bail!("--ploidy must be at least 1");
    }
    if !Path::new(&config.alignment).exists() {
        bail!("Alignment file not found: {}", config.alignment);
    }
    if !Path::new(&config.graph).exists() {
        bail!("Graph file not found: {}", config.graph);
    }
    if let Some(path) = &config.index_sequence {
        if !Path::new(path).exists() {
            bail!("Index sequence FASTA not found: {}", path);
        }
    }
    if let Some(path) = &config.alignment_reference {
        if !Path::new(path).exists() {
            bail!("Alignment reference FASTA not found: {}", path);
        }
    }
    if let Some(path) = &config.reference_coverage {
        if !Path::new(path).exists() {
            bail!("Reference coverage file not found: {}", path);
        }
    }
    match config.aligner.as_str() {
        "minimap2" | "bwa-mem" | "bwa" => {}
        other => bail!("Unsupported aligner '{}'. Use minimap2 or bwa-mem.", other),
    }

    require_tool("samtools")?;
    require_tool("gfainject")?;
    require_tool("gafpack")?;
    match config.aligner.as_str() {
        "minimap2" => {
            require_tool("minimap2")?;
        }
        "bwa-mem" | "bwa" => {
            require_tool("bwa")?;
        }
        _ => {}
    }
    if config.index_sequence.is_none() || config.reference_coverage.is_none() {
        require_tool("odgi")?;
    }
    Ok(())
}

fn make_work_files(output_dir: &Path, sample_id: &str) -> WorkFiles {
    WorkFiles {
        reference_coverage: output_dir.join("reference_coverage.tsv.gz"),
        index_sequence: output_dir.join("graph_paths.fa"),
        region_bam: output_dir.join(format!("{}.region.bam", sample_id)),
        name_sorted_bam: output_dir.join(format!("{}.region.name.bam", sample_id)),
        reads_fastq: output_dir.join(format!("{}.region.fq", sample_id)),
        realigned_sam: output_dir.join(format!("{}.graph.sam", sample_id)),
        realigned_bam: output_dir.join(format!("{}.graph.bam", sample_id)),
        sample_gaf: output_dir.join(format!("{}.graph.gaf", sample_id)),
        sample_coverage: output_dir.join(format!("{}.sample_coverage.tsv.gz", sample_id)),
        genotype_tsv: output_dir.join(format!("{}.genotype.tsv", sample_id)),
        sorted_combinations: output_dir.join(format!("{}.sorted_combos.tsv.gz", sample_id)),
        summary_json: output_dir.join(format!("{}.geno.summary.json", sample_id)),
    }
}

fn prepare_reference_coverage(config: &GenoConfig, files: &WorkFiles) -> Result<PathBuf> {
    if let Some(path) = &config.reference_coverage {
        if config.verbose {
            println!("Using reference coverage: {}", path);
        }
        return Ok(PathBuf::from(path));
    }

    if config.verbose {
        println!("Generating reference coverage from graph paths");
    }

    let odgi = require_tool("odgi")?;
    let mut command = Command::new(odgi);
    command
        .args(["paths", "-i", &config.graph, "-H"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let output = command
        .output()
        .context("Failed to run odgi paths for reference coverage")?;
    if !output.status.success() {
        bail!(
            "odgi paths failed while generating reference coverage: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }

    write_gzip_bytes(&files.reference_coverage, &output.stdout)?;
    Ok(files.reference_coverage.clone())
}

fn prepare_index_sequence(config: &GenoConfig, files: &WorkFiles) -> Result<PathBuf> {
    if let Some(path) = &config.index_sequence {
        if config.verbose {
            println!("Using graph path FASTA/index sequence: {}", path);
        }
        return Ok(PathBuf::from(path));
    }

    if config.verbose {
        println!("Extracting graph paths as FASTA for realignment");
    }

    let odgi = require_tool("odgi")?;
    let mut command = Command::new(odgi);
    command
        .args(["paths", "-i", &config.graph, "-f"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let output = command
        .output()
        .context("Failed to run odgi paths for graph path FASTA")?;
    if !output.status.success() {
        bail!(
            "odgi paths failed while extracting graph path FASTA: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    fs::write(&files.index_sequence, &output.stdout)
        .with_context(|| format!("Failed to write {}", files.index_sequence.display()))?;
    Ok(files.index_sequence.clone())
}

fn resolve_panplexity_files(config: &GenoConfig) -> Result<(Option<PathBuf>, Option<PathBuf>)> {
    let mut mask = config.panplexity_mask.as_ref().map(PathBuf::from);
    let mut weights = config.panplexity_weights.as_ref().map(PathBuf::from);

    if config.panplexity {
        let prefix = config
            .graph
            .strip_suffix(".gfa")
            .unwrap_or(config.graph.as_str());

        if mask.is_none() {
            let candidate = PathBuf::from(format!("{}.panplexity.mask", prefix));
            if candidate.exists() {
                mask = Some(candidate);
            }
        }
        if weights.is_none() {
            let candidate = PathBuf::from(format!("{}.panplexity.weights.txt", prefix));
            if candidate.exists() {
                weights = Some(candidate);
            }
        }
        if mask.is_none() && weights.is_none() {
            bail!(
                "--panplexity was requested, but no sibling .panplexity.mask or .panplexity.weights.txt file was found for {}",
                config.graph
            );
        }
    }

    for path in mask.iter().chain(weights.iter()) {
        if !path.exists() {
            bail!("Panplexity filter file not found: {}", path.display());
        }
    }

    Ok((mask, weights))
}

fn extract_region_bam(config: &GenoConfig, output_bam: &Path) -> Result<()> {
    if config.verbose {
        println!("Extracting BAM/CRAM region {}", config.region);
    }

    let samtools = require_tool("samtools")?;
    let mut command = Command::new(samtools);
    command
        .arg("view")
        .arg("-b")
        .arg("-@")
        .arg(config.threads.to_string())
        .arg("-F")
        .arg(&config.exclude_flags)
        .arg("-o")
        .arg(output_bam);

    if config.min_mapq > 0 {
        command.arg("-q").arg(config.min_mapq.to_string());
    }
    if let Some(reference) = &config.alignment_reference {
        command.arg("-T").arg(reference);
    }

    command.arg(&config.alignment).arg(&config.region);
    run_checked(
        &mut command,
        "samtools view failed while extracting the requested region",
    )
}

fn bam_to_fastq(
    region_bam: &Path,
    name_sorted_bam: &Path,
    reads_fastq: &Path,
    threads: usize,
) -> Result<()> {
    let samtools = require_tool("samtools")?;
    let mut sort = Command::new(&samtools);
    sort.args(["sort", "-n", "-@"])
        .arg(threads.to_string())
        .arg("-o")
        .arg(name_sorted_bam)
        .arg(region_bam);
    run_checked(&mut sort, "samtools sort -n failed before FASTQ extraction")?;

    let output_file = File::create(reads_fastq)
        .with_context(|| format!("Failed to create {}", reads_fastq.display()))?;
    let mut fastq = Command::new(&samtools);
    fastq
        .args(["fastq", "-@"])
        .arg(threads.to_string())
        .arg("-n")
        .arg(name_sorted_bam)
        .stdout(Stdio::from(output_file))
        .stderr(Stdio::piped());
    run_checked(
        &mut fastq,
        "samtools fastq failed while converting region reads",
    )
}

fn realign_reads(
    config: &GenoConfig,
    target: &Path,
    reads_fastq: &Path,
    output_sam: &Path,
) -> Result<()> {
    if config.verbose {
        println!(
            "Realigning {} reads to {} with {}",
            reads_fastq.display(),
            target.display(),
            config.aligner
        );
    }

    let output_file = File::create(output_sam)
        .with_context(|| format!("Failed to create {}", output_sam.display()))?;

    match config.aligner.as_str() {
        "minimap2" => {
            let minimap2 = require_tool("minimap2")?;
            let mut command = Command::new(minimap2);
            command
                .args(["-ax", &config.preset, "-t"])
                .arg(config.threads.to_string())
                .arg("--secondary=no")
                .arg(target)
                .arg(reads_fastq)
                .stdout(Stdio::from(output_file))
                .stderr(Stdio::piped());
            run_checked(
                &mut command,
                "minimap2 failed while realigning reads to graph paths",
            )
        }
        "bwa-mem" | "bwa" => {
            ensure_bwa_index(target, !config.no_build_index, config.verbose)?;
            let bwa = require_tool("bwa")?;
            let mut command = Command::new(bwa);
            command
                .args(["mem", "-t"])
                .arg(config.threads.to_string())
                .arg(target)
                .arg(reads_fastq)
                .stdout(Stdio::from(output_file))
                .stderr(Stdio::piped());
            run_checked(
                &mut command,
                "bwa mem failed while realigning reads to graph paths",
            )
        }
        other => bail!("Unsupported aligner '{}'", other),
    }
}

fn sam_to_bam(input_sam: &Path, output_bam: &Path, threads: usize) -> Result<()> {
    let samtools = require_tool("samtools")?;
    let mut command = Command::new(samtools);
    command
        .arg("view")
        .arg("-b")
        .arg("-@")
        .arg(threads.to_string())
        .arg("-o")
        .arg(output_bam)
        .arg(input_sam);
    run_checked(
        &mut command,
        "samtools view failed while converting graph SAM to BAM",
    )
}

fn inject_bam_to_graph(graph: &str, bam: &Path, output_gaf: &Path) -> Result<()> {
    let output_file = File::create(output_gaf)
        .with_context(|| format!("Failed to create {}", output_gaf.display()))?;
    let gfainject = require_tool("gfainject")?;
    let mut command = Command::new(gfainject);
    command
        .args(["--gfa", graph, "--bam"])
        .arg(bam)
        .stdout(Stdio::from(output_file))
        .stderr(Stdio::piped());
    run_checked(
        &mut command,
        "gfainject failed while projecting alignments into the graph",
    )
}

fn generate_sample_coverage(graph: &str, gaf: &Path, output_coverage: &Path) -> Result<()> {
    let gafpack = require_tool("gafpack")?;
    let mut command = Command::new(&gafpack);
    if gafpack_supports_current_cli(&gafpack) {
        command.args(["--gfa", graph, "--gaf"]);
    } else {
        command.args(["--graph", graph, "--alignments"]);
    }
    command
        .arg(gaf)
        .arg("--len-scale")
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());

    let output = command.output().context("Failed to run gafpack")?;
    if !output.status.success() {
        bail!(
            "gafpack failed while calculating sample coverage: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    write_gzip_bytes(output_coverage, &output.stdout)
}

fn gafpack_supports_current_cli(gafpack: &Path) -> bool {
    let Ok(output) = Command::new(gafpack).arg("--help").output() else {
        return false;
    };
    let help = String::from_utf8_lossy(&output.stdout);
    help.contains("--gfa <GFA>") && help.contains("--gaf <GAF>")
}

fn build_genotype_inputs(
    reference_path: &Path,
    sample_path: &Path,
    filters: &NodeFilters,
    exclude_haplotype_patterns: &[String],
) -> Result<GenotypeInputs> {
    let reference = read_coverage_table(reference_path)?;
    let sample = read_coverage_table(sample_path)?;
    if sample.values.is_empty() {
        bail!("No sample coverage row found in {}", sample_path.display());
    }

    let sample_lookup = build_column_lookup(&sample.columns);
    let mut node_pairs = Vec::new();
    let mut common_nodes = 0usize;
    let mut excluded_by_panplexity = 0usize;

    for (ref_idx, ref_col) in reference.columns.iter().enumerate() {
        if is_metadata_column(ref_col) {
            continue;
        }

        let Some(sample_idx) = find_matching_column(ref_col, &sample_lookup) else {
            continue;
        };

        common_nodes += 1;
        let normalized = normalize_node_id(ref_col);
        if filters.excluded.contains(&normalized) {
            excluded_by_panplexity += 1;
            continue;
        }
        node_pairs.push((ref_idx, sample_idx, normalized));
    }

    if common_nodes == 0 {
        bail!(
            "No shared graph-node columns between reference coverage {} and sample coverage {}",
            reference_path.display(),
            sample_path.display()
        );
    }
    if node_pairs.is_empty() {
        bail!("All {} common graph nodes were filtered out", common_nodes);
    }

    let keep_rows: Vec<usize> = reference
        .ids
        .iter()
        .enumerate()
        .filter(|(_, id)| {
            !exclude_haplotype_patterns
                .iter()
                .any(|pattern| id.contains(pattern))
        })
        .map(|(idx, _)| idx)
        .collect();

    if keep_rows.is_empty() {
        bail!("No reference haplotypes remain after haplotype filtering");
    }

    let ids = keep_rows
        .iter()
        .map(|&idx| reference.ids[idx].clone())
        .collect::<Vec<_>>();
    let references = keep_rows
        .iter()
        .map(|&idx| {
            node_pairs
                .iter()
                .map(|(ref_idx, _, _)| reference.values[idx][*ref_idx])
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let sample_vec = node_pairs
        .iter()
        .map(|(_, sample_idx, _)| sample.values[0][*sample_idx])
        .collect::<Vec<_>>();

    let weights = if filters.weights.is_empty() {
        None
    } else {
        Some(
            node_pairs
                .iter()
                .map(|(_, _, node)| filters.weights.get(node).copied().unwrap_or(1.0))
                .collect::<Vec<_>>(),
        )
    };

    Ok(GenotypeInputs {
        ids,
        references,
        sample: sample_vec,
        weights,
        common_nodes,
        used_nodes: node_pairs.len(),
        panplexity_nodes_excluded: excluded_by_panplexity,
    })
}

fn evaluate_genotypes(
    inputs: &GenotypeInputs,
    ploidy: usize,
    threads: usize,
) -> Result<Vec<GenoCombination>> {
    let indices = (0..inputs.ids.len()).collect::<Vec<_>>();
    let combinations = indices
        .iter()
        .copied()
        .combinations_with_replacement(ploidy)
        .collect::<Vec<_>>();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("Failed to build Rayon thread pool")?;

    let mut results = pool.install(|| {
        combinations
            .par_iter()
            .map(|combo| {
                let coverage_refs = combo
                    .iter()
                    .map(|&idx| inputs.references[idx].as_slice())
                    .collect::<Vec<_>>();
                let combined = math::sum_vectors(&coverage_refs);
                let similarity = match &inputs.weights {
                    Some(weights) => weighted_cosine_similarity(&combined, &inputs.sample, weights),
                    None => math::cosine_similarity(&combined, &inputs.sample),
                };
                let graph_qv = validation::calculate_qv(1.0, similarity);
                let haplotypes = combo
                    .iter()
                    .map(|&idx| inputs.ids[idx].clone())
                    .collect::<Vec<_>>();
                GenoCombination {
                    rank: 0,
                    haplotypes,
                    cosine_similarity: similarity,
                    graph_qv,
                }
            })
            .collect::<Vec<_>>()
    });

    results.sort_by(|a, b| {
        b.cosine_similarity
            .partial_cmp(&a.cosine_similarity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (idx, result) in results.iter_mut().enumerate() {
        result.rank = idx + 1;
    }
    Ok(results)
}

fn read_coverage_table(path: &Path) -> Result<CoverageTable> {
    let file = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let reader: Box<dyn Read> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut lines = BufReader::new(reader).lines();

    let header = loop {
        let Some(line) = lines.next() else {
            bail!("Coverage file is empty: {}", path.display());
        };
        let line = line?;
        if !line.trim().is_empty() {
            break line;
        }
    };

    let header_fields = header
        .split('\t')
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    if header_fields.len() < 2 {
        bail!(
            "Coverage header has fewer than 2 columns in {}",
            path.display()
        );
    }
    let columns = header_fields[1..].to_vec();
    let mut ids = Vec::new();
    let mut values = Vec::new();

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != columns.len() + 1 {
            bail!(
                "Coverage row has {} fields, expected {} in {}",
                fields.len(),
                columns.len() + 1,
                path.display()
            );
        }
        ids.push(fields[0].to_string());
        let row = fields[1..]
            .iter()
            .map(|value| {
                value
                    .parse::<f64>()
                    .with_context(|| format!("Failed to parse coverage value '{}'", value))
            })
            .collect::<Result<Vec<_>>>()?;
        values.push(row);
    }

    Ok(CoverageTable {
        columns,
        ids,
        values,
    })
}

fn build_column_lookup(columns: &[String]) -> HashMap<String, usize> {
    let mut lookup = HashMap::new();
    for (idx, col) in columns.iter().enumerate() {
        lookup.insert(col.clone(), idx);
        lookup.insert(normalize_node_id(col), idx);
    }
    lookup
}

fn find_matching_column(column: &str, lookup: &HashMap<String, usize>) -> Option<usize> {
    lookup
        .get(column)
        .copied()
        .or_else(|| lookup.get(&normalize_node_id(column)).copied())
}

fn is_metadata_column(column: &str) -> bool {
    matches!(
        column,
        "group.name" | "path.length" | "path.step.count" | "sample" | "#sample"
    )
}

fn load_node_filters(mask: Option<&Path>, weights: Option<&Path>) -> Result<NodeFilters> {
    let mut filters = NodeFilters::default();
    if let Some(path) = mask {
        filters.excluded = parse_panplexity_mask(path)?;
    }
    if let Some(path) = weights {
        filters.weights = parse_panplexity_weights(path)?;
    }
    Ok(filters)
}

fn parse_panplexity_mask(path: &Path) -> Result<HashSet<String>> {
    let mut excluded = HashSet::new();
    let mut line_no = 0usize;
    for line in read_text_lines(path)? {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = split_filter_line(trimmed);
        if fields.is_empty() || looks_like_header(&fields) {
            continue;
        }

        line_no += 1;

        // panplexity writes one bare 0/1 per line, in node order, with NO node-ID
        // column (0 = low-complexity/exclude, 1 = normal/keep). A single numeric
        // field is that positional format -- the node ID is the line's ordinal
        // position, not the literal value on the line.
        if fields.len() == 1 {
            if let Ok(v) = fields[0].parse::<f64>() {
                if v == 0.0 {
                    excluded.insert(normalize_node_id(&line_no.to_string()));
                }
                continue;
            }
            excluded.insert(normalize_node_id(fields[0]));
            continue;
        }

        if let Some(mask_value) = fields
            .iter()
            .rev()
            .find_map(|field| field.parse::<f64>().ok())
        {
            if mask_value > 0.0 {
                excluded.insert(normalize_node_id(fields[0]));
            }
        } else {
            excluded.insert(normalize_node_id(fields[0]));
        }
    }
    Ok(excluded)
}

fn parse_panplexity_weights(path: &Path) -> Result<HashMap<String, f64>> {
    let mut weights = HashMap::new();
    let mut line_no = 0usize;
    for line in read_text_lines(path)? {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = split_filter_line(trimmed);
        if fields.is_empty() || looks_like_header(&fields) {
            continue;
        }
        line_no += 1;

        if fields.len() == 1 {
            if let Ok(w) = fields[0].parse::<f64>() {
                weights.insert(normalize_node_id(&line_no.to_string()), w.max(0.0));
            }
            continue;
        }

        if let Some(weight) = fields
            .iter()
            .rev()
            .find_map(|field| field.parse::<f64>().ok())
        {
            weights.insert(normalize_node_id(fields[0]), weight.max(0.0));
        }
    }
    Ok(weights)
}

fn split_filter_line(line: &str) -> Vec<&str> {
    line.split(|c: char| c == '\t' || c == ',' || c.is_whitespace())
        .filter(|field| !field.is_empty())
        .collect()
}

fn looks_like_header(fields: &[&str]) -> bool {
    fields.iter().any(|field| {
        let lower = field.to_ascii_lowercase();
        lower == "node" || lower == "node_id" || lower == "mask" || lower == "weight"
    })
}

fn read_text_lines(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let reader: Box<dyn Read> = if path.extension().and_then(|s| s.to_str()) == Some("gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    BufReader::new(reader)
        .lines()
        .collect::<std::io::Result<Vec<_>>>()
        .with_context(|| format!("Failed to read {}", path.display()))
}

fn normalize_node_id(raw: &str) -> String {
    raw.trim()
        .trim_start_matches("node.")
        .trim_start_matches('>')
        .trim_start_matches('<')
        .trim_end_matches('+')
        .trim_end_matches('-')
        .to_string()
}

fn weighted_cosine_similarity(a: &[f64], b: &[f64], weights: &[f64]) -> f64 {
    let mut dot = 0.0;
    let mut norm_a = 0.0;
    let mut norm_b = 0.0;

    for ((x, y), weight) in a.iter().zip(b.iter()).zip(weights.iter()) {
        if *weight <= 0.0 {
            continue;
        }
        dot += weight * x * y;
        norm_a += weight * x * x;
        norm_b += weight * y * y;
    }

    if norm_a == 0.0 || norm_b == 0.0 {
        0.0
    } else {
        dot / (norm_a.sqrt() * norm_b.sqrt())
    }
}

fn write_genotype_tsv(
    path: &Path,
    sample_id: &str,
    ploidy: usize,
    best: &GenoCombination,
) -> Result<()> {
    let mut file =
        File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let mut header = vec!["sample.id".to_string(), "rank".to_string()];
    for idx in 1..=ploidy {
        header.push(format!("haplotype.{}", idx));
    }
    header.push("cosine.similarity".to_string());
    header.push("graph.qv".to_string());
    writeln!(file, "{}", header.join("\t"))?;

    let mut row = vec![sample_id.to_string(), best.rank.to_string()];
    for idx in 0..ploidy {
        row.push(
            best.haplotypes
                .get(idx)
                .cloned()
                .unwrap_or_else(|| "".to_string()),
        );
    }
    row.push(format!("{:.16}", best.cosine_similarity));
    row.push(format!("{:.6}", best.graph_qv));
    writeln!(file, "{}", row.join("\t"))?;
    Ok(())
}

fn write_sorted_combinations(
    path: &Path,
    ploidy: usize,
    combinations: &[GenoCombination],
    count: usize,
) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());

    let mut header = vec!["rank".to_string()];
    for idx in 1..=ploidy {
        header.push(format!("haplotype.{}", idx));
    }
    header.push("cosine.similarity".to_string());
    header.push("graph.qv".to_string());
    writeln!(encoder, "{}", header.join("\t"))?;

    for genotype in combinations.iter().take(count) {
        let mut row = vec![genotype.rank.to_string()];
        for idx in 0..ploidy {
            row.push(
                genotype
                    .haplotypes
                    .get(idx)
                    .cloned()
                    .unwrap_or_else(|| "".to_string()),
            );
        }
        row.push(format!("{:.16}", genotype.cosine_similarity));
        row.push(format!("{:.6}", genotype.graph_qv));
        writeln!(encoder, "{}", row.join("\t"))?;
    }

    encoder.finish()?;
    Ok(())
}

fn print_result(result: &GenoRunResult, format: &str) -> Result<()> {
    match format {
        "json" => {
            println!("{}", serde_json::to_string_pretty(result)?);
        }
        "tsv" => {
            println!("sample\trank\tgenotype\tcosine_similarity\tgraph_qv\treads_extracted\treads_realigned\tnodes_used");
            for genotype in &result.top_genotypes {
                println!(
                    "{}\t{}\t{}\t{:.8}\t{:.4}\t{}\t{}\t{}",
                    result.sample_id,
                    genotype.rank,
                    genotype.haplotypes.join("+"),
                    genotype.cosine_similarity,
                    genotype.graph_qv,
                    result.reads_extracted,
                    result.reads_realigned,
                    result.graph_nodes_used
                );
            }
        }
        "text" | "table" => {
            let best = &result.top_genotypes[0];
            println!("GENOTYPE RESULT");
            println!("sample: {}", result.sample_id);
            println!("region: {}", result.region);
            println!("best: {}", best.haplotypes.join(" + "));
            println!("similarity: {:.6}", best.cosine_similarity);
            println!("graph_qv: {:.2}", best.graph_qv);
            println!(
                "reads: {} extracted, {} realigned",
                result.reads_extracted, result.reads_realigned
            );
            println!(
                "nodes: {} used of {} common ({} Panplexity-excluded)",
                result.graph_nodes_used,
                result.graph_nodes_common,
                result.panplexity_nodes_excluded
            );
            println!("top genotypes:");
            for genotype in &result.top_genotypes {
                println!(
                    "  {}\t{:.6}\t{}",
                    genotype.rank,
                    genotype.cosine_similarity,
                    genotype.haplotypes.join(" + ")
                );
            }
        }
        other => bail!(
            "Unsupported output format '{}'. Use text, table, tsv, or json.",
            other
        ),
    }
    Ok(())
}

fn cleanup_intermediates(files: &WorkFiles, config: &GenoConfig) -> Result<()> {
    let mut paths = vec![
        files.region_bam.clone(),
        files.name_sorted_bam.clone(),
        files.reads_fastq.clone(),
        files.realigned_sam.clone(),
        files.realigned_bam.clone(),
        files.sample_gaf.clone(),
    ];

    if config.reference_coverage.is_some() {
        paths.retain(|path| path != &files.reference_coverage);
    }
    if config.index_sequence.is_some() {
        paths.retain(|path| path != &files.index_sequence);
    }

    for path in paths {
        if path.exists() {
            fs::remove_file(&path).with_context(|| {
                format!("Failed to remove intermediate file {}", path.display())
            })?;
        }
    }
    Ok(())
}

fn count_bam_records(path: &Path, mapped_only: bool, threads: usize) -> Result<usize> {
    let samtools = require_tool("samtools")?;
    let mut command = Command::new(samtools);
    command
        .arg("view")
        .arg("-c")
        .arg("-@")
        .arg(threads.to_string());
    if mapped_only {
        command.args(["-F", "4"]);
    }
    command.arg(path);

    let output = command.output().context("Failed to run samtools view -c")?;
    if !output.status.success() {
        bail!(
            "samtools view -c failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let count = String::from_utf8(output.stdout)?
        .trim()
        .parse::<usize>()
        .context("Failed to parse samtools view -c output")?;
    Ok(count)
}

fn count_fastq_reads(path: &Path) -> Result<usize> {
    let file = File::open(path).with_context(|| format!("Failed to open {}", path.display()))?;
    let line_count = BufReader::new(file).lines().count();
    Ok(line_count / 4)
}

fn ensure_bwa_index(target: &Path, build_if_missing: bool, verbose: bool) -> Result<()> {
    let required = ["amb", "ann", "bwt", "pac", "sa"];
    let missing = required
        .iter()
        .any(|suffix| !PathBuf::from(format!("{}.{}", target.to_string_lossy(), suffix)).exists());

    if !missing {
        return Ok(());
    }

    if !build_if_missing {
        bail!(
            "BWA index is missing for {} and --no-build-index was supplied",
            target.display()
        );
    }

    if verbose {
        println!(
            "BWA index missing for {}; building it once",
            target.display()
        );
    }

    let bwa = require_tool("bwa")?;
    let mut command = Command::new(bwa);
    command.arg("index").arg(target);
    run_checked(&mut command, "bwa index failed")
}

fn write_gzip_bytes(path: &Path, bytes: &[u8]) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(bytes)?;
    encoder.finish()?;
    Ok(())
}

fn run_checked(command: &mut Command, context: &str) -> Result<()> {
    let output = command.output().with_context(|| context.to_string())?;
    if !output.status.success() {
        bail!("{}: {}", context, String::from_utf8_lossy(&output.stderr));
    }
    Ok(())
}

fn require_tool(name: &str) -> Result<PathBuf> {
    which::which(name).with_context(|| format!("Required tool '{}' was not found in PATH", name))
}

fn derive_sample_id(alignment: &str) -> String {
    Path::new(alignment)
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("sample")
        .trim_end_matches(".bam")
        .trim_end_matches(".cram")
        .trim_end_matches(".sam")
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_parse_panplexity_mask_and_weights() {
        let dir = TempDir::new().unwrap();
        let mask = dir.path().join("test.mask");
        let weights = dir.path().join("test.weights.txt");

        fs::write(&mask, "node\tmask\nnode.2\t1\n3\t1\nnode.4\t0\nnode.5\n").unwrap();
        fs::write(&weights, "node\tweight\nnode.1\t1.0\n2\t0.25\nnode.3\t0\n").unwrap();

        let excluded = parse_panplexity_mask(&mask).unwrap();
        assert!(excluded.contains("2"));
        assert!(excluded.contains("3"));
        assert!(excluded.contains("5"));
        assert!(!excluded.contains("4"));

        let parsed_weights = parse_panplexity_weights(&weights).unwrap();
        assert_eq!(parsed_weights.get("1"), Some(&1.0));
        assert_eq!(parsed_weights.get("2"), Some(&0.25));
        assert_eq!(parsed_weights.get("3"), Some(&0.0));
    }

    #[test]
    fn test_parse_panplexity_positional_format() {
        // panplexity writes one bare value per line, in node order, with no ID
        // column: the line's ordinal position is the node ID.
        let dir = TempDir::new().unwrap();
        let mask = dir.path().join("positional.mask");
        let weights = dir.path().join("positional.weights.txt");

        fs::write(&mask, "1\n1\n0\n1\n0\n").unwrap();
        fs::write(&weights, "1.0\n0.5\n0.0\n0.25\n").unwrap();

        let excluded = parse_panplexity_mask(&mask).unwrap();
        assert_eq!(excluded.len(), 2);
        assert!(excluded.contains("3"));
        assert!(excluded.contains("5"));
        assert!(!excluded.contains("1"));

        let parsed_weights = parse_panplexity_weights(&weights).unwrap();
        assert_eq!(parsed_weights.get("1"), Some(&1.0));
        assert_eq!(parsed_weights.get("2"), Some(&0.5));
        assert_eq!(parsed_weights.get("3"), Some(&0.0));
        assert_eq!(parsed_weights.get("4"), Some(&0.25));
    }

    #[test]
    fn test_genotyping_with_panplexity_mask() {
        let dir = TempDir::new().unwrap();
        let reference = dir.path().join("reference.tsv.gz");
        let sample = dir.path().join("sample.tsv.gz");
        let mask = dir.path().join("mask.txt");

        write_gzip_bytes(
            &reference,
            b"path.name\tpath.length\tpath.step.count\tnode.1\tnode.2\tnode.3\nhapA\t3\t3\t1\t100\t0\nhapB\t3\t3\t0\t100\t1\nhapC\t3\t3\t10\t0\t10\n",
        )
        .unwrap();
        write_gzip_bytes(
            &sample,
            b"#sample\tnode.1\tnode.2\tnode.3\nsample\t1\t500\t1\n",
        )
        .unwrap();
        fs::write(&mask, "node.2\n").unwrap();

        let filters = load_node_filters(Some(&mask), None).unwrap();
        let inputs = build_genotype_inputs(&reference, &sample, &filters, &[]).unwrap();
        assert_eq!(inputs.common_nodes, 3);
        assert_eq!(inputs.used_nodes, 2);
        assert_eq!(inputs.panplexity_nodes_excluded, 1);

        let results = evaluate_genotypes(&inputs, 2, 2).unwrap();
        assert_eq!(results[0].haplotypes, vec!["hapA", "hapB"]);
        assert!((results[0].cosine_similarity - 1.0).abs() < 1e-10);
    }

    #[tokio::test]
    async fn test_geno_pipeline_with_external_bam_if_tools_available() {
        for tool in ["samtools", "minimap2", "gfainject", "gafpack"] {
            if which::which(tool).is_err() {
                eprintln!(
                    "Skipping external BAM geno test because {} is not installed",
                    tool
                );
                return;
            }
        }

        let dir = TempDir::new().unwrap();
        let seq_a = make_sequence(11, 240);
        let seq_b = make_sequence(29, 240);
        let graph = dir.path().join("test.gfa");
        let graph_paths = dir.path().join("graph_paths.fa");
        let reference_coverage = dir.path().join("reference_coverage.tsv.gz");
        let external_ref = dir.path().join("external.fa");
        let reads = dir.path().join("reads.fq");
        let external_sam = dir.path().join("external.sam");
        let external_bam = dir.path().join("external.bam");
        let sorted_bam = dir.path().join("external.sorted.bam");
        let output_dir = dir.path().join("geno");

        fs::write(
            &graph,
            format!(
                "H\tVN:Z:1.0\nS\t1\t{}\nS\t2\t{}\nP\thapA\t1+\t*\nP\thapB\t2+\t*\n",
                seq_a, seq_b
            ),
        )
        .unwrap();
        fs::write(
            &graph_paths,
            format!(">hapA\n{}\n>hapB\n{}\n", seq_a, seq_b),
        )
        .unwrap();
        write_gzip_bytes(
            &reference_coverage,
            b"path.name\tpath.length\tpath.step.count\tnode.1\tnode.2\nhapA\t240\t1\t240\t0\nhapB\t240\t1\t0\t240\n",
        )
        .unwrap();
        fs::write(&external_ref, format!(">chr1\n{}\n", seq_b)).unwrap();

        let mut fastq = String::new();
        for idx in 0..6 {
            let start = idx * 20;
            let read = &seq_b[start..start + 100];
            fastq.push_str(&format!(
                "@r{}\n{}\n+\n{}\n",
                idx,
                read,
                "I".repeat(read.len())
            ));
        }
        fs::write(&reads, fastq).unwrap();

        run_stdout_file(
            Command::new(require_tool("minimap2").unwrap())
                .args(["-ax", "sr"])
                .arg(&external_ref)
                .arg(&reads),
            &external_sam,
        )
        .unwrap();
        run_ok(
            Command::new(require_tool("samtools").unwrap())
                .args(["view", "-b", "-o"])
                .arg(&external_bam)
                .arg(&external_sam),
        )
        .unwrap();
        run_ok(
            Command::new(require_tool("samtools").unwrap())
                .args(["sort", "-o"])
                .arg(&sorted_bam)
                .arg(&external_bam),
        )
        .unwrap();
        run_ok(
            Command::new(require_tool("samtools").unwrap())
                .arg("index")
                .arg(&sorted_bam),
        )
        .unwrap();

        let result = run_geno(GenoConfig {
            alignment: sorted_bam.to_string_lossy().to_string(),
            region: "chr1:1-240".to_string(),
            graph: graph.to_string_lossy().to_string(),
            index_sequence: Some(graph_paths.to_string_lossy().to_string()),
            alignment_reference: None,
            reference_coverage: Some(reference_coverage.to_string_lossy().to_string()),
            output_dir: output_dir.to_string_lossy().to_string(),
            sample_id: Some("sampleB".to_string()),
            ploidy: 2,
            threads: 2,
            aligner: "minimap2".to_string(),
            preset: "sr".to_string(),
            no_build_index: false,
            min_mapq: 0,
            exclude_flags: "0x904".to_string(),
            panplexity: false,
            panplexity_mask: None,
            panplexity_weights: None,
            exclude_haplotype_patterns: Vec::new(),
            top_n: 3,
            keep_files: true,
            format: "json".to_string(),
            verbose: false,
        })
        .await
        .unwrap();

        assert!(result.reads_extracted > 0);
        assert!(result.reads_realigned > 0);
        assert_eq!(
            result.top_genotypes[0].haplotypes,
            vec!["hapB".to_string(), "hapB".to_string()]
        );
        assert!(result.top_genotypes[0].cosine_similarity > 0.99);
        assert!(Path::new(&result.output_files.genotype_tsv).exists());
    }

    fn make_sequence(seed: u32, len: usize) -> String {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut x = seed;
        let mut seq = Vec::with_capacity(len);
        for _ in 0..len {
            x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            seq.push(bases[((x >> 30) & 3) as usize]);
        }
        String::from_utf8(seq).unwrap()
    }

    fn run_ok(command: &mut Command) -> Result<()> {
        let output = command.output()?;
        if !output.status.success() {
            bail!("{}", String::from_utf8_lossy(&output.stderr));
        }
        Ok(())
    }

    fn run_stdout_file(command: &mut Command, output_path: &Path) -> Result<()> {
        let output = command.output()?;
        if !output.status.success() {
            bail!("{}", String::from_utf8_lossy(&output.stderr));
        }
        fs::write(output_path, output.stdout)?;
        Ok(())
    }
}
