use clap::{Parser, Subcommand};
use anyhow::Result;

#[derive(Parser)]
#[command(name = "likegt")]
#[command(about = "A tool for pangenome graph-based genotyping validation")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run validation tests (hold-0-out, hold-2-out, reference bias)
    Validate {
        /// Input FASTA file with sequences
        #[arg(short, long)]
        fasta: String,
        
        /// Output directory for results
        #[arg(short, long, default_value = "validation_results")]
        output: String,
        
        /// Number of threads to use
        #[arg(short, long, default_value = "4")]
        threads: usize,
        
        /// K-mer size for graph construction
        #[arg(short, long, default_value = "51")]
        kmer_size: usize,
    },
    
    /// Check if a graph is suitable for genotyping
    Check {
        /// Input GFA file
        #[arg(short, long)]
        gfa: String,
        
        /// Output report file
        #[arg(short, long, default_value = "genotyping_report.txt")]
        output: String,
    },
    
    /// Build a pangenome graph from FASTA sequences using allwave + seqwish + odgi
    Build {
        /// Input FASTA file (.fa, .fa.gz)
        #[arg(short, long)]
        fasta: String,
        
        /// Output GFA file (or prefix for multiple k-mer sizes)
        #[arg(short, long)]
        output: String,
        
        /// K-mer size for seqwish (or comma-separated list: 0,25,51,101)
        #[arg(short, long, default_value = "51")]
        kmer_sizes: String,
        
        /// Number of threads
        #[arg(short, long, default_value = "8")]
        threads: usize,
        
        /// Allwave pruning mode (none, low, medium, high)
        #[arg(long, default_value = "none")]
        pruning: String,
        
        /// Generate PNG visualizations with odgi viz
        #[arg(long)]
        visualize: bool,
        
        /// Keep intermediate files (PAF, seqwish GFA)
        #[arg(long)]
        keep_intermediates: bool,

        /// Graph builder backend (allwave-seqwish, impg)
        #[arg(long, default_value = "allwave-seqwish")]
        builder: String,

        /// impg GFA engine when using --builder impg (pggb, seqwish, poa, or partitioned forms like pggb:10000)
        #[arg(long, default_value = "pggb")]
        impg_gfa_engine: String,

        /// Run panplexity on each final GFA
        #[arg(long)]
        panplexity: bool,

        /// Panplexity window size
        #[arg(long, default_value = "100")]
        panplexity_window_size: usize,

        /// Panplexity threshold, either a number or "auto"
        #[arg(long, default_value = "auto")]
        panplexity_threshold: String,

        /// Panplexity IQR multiplier for auto thresholding
        #[arg(long, default_value = "1.5")]
        panplexity_iqr_multiplier: f64,

        /// Panplexity complexity metric (linguistic or entropy)
        #[arg(long, default_value = "linguistic")]
        panplexity_complexity: String,
    },
    
    /// Run complete hold-out validation pipeline
    HoldOut {
        /// Input FASTA sequences  
        #[arg(short, long)]
        fasta: String,
        
        /// Input graph file (.gfa or .og)
        #[arg(short, long)]
        graph: String,
        
        /// Output directory for results
        #[arg(short, long, default_value = "hold2out_results")]
        output: String,
        
        /// Test individual to hold out (e.g., "HG00096" or "all" for all samples)
        #[arg(short, long, default_value = "all")]
        individual: String,
        
        /// Number of haplotypes to hold out (1 or 2)
        #[arg(long, default_value = "2")]
        hold: usize,
        
        /// Ploidy (number of haplotypes per individual)
        #[arg(short, long, default_value = "2")]
        ploidy: usize,
        
        /// Number of threads
        #[arg(short, long, default_value = "4")]
        threads: usize,
        
        /// K-mer size
        #[arg(short, long, default_value = "51")]
        kmer_size: usize,
        
        /// Read simulator (wgsim, mason, pbsim3)
        #[arg(long, default_value = "wgsim")]
        simulator: String,
        
        /// Read length for simulation
        #[arg(long, default_value = "150")]
        read_length: usize,
        
        /// Coverage depth for read simulation
        #[arg(long, default_value = "30")]
        coverage_depth: usize,
        
        /// Fragment length mean (paired-end reads)
        #[arg(long, default_value = "500")]
        fragment_length: usize,
        
        /// Fragment length std dev
        #[arg(long, default_value = "50")]
        fragment_std: usize,
        
        /// Aligner (minimap2, bwa-mem, strobealign)
        #[arg(long, default_value = "minimap2")]
        aligner: String,
        
        /// Alignment preset (sr for short reads, map-ont for long reads)
        #[arg(long, default_value = "sr")]
        preset: String,
        
        /// Keep intermediate files
        #[arg(long)]
        keep_files: bool,
        
        /// Extract sequences matching prefix from pangenome FASTA for bias filtering
        /// e.g., "CHM13" or "grch38" - filters reads that align to sequences matching prefix  
        #[arg(long)]
        bias_prefix: Option<String>,
        
        /// External reference FASTA file for bias filtering (e.g., grch38_chr6_pansn.fa.gz)
        /// Uses complete reference FASTA file for initial alignment before filtering
        #[arg(long)]
        bias_fasta: Option<String>,
        
        /// Output sequence-level QV validation
        #[arg(long)]
        sequence_qv: bool,
        
        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
        
        /// Output format (text, json, csv, tsv, table)
        #[arg(long, default_value = "text")]
        format: String,
    },

    /// Genotype reads from an external alignment (BAM/CRAM) or from raw sequences (FASTA/FASTQ)
    Geno {
        /// Input BAM/CRAM aligned to an external reference, or FASTA/FASTQ with raw reads
        #[arg(short = 'b', long = "alignment", visible_alias = "bam", visible_alias = "cram")]
        alignment: String,

        /// Region to extract from the alignment, in samtools format (for example chr6:29600000-29720000); required for BAM/CRAM, ignored for FASTA/FASTQ
        #[arg(short, long)]
        region: Option<String>,

        /// Pangenome graph to genotype against
        #[arg(short, long)]
        graph: String,

        /// FASTA of graph paths to realign reads against; existing BWA indexes are reused
        #[arg(short = 'x', long = "index-sequence", visible_alias = "graph-fasta", visible_alias = "paths-fasta")]
        index_sequence: Option<String>,

        /// Reference FASTA for CRAM decoding or external-alignment region extraction
        #[arg(long = "alignment-reference", visible_alias = "cram-reference")]
        alignment_reference: Option<String>,

        /// Existing odgi/gafpack-compatible reference coverage TSV(.gz)
        #[arg(long)]
        reference_coverage: Option<String>,

        /// Output directory for genotyping results
        #[arg(short, long, default_value = "geno_results")]
        output: String,

        /// Sample ID to use in outputs
        #[arg(short, long)]
        sample: Option<String>,

        /// Ploidy for genotype combinations
        #[arg(short, long, default_value = "2")]
        ploidy: usize,

        /// Number of threads
        #[arg(short, long, default_value = "4")]
        threads: usize,

        /// Aligner for realigning region reads to graph path sequences (minimap2, bwa-mem)
        #[arg(long, default_value = "minimap2")]
        aligner: String,

        /// Aligner preset, for example sr for short reads or map-ont for ONT
        #[arg(long, default_value = "sr")]
        preset: String,

        /// Do not build a missing BWA index for --index-sequence
        #[arg(long)]
        no_build_index: bool,

        /// Minimum MAPQ when extracting reads from the external alignment
        #[arg(long, default_value = "0")]
        min_mapq: u32,

        /// SAM flags to exclude when extracting reads from the external alignment
        #[arg(long, default_value = "0x904")]
        exclude_flags: String,

        /// Auto-load sibling Panplexity mask/weights files for the graph
        #[arg(long)]
        panplexity: bool,

        /// Panplexity node mask/list to exclude from genotyping
        #[arg(long, visible_alias = "pamplexity-mask")]
        panplexity_mask: Option<String>,

        /// Panplexity node weights to use in weighted cosine similarity
        #[arg(long, visible_alias = "pamplexity-weights")]
        panplexity_weights: Option<String>,

        /// Exclude haplotypes whose names contain this pattern; may be repeated
        #[arg(long)]
        exclude_haplotype: Vec<String>,

        /// Number of top genotype combinations to write and print; 0 writes all
        #[arg(long, default_value = "10")]
        top: usize,

        /// Keep intermediate BAM/SAM/GAF/FASTQ files
        #[arg(long)]
        keep_files: bool,

        /// Output format (text, table, tsv, json)
        #[arg(long, default_value = "text")]
        format: String,

        /// Verbose progress output
        #[arg(short, long)]
        verbose: bool,
    },
    
    /// Compute maximum attainable QV using sequence alignment
    MaxQv {
        /// Input FASTA file with sequences
        #[arg(short, long)]
        fasta: String,
        
        /// Individual to analyze (name or "all")
        #[arg(short, long)]
        individual: Option<String>,
        
        /// Output file (optional, defaults to stdout)
        #[arg(short, long)]
        output: Option<String>,
        
        /// Number of threads for alignment
        #[arg(short, long, default_value = "4")]
        threads: usize,
        
        /// Alignment method: "allwave" (exact) or "wfmash" (approximate, faster)
        #[arg(short = 'm', long, default_value = "allwave")]
        method: String,
        
        /// Sparsification strategy for allwave (default: tree:5:0:0)
        /// Options: "none" for exact all-vs-all, "tree:N:F:R" for neighbor-based
        /// (Only used with allwave method)
        #[arg(short = 'p', long, default_value = "tree:5:0:0")]
        sparsification: String,
        
        /// Verbose output
        #[arg(short, long)]
        verbose: bool,
    },
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::init();
    
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Validate { fasta, output, threads, kmer_size } => {
            likegt::commands::validate::run_validation(
                &fasta,
                &output, 
                threads,
                kmer_size,
            ).await
        }
        
        Commands::Check { gfa, output } => {
            likegt::commands::check::check_graph_genotyping_suitability(&gfa, &output)
        }
        
        Commands::Build {
            fasta,
            output,
            kmer_sizes,
            threads,
            pruning,
            visualize,
            keep_intermediates,
            builder,
            impg_gfa_engine,
            panplexity,
            panplexity_window_size,
            panplexity_threshold,
            panplexity_iqr_multiplier,
            panplexity_complexity,
        } => {
            let panplexity_config = panplexity.then(|| {
                likegt::pipeline::build::PanplexityConfig {
                    window_size: panplexity_window_size,
                    threshold: panplexity_threshold,
                    iqr_multiplier: panplexity_iqr_multiplier,
                    complexity: panplexity_complexity,
                }
            });

            match builder.as_str() {
                "allwave-seqwish" | "allwave" | "seqwish" => {
                    likegt::pipeline::build::build_graph_allwave_seqwish(
                        &fasta,
                        &output,
                        &kmer_sizes,
                        threads,
                        &pruning,
                        visualize,
                        keep_intermediates,
                        panplexity_config.as_ref(),
                    ).await
                }
                "impg" => {
                    let output_gfa = if output.ends_with(".gfa") {
                        output
                    } else {
                        format!("{}.gfa", output)
                    };
                    likegt::pipeline::build::build_graph_with_impg(
                        &fasta,
                        &output_gfa,
                        threads,
                        &impg_gfa_engine,
                        panplexity_config.as_ref(),
                    ).await
                }
                other => Err(anyhow::anyhow!(
                    "Unsupported build backend '{}'. Use 'allwave-seqwish' or 'impg'.",
                    other
                )),
            }
        }
        
        Commands::HoldOut { 
            fasta,
            graph,
            output, 
            individual,
            hold,
            ploidy, 
            threads, 
            kmer_size,
            simulator,
            read_length,
            coverage_depth,
            fragment_length,
            fragment_std,
            aligner,
            preset,
            keep_files,
            bias_prefix,
            bias_fasta,
            sequence_qv,
            verbose,
            format 
        } => {
            // Check if batch processing
            if individual == "all" || individual.contains(',') {
                likegt::commands::hold2out::run_batch_hold2out(
                    &fasta,
                    &graph,
                    &output,
                    &individual,
                    hold,
                    ploidy,
                    threads,
                    kmer_size,
                    &simulator,
                    read_length,
                    coverage_depth,
                    fragment_length,
                    fragment_std,
                    &aligner,
                    &preset,
                    keep_files,
                    bias_prefix.as_deref(),
                    bias_fasta.as_deref(),
                    sequence_qv,
                    verbose,
                    &format,
                ).await
            } else {
                likegt::commands::hold2out::run_complete_hold2out_pipeline(
                    &fasta,
                    &graph,
                    &output,
                    &individual,
                    hold,
                    ploidy,
                    threads,
                    kmer_size,
                    &simulator,
                    read_length,
                    coverage_depth,
                    fragment_length,
                    fragment_std,
                    &aligner,
                    &preset,
                    keep_files,
                    bias_prefix.as_deref(),
                    bias_fasta.as_deref(),
                    sequence_qv,
                    verbose,
                    &format,
                ).await
            }
        }

        Commands::Geno {
            alignment,
            region,
            graph,
            index_sequence,
            alignment_reference,
            reference_coverage,
            output,
            sample,
            ploidy,
            threads,
            aligner,
            preset,
            no_build_index,
            min_mapq,
            exclude_flags,
            panplexity,
            panplexity_mask,
            panplexity_weights,
            exclude_haplotype,
            top,
            keep_files,
            format,
            verbose,
        } => {
            likegt::commands::geno::run_geno(likegt::commands::geno::GenoConfig {
                alignment,
                region,
                graph,
                index_sequence,
                alignment_reference,
                reference_coverage,
                output_dir: output,
                sample_id: sample,
                ploidy,
                threads,
                aligner,
                preset,
                no_build_index,
                min_mapq,
                exclude_flags,
                panplexity,
                panplexity_mask,
                panplexity_weights,
                exclude_haplotype_patterns: exclude_haplotype,
                top_n: top,
                keep_files,
                format,
                verbose,
            }).await.map(|_| ())
        }
        
        Commands::MaxQv { fasta, individual, output, threads, method, sparsification, verbose } => {
            likegt::commands::max_qv::run_max_qv_analysis(
                &fasta,
                individual.as_deref(),
                output.as_deref(),
                threads,
                &method,
                Some(&sparsification),
                verbose,
            ).await
        }
    }
}
