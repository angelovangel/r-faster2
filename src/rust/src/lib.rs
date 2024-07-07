use extendr_api::prelude::*;
use kseq::parse_path;


/// FASTQ read lengths
/// 
/// Obtain read lengths from a fastq/fastq.gz file 
/// 
/// @param infile Path to fastq/fastq.gz file
/// @return Numeric vector with fastq records lengths
/// @export
#[extendr]
fn fq_lengths(infile: &str) -> Vec<i64> {
    let mut len_vector: Vec<i64> = Vec::new();
    //let mut bases: i64 = 0;
    let mut records = parse_path(infile).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        //bases += record.len() as i64;
        let len = record.len() as i64;
        len_vector.push(len);
    }
    len_vector
} 

// helper function for fq_quals
fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum = 0.0;
    for &item in q {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        qprob_sum += prob
    }
    qprob_sum
}

/// FASTQ 'mean' quality scores
/// 
/// Obtain'mean' fastq read quality values. Mean Phred scores are calculated by converting to probabiliy, calculating mean quality, then converting back to Phred.
/// 
/// @param infile Path to fastq/fastq.gz file
/// @param phred Logical, report Phred score (error probability otherwise)
/// @return Numeric vector with 'mean' quality score per read
/// @export
#[extendr]
fn fq_quals(infile: &str, phred: bool) -> Vec<f32> {
    let mut qual_vector: Vec<f32> = Vec::new();
    let mut records = parse_path(infile).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        let mean_prob = qscore_probs( record.qual().as_bytes() ) / record.len() as f32;
        if phred {
            qual_vector.push(-10.0 * mean_prob.log10()); // convert mean_prob to phred score
        } else {
            qual_vector.push(mean_prob)
        }
    }
    qual_vector
}

/// FASTQ number of reads
/// 
/// Obtain the number of reads in a fastq/fastq.gz file
/// 
/// @param infile Path to fastq/fastq.gz file
/// @return Numeric 
/// @export
#[extendr]
fn fq_reads(infile: &str) -> i64 {
    let mut reads: i64 = 0;
    let mut records = parse_path(infile).unwrap();
    while let Some(_) = records.iter_record().unwrap() {
        reads += 1;
    }
    reads
}

/// FASTQ number of bases
/// 
/// Obtain the number of bases in a fastq/fastq.gz file
/// 
/// @param infile Path to fastq/fastq.gz file
/// @return Numeric
/// @export
#[extendr]
fn fq_bases(infile: &str) -> i64 {
    let mut bases: i64 = 0;
    let mut records = parse_path(infile).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        bases += record.len() as i64;
    }
    bases
}


// helper functions for summary
fn get_n_bases(seq: &[u8]) -> i64 {
    let mut n = 0;
    for s in seq {
        if matches!(s, &78u8 | &110u8) { // N or n
            n += 1;
        }
    }
    n
}

fn get_gc_bases(seq: &[u8]) -> i64 {
    let mut n: i64 = 0;
    for s in seq {
        if matches!(s, &103u8 | &71u8 |  &99u8 | &67u8) { //G, g, C or c
            n += 1;
        }
    }
        n
}

fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
    let mut n = 0;
    // functional style is much slower?
    //q.iter().filter(|x| **x >= qx).collect::<Vec<&u8>>().len() as i64
    for item in q {
        if *item >= qx {
            n += 1
        }
    }
    n
}

fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {
    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.iter().sum::<i64>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers
        .iter()
        .scan(0, |sum, i| {
            *sum += i;
            Some(*sum)
        })
        .collect::<Vec<_>>();
    let n50_index = cumsum.iter().position(|&x| x > halfsum as i64).unwrap();

    numbers[n50_index]
}

/// Return vector with summary values for the fastq file, this is not exported in the R package namespace
#[extendr]
//filename, reads, bases, num_n, min_len, max_len, n50, gc_content, q20
fn rust_summary(infile: &str) -> Vec<i64> {
    let mut records = parse_path(infile).unwrap();
    
    let mut len_vector: Vec<i64> = Vec::new();
    let mut reads: i64 = 0;
    let mut bases: i64 = 0;
    let mut num_n: i64 = 0;
    let mut gc_bases: i64 = 0;
    let mut qual20: i64 = 0;
    let mut qual30: i64 = 0;

    while let Some(record) = records.iter_record().unwrap() {
        bases += record.len() as i64;
        reads += 1;
        num_n += get_n_bases(record.seq().as_bytes() );
        len_vector.push(record.len() as i64);
        // gc bases are written in the vector, content is calculated in R
        gc_bases += get_gc_bases(record.seq().as_bytes());
        qual20 += get_qual_bases(record.qual().as_bytes(), 53); // 33 offset
        qual30 += get_qual_bases(record.qual().as_bytes(), 63); // 33 offseq
    }
    let min_len = len_vector.iter().min().unwrap().clone();
    let max_len = len_vector.iter().max().unwrap().clone();
    let n50 = get_nx(&mut len_vector, 0.5);

    let summary = vec![
        reads, bases, num_n, min_len, max_len, n50, gc_bases, qual20, qual30
    ];
    summary
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rfaster2;
    fn fq_lengths;
    fn fq_quals;
    fn fq_reads;
    fn fq_bases;
    fn rust_summary;
}
