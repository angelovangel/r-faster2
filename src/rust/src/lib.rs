use extendr_api::prelude::*;
use kseq::parse_path;


/// Return vector of fastq read lengths
///@export
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

/// Return vector of 'mean' fastq read quality values
/// @export
//#[extendr]
//fn fq_qual(infile: &str) -> Vec<i64> {
//}

/// Return number of reads in a fastq file
/// @export
#[extendr]
fn fq_nreads(infile: &str) -> i64 {
    let mut reads: i64 = 0;
    let mut records = parse_path(infile).unwrap();
    while let Some(_) = records.iter_record().unwrap() {
        reads += 1;
    }
    reads
}

/// Return number of bases in a fastq file
/// @export
#[extendr]
fn fq_nbases(infile: &str) -> i64 {
    let mut bases: i64 = 0;
    let mut records = parse_path(infile).unwrap();
    while let Some(record) = records.iter_record().unwrap() {
        bases += record.len() as i64;
    }
    bases
}


// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod rfaster2;
    fn fq_lengths;
    fn fq_nreads;
    fn fq_nbases;
}
