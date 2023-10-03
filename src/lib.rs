use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;
use std::fs::File;
use std::cmp::{max, min};
use std::time::Instant;
use std::error::Error;

use natord::compare;

use chrono::Datelike;

use colored::Colorize;
use peak_alloc::PeakAlloc;

use log::Level;

use indoc::indoc;

mod bed;
use bed::BedRecord;

mod codon;
use codon::Codon;

mod error;
use error::ParseError;


const SOURCE: &str = "bed2gff";
const VERSION: &str = "0.1.0";
const GFF3: &str = "##gff-version 3";
const PROVIDER: &str = "bed2gff";
const REPOSITORY: &str = "github.com/alejandrogzi/bed2gff";


#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;



fn get_isoforms(path: PathBuf) -> Result<HashMap<String, String>, ParseError> {
    let file: File = File::open(path).unwrap();
    let reader: BufReader<File> = BufReader::new(file);
    let mut isoforms: HashMap<String, String> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let content: Vec<&str> = line.split("\t").collect();

        let gene: &str = content[0];
        let isoform: &str = content[1];
        isoforms.insert(isoform.to_string(), gene.to_string());
    }

    return Ok(isoforms);
} 



/// Get the coordinates of the first codon (start/stop).
/// If not in frame, return an empty codon.
fn find_first_codon(record: &BedRecord) -> Codon {
    let mut codon = Codon::new();
    let mut exon = 0;

    for k in 0..record.get_exon_frames().len() {
        if record.get_exon_frames()[k] >= 0 {
            exon = k;
            break;
        } else {
            return codon;
        }
    }

    let cds_exon_start = max(record.exon_start()[exon], record.cds_start());
    let cds_exon_end = min(record.exon_end()[exon], record.cds_end());

    let frame = if record.strand() == "+" {
        record.get_exon_frames()[exon]
    } else {
        (record.get_exon_frames()[exon] + (cds_exon_end - cds_exon_start)) % 3
    };

    if frame != 0 {
        return codon;
    }

    codon.start = record.cds_start();
    codon.end = record.cds_start() + (record.cds_end() - record.cds_start()).min(3);
    codon.index = exon as i32;

    if codon.end - codon.start < 3 {
        exon = exon + 1;
        if exon == record.exon_count() as usize {
            return codon;
        };
        let need = 3 - (codon.end - codon.start);
        if (record.cds_end() - record.cds_start()) < need {
            return codon;
        }
        codon.start2 = record.cds_start();
        codon.end2 = record.cds_start() + need;
    }
    codon
}



/// Get the coordinates of the last codon (start/stop).
/// If not in frame, return an empty codon.
fn find_last_codon(record: &BedRecord) -> Codon {
    let mut codon = Codon::new();
    let mut exon = 0;

    for k in (0..record.get_exon_frames().len()).step_by(1).rev() {
        if record.get_exon_frames()[k] >= 0 {
            exon = k;
            break;
        } else {
            return codon;
        }
    }

    let cds_exon_start = max(record.exon_start()[exon], record.cds_start());
    let cds_exon_end = min(record.exon_end()[exon], record.cds_end());

    let frame = if record.strand() == "-" {
        record.get_exon_frames()[exon]
    } else {
        (record.get_exon_frames()[exon] + (cds_exon_end - cds_exon_start)) % 3
    };

    if frame != 0 {
        return codon;
    }

    codon.start = max(record.cds_start(), record.cds_end() - 3);
    codon.end = record.cds_end();
    codon.index = exon as i32;

    if codon.end - codon.start < 3 {
        exon = exon + 1;
        if exon == record.exon_count() as usize {
            return codon;
        };
    
        let need = 3 - (codon.end - codon.start);
        if (record.cds_end() - record.cds_start()) < need {
            return codon;
        }
        codon.start2 = record.cds_start();
        codon.end2 = record.cds_start() + need;
    }
    codon
}



/// Check if all the bases of a codon are defined.
fn codon_complete(codon: &Codon) -> bool {
    ((codon.end - codon.start) + (codon.end2 - codon.start2)) == 3
}



/// Check if a given coordinate is within exon boundaries.
fn in_exon(record: &BedRecord, pos:i32, exon: usize) -> bool {
    (record.exon_start()[exon] <= pos) && (pos <= record.exon_end()[exon])
}



/// Move a position in an exon by a given distance, which is positive
/// to move forward and negative to move backwards. Introns are not
/// considered. If can't move distance and stay within exon boundaries,
/// panic.
fn move_pos(record: &BedRecord, pos: i32, dist: i32) -> i32 {
    let mut pos = pos;
    assert!(record.tx_start() <= pos && pos <= record.tx_end());
    let mut exon: Option<i16> = None;
    for i in 0..record.exon_count() {
        if in_exon(record, pos, i as usize) {
            exon = Some(i);
            break;
        }
    } 

    if exon.is_none() {
        panic!("Position {} not in exons", pos);
    }

    let mut steps = dist.abs();
    let direction = if dist >= 0 { 1 } else { -1 };

    while (0..record.exon_count()).contains(&(exon.unwrap())) && steps > 0 {
        if in_exon(record, pos + direction, exon.unwrap() as usize) {
            pos += direction;
            steps -= 1;
        } else if direction >= 0 {
            exon = Some(exon.unwrap() + 1);
            if let Some(ex) = exon {
                if (ex as usize) < record.exon_count() as usize {
                    pos = record.exon_start()[ex as usize];
                }
            }
        } else {
            exon = Some(exon.unwrap() - 1);
            if let Some(ex) = exon {
                if ex >= 0 {
                    pos = record.exon_end()[ex as usize] - 1;
                    steps -= 1;
                }
            }
        }
    }

    if steps > 0 {
        panic!("can't move {} by {}", pos, dist);
    }

    pos
} 



/// Build a "gene" feature line for a given group of transcripts.
/// Each line is unique for a given group.
fn build_gene_line(gene_name: &str, record: &BedRecord, file: &mut File) {
    assert!(gene_name.len() > 0);
    let gene_line = format!("{}\t{}\tgene\t{}\t{}\t.\t{}\t.\tID={};gene_id={}\n",
        record.chrom(),
        SOURCE,
        record.tx_start() + 1,
        record.tx_end(),
        record.strand(),
        gene_name,
        gene_name
    );
    file.write_all(gene_line.as_bytes()).unwrap();
}



/// Build a GTF line for a given feature (transcript, exon, CDS, five_prime_utr, three_prime_utr).
fn build_gtf_line(record: &BedRecord, 
    gene_name: &str, 
    feat_type: &str, 
    exon_start: i32, 
    exon_end: i32, 
    frame: i32, 
    exon: i16, 
    file: &mut File) {
    
    assert!(record.tx_start() < record.tx_end());

    let phase = match frame {
        -1 | -2 => ".",
        0 => "0",
        1 => "2",
        _ => "1",
    };
    

    let mut gtf_line = format!("{}\t{}\t{}\t{}\t{}\t.\t{}\t{}\t",
        record.chrom(),
        SOURCE,
        feat_type,
        exon_start + 1,
        exon_end,
        record.strand(),
        phase,
    );

    if feat_type == "transcript" {
        gtf_line += &format!("ID={};Parent={};gene_id={};transcript_id={}\n",
            record.name(), 
            gene_name, 
            gene_name, 
            record.name());
    } else {
        let prefix = match feat_type {
            "exon" => "exon",
            "CDS" => "CDS",
            "five_prime_utr" => "UTR5",
            "three_prime_utr" => "UTR3",
            "start_codon" => "start_codon",
            "stop_codon" => "stop_codon",
            _ => panic!("Unknown feature type {}", feat_type)
        };

        // Excludes UTRs
        if exon >= 0 {
            match record.strand() {
                "-" => {
                    gtf_line += &format!("ID={}:{}.{};Parent={};gene_id={};transcript_id={},exon_number={}\n", 
                    prefix,
                    record.name(), 
                    record.exon_count() - exon, 
                    record.name(), 
                    gene_name, 
                    record.name(),
                    record.exon_count() - exon);
                },
                "+" => {
                    gtf_line += &format!("ID={}:{}.{};Parent={};gene_id={};transcript_id={},exon_number={}\n", 
                    prefix,
                    record.name(),
                    exon + 1, 
                    record.name(), 
                    gene_name, 
                    record.name(),
                    exon + 1);
                },
                _ => panic!("Invalid strand {}", record.strand())
            }
        } else {
            gtf_line += &format!("ID={}:{};Parent={};gene_id={};transcript_id={}\n", 
            prefix,
            record.name(), 
            record.name(), 
            gene_name, 
            record.name());
        }
    }
    let _ = file.write_all(gtf_line.as_bytes());
}



/// Write the features of a given exon, including UTRs and CDS.
fn write_features(i: usize, 
    record: &BedRecord, 
    gene_name: &str, 
    first_utr_end: i32, 
    cds_start: i32, 
    cds_end: i32, 
    last_utr_start: i32, 
    frame: i32, 
    file: &mut File) {

    let exon_start = record.exon_start()[i];
    let exon_end = record.exon_end()[i];

    if exon_start < first_utr_end {
        let end = min(exon_end, first_utr_end);
        let utr_type = if record.strand() == "+" { "five_prime_utr" } else { "three_prime_utr" };
        build_gtf_line(record, gene_name, utr_type, exon_start, end, frame, -1, file);
    }

    if record.cds_start() < exon_end && exon_start < record.cds_end() {
        let start = max(exon_start, cds_start);
        let end = min(exon_end, cds_end);
        build_gtf_line(record, gene_name, "CDS", start, end, frame, i as i16, file);
    }

    if exon_end > last_utr_start {
        let start = max(exon_start, last_utr_start);
        let utr_type = if record.strand() == "+" { "three_prime_utr" } else { "five_prime_utr" };
        build_gtf_line(record, gene_name, utr_type, start, exon_end, frame, -1, file);
    }
}



/// Write the codon features (start/stop) for a given exon.
fn write_codon(record: &BedRecord, 
    gene_name: &str, 
    gene_type: &str, 
    codon: Codon, 
    file: &mut File) {

    build_gtf_line(record, 
        gene_name,
        gene_type,
        codon.start, 
        codon.end,
        0,
        codon.index as i16,
        file);

    if codon.start2 < codon.end2 {
        build_gtf_line(record, 
            gene_name, 
            gene_type, 
            codon.start, 
            codon.end, 
            codon.start2, 
            (codon.end - codon.start) as i16, 
            file);
    }
}



/// Convert a BED record to a GTF record.
fn to_gtf(record: &BedRecord, isoforms: &HashMap<String, String>, file: &mut File, gene_line: bool) {

    let gene_name = isoforms.get(record.name()).unwrap();
    let first_codon = find_first_codon(record);
    let last_codon = find_last_codon(record);

    let first_utr_end = record.cds_start();
    let last_utr_start = record.cds_end();

    let cds_end: i32 = if record.strand() == "+" && codon_complete(&last_codon) {
        move_pos(record, last_codon.end, -3)
    } else {record.cds_end()};

    let cds_start = if record.strand() == "-" && codon_complete(&first_codon) {
        move_pos(record, first_codon.start, 3)
    } else {record.cds_start()};

    if gene_line {build_gene_line(gene_name, record, file)};

    let _ = build_gtf_line(record, gene_name, "transcript", record.tx_start(), record.tx_end(), -1, -1, file);

    for i in 0..record.exon_count() as usize {
        build_gtf_line(record, gene_name, "exon", record.exon_start()[i], record.exon_end()[i], -1, i as i16, file);
        if cds_start < cds_end {
            write_features(i, record, gene_name, first_utr_end, cds_start, cds_end, last_utr_start, record.get_exon_frames()[i], file);
        }
    }

    match record.strand() {
        "+" => {
            if codon_complete(&first_codon) {
                write_codon(record, gene_name, "start_codon", first_codon, file);
            }
            if codon_complete(&last_codon) {
                write_codon(record, gene_name, "stop_codon", last_codon, file);
            }
        },
        "-" => {
            if codon_complete(&last_codon) {
                write_codon(record, gene_name, "start_codon", last_codon, file);
            }
            if codon_complete(&first_codon) {
                write_codon(record, gene_name, "stop_codon", first_codon, file);
            }
        },
        _ => panic!("Invalid strand {}", record.strand())
    }
}




fn bedsort(bed: &String) -> Result<Vec<(String, i32, String)>, ParseError> {

    let bedfile = File::open(PathBuf::from(bed)).unwrap();
    let reader = BufReader::new(bedfile);
    let mut tmp: Vec<(String, i32, String)> = Vec::new();

    for line in reader.lines() {
        let record = BedRecord::layer(&line?)?;
        tmp.push(record);
    }

    tmp.sort_by(|a, b| {
        let cmp_chr = compare(&a.0, &b.0);
        if cmp_chr == std::cmp::Ordering::Equal {
            a.1.cmp(&b.1)
        } else {
            cmp_chr
        }
    });

    Ok(tmp)
}



/// Convert a BED file to a GFF file.
/// ```
/// use bed2gff::bed2gff;
/// bed2gff("input.bed", "isoforms.txt", "output.gtf");
/// ```
pub fn bed2gff(input: &String, isoforms: &String, output: &String) -> Result<(), Box<dyn Error>> {

    msg();
    simple_logger::init_with_level(Level::Info)?;

    let start = Instant::now();
    
    let bed = bedsort(input).unwrap();
    let isoforms = get_isoforms(isoforms.into()).unwrap();
    let mut output = File::create(PathBuf::from(output)).unwrap();
    let mut seen_genes: HashSet<String> = HashSet::new();

    let _ = comments(&mut output);

    for line in bed {
        let record = BedRecord::new(&line.2);

        if let Ok(record) = record {
            let key = match isoforms.get(record.name()) {
                Some(gene) => Ok(gene),
                None => {
                    log::error!("Isoform {} not found in isoforms file.", &record.name().bright_red().bold());
                    Err("Isoform not found in isoforms file")
                }
            };
            
            if key.is_err() {
                println!("{} {}", 
                "Fail:".bright_red().bold(),
                "BED file could not be converted. Please check your isoforms file.");
                std::process::exit(1);
            }

            if !seen_genes.contains(key?) {
                seen_genes.insert(key?.to_string());
                let _ = to_gtf(&record, &isoforms, &mut output, true);
            } else {
                let _ = to_gtf(&record, &isoforms, &mut output, false);
            };
        } else {
            log::error!("Failed to parse a BedRecord.");
        };
    }

    let peak_mem = PEAK_ALLOC.peak_usage_as_mb();

    log::info!("Memory usage: {} MB", peak_mem);
    log::info!("Elapsed: {:.4?} secs", start.elapsed().as_secs_f32());

    Ok(())
}



fn msg() {
    println!("{}\n{}",
        "\n##### BED2GFF #####".bright_blue().bold(),
        indoc!("A Rust BED-to-GFF translator.
        Repository: https://github.com/alejandrogzi/bed2gff
        Feel free to contact the developer if any issue/suggest/bug is found.
        "));
}


fn get_date() -> String {
    let now = chrono::Utc::now();
    let year = now.year();
    let month = now.month();
    let day = now.day();

    format!("#{}-{}-{}", year, month, day)
}


fn comments(file: &mut File) {
    let _ = file.write_all(format!("{}\n", GFF3).as_bytes());
    let _ = file.write_all(format!("#provider: {}\n", PROVIDER).as_bytes());
    let _ = file.write_all(format!("#version: {}\n", VERSION).as_bytes());
    let _ = file.write_all(format!("#contact: {}\n", REPOSITORY).as_bytes());
    let _ = file.write_all(format!("#date: {}\n", get_date()).as_bytes());
}