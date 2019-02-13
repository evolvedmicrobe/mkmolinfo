extern crate clap;
use clap::{Arg, App};


use std::path::Path;
use std::collections::{HashMap};

extern crate ndarray;
//use ndarray::Array;
//use ndarray::Zip;
//use ndarray::Array1;

use rust_htslib::bam;
use rust_htslib::bam::Read;
//use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::{Aux,Record};

//extern crate failure;
//#[macro_use] extern crate failure_derive;

extern crate h5;
use h5::types::FixedAscii;
//use h5::File;

//use hdf5_rs::prelude::*;

#[derive(Debug, Default)]
struct barcode_counts {
    //barcode: u64,
    //umi: u64,
    // feature_mapped: String,
    total_reads: u32, // these are total confidentally mapped
    barcode_corrected_reads: u32,
    conf_mapped: u32,
    nonconf_mapped_reads: u32,
    umi_corrected_reads: u32,
    unmapped_reads: u32
}

fn barcode_str_to_u64(arr : &[u8]) -> u64 {
    let mut result : u64 = 0;
    for nuc in arr.iter() {
        result = result << 2;
        let coding = match nuc {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                println!("{}", String::from_utf8_lossy(arr).to_string());
                panic!("Barcode contained non-ACGT character");
            }
        };
        result |= coding;
    }
    return result
}

fn combine_bc_and_umi(bc: u64, umi:u64) -> u64 {
    let high_end = bc << 25;
    let combined = high_end | umi;
    return combined;
}




fn parse_bam_for_metrics(bamname : &Path, bc_recs: &mut HashMap<u64, barcode_counts>, verbose : bool) -> () {
    let mut in_bam = bam::Reader::from_path(bamname).expect("Failed to open BAM");
    let mut cnt: i32 = 1;
    let mut record: Record = Record::new();


    while let Ok(_) = in_bam.read(&mut record) {
        cnt += 1;
        if (cnt % 10000000 == 0) & verbose {
            println!("Parsed {} BAM records", cnt);
        }
       
        if let Some(Aux::String(barcode)) = record.aux(b"CB") {
            if barcode[barcode.len() - 1] != b'1' {
                println!("Error: Barcode is from multiple GEM groups, only one is supported. Please email the author to enable this. BC: {}", String::from_utf8_lossy(barcode).to_string());
                std::process::exit(1);
            }
            if let Some(Aux::String(umi)) = record.aux(b"UB") {
                let trimmed_bc = &barcode[..barcode.len() - 2];
                let bc_i = barcode_str_to_u64(&trimmed_bc); // cut-off "-1"
                let umi_i = barcode_str_to_u64(&umi);
                let combined = combine_bc_and_umi(bc_i, umi_i);
                if bc_recs.contains_key(&combined) {

                    // if this wasn't in our list, forget about it
                    // if let Some(ref cnts) = bc_recs.get_mut(&bc_i) {
                   let cnts = bc_recs.entry(combined).or_default();

/*                    // See if feature matches
                    if let Some(Aux::String(fx)) = record.aux(b"fx") {
                        let fxs = String::from_utf8_lossy(fx).to_string();
                        //println!("{} - {}", fxs, cnts.feature_mapped);
                        if fxs == cnts.feature_mapped {
                            cnts.total_reads += 1;
                        }
                    }*/

                    cnts.total_reads += 1;
                    //if !record.is_duplicate() && !record.is_secondary() {
                    // Check for barcode correction
                    if let Some(Aux::String(raw_barcode)) = record.aux(b"CR") {
                        let bc_matches = trimmed_bc.iter().zip(raw_barcode.iter()).all(|(x, y)| x == y);
                        if !bc_matches {
                            cnts.barcode_corrected_reads += 1;
                        }
                    }

                    // Check for UMI correction
                    if let Some(Aux::String(raw_umi)) = record.aux(b"UR") {
                        let r_umi_i = barcode_str_to_u64(raw_umi);
                        if r_umi_i != umi_i {
                            cnts.umi_corrected_reads += 1;
                        }
                    }

                    // conf mapped, unmapped and nonconf mapped
                    if record.tid() == -1 {
                        cnts.unmapped_reads += 1;
                    } else if record.mapq() != 255 {
                        cnts.nonconf_mapped_reads += 1;
                    } else if let Some(_) = record.aux(b"TX") {
                        if record.mapq() == 255 {
                            cnts.conf_mapped += 1;
                        }
                    }
                }
            }
        }
    }
}

fn write_array(arr :&Vec<u32>, name : &str, h5f : &h5::Group) {
    let ds1 = h5f.new_dataset::<u32>()
    .shuffle(true).gzip(1).chunk_infer()
    .create(name, arr.len()).unwrap();
    //let arr = Array1::from(arr);
    ds1.write(arr.as_slice()).unwrap();
}

fn write_array64(arr :&Vec<u64>, name : &str, h5f : &h5::Group) {
    let ds1 = h5f.new_dataset::<u64>()
        .shuffle(true).gzip(1).chunk_infer()
        .create(name, arr.len()).unwrap();
    //let arr = Array1::from(arr);
    ds1.write(arr.as_slice()).unwrap();
}


fn main() {
    let matches = App::new("mkmolinfo")
        .version("1.0")
        .author("Nigel Delaney <nigel.delaney@10xgenomics.com>")
        .about("Takes a path to a Cell Ranger Output directory and creates an older (Cell Ranger < 3.0) style molecule_info.h5 file.")
        .arg(Arg::with_name("DIR")
            .value_name("OUT_DIRECTORY")
            .help("Specify the output directory produced by Cell Ranger")
            .required(true)
            .index(1))
        .arg(Arg::with_name("OUTPUT")
            .short("o")
            .long("output")
            .value_name("FILE")
            .help("Specify an output file name (default molecule_info_new.h5)")
            .required(false))
        .arg(Arg::with_name("v")
            .short("v")
            .multiple(true)
            .help("Print progress as program parses through BAM file."))
        .get_matches();

    let outfile = matches.value_of("OUTPUT").unwrap_or("molecule_info_new.h5");
    let out_dir = matches.value_of("DIR").unwrap();
    let top = Path::new(out_dir).join("outs");
    let vcount = matches.occurrences_of("v");
    let verbose = vcount > 0;
    let bam_counts_out = vcount > 1;
    let bam_path = top.join("possorted_genome_bam.bam");
    let h5_path = top.join("molecule_info.h5");

    std::fs::copy(&h5_path, &outfile).expect("Could not copy original h5 file");
    let h5f = h5::File::open(outfile, "r+").expect("Could not open hdf5 output file (copy of original).");
    //let a = f.group("/").unwrap();



    // Read through the hdf5 and define each feature a read maps to.

    let bcs_i =  {
        let bcs = h5f.dataset("/barcode_idx").expect("h5 file did not contain barcode_idx");
        bcs.read_1d::<u64>().expect("Could not read in barcodes.")
    };
    let umi_i =  {
        let umis = h5f.dataset("/umi").expect("H5 file did not contain umi");
        umis.read_1d::<u64>().expect("Could not read in UMIs.")
    };

    let n_barcodes : usize = bcs_i.len();

    // load barcode strings
    let barcode_str =  {
        const FBC : usize = 18;
        let bc_strings = h5f.dataset("/barcodes").expect("H5 file did not contain barcodes");
        bc_strings.read_1d::<FixedAscii<[u8;FBC]>>().expect("Could not read in barcode strings")
    };

/*    // Load feature strings
    let feature_idx = {
        let feature_i = h5f.dataset("/feature_idx").expect("Could not find feature_idx in H5 file");
        feature_i.read_1d::<u32>().expect("Could not read in feature idx")
    };

    let feature_str =  {
        let feature_h5 = h5f.dataset("/features/id").expect("Could not find feature IDs in H5 file.");
        feature_h5.read_1d::<FixedAscii<[u8;25]>>().expect("Could not read in feature strings")
    };*/

    let mut bc_translate = {
        let mut maker = vec![0u64; barcode_str.len()];
        for i in 0..barcode_str.len() {
            let bcbytes = barcode_str[i].as_bytes();
            if bcbytes.len() > 16 {
                println!("Error: Barcode is from multiple GEM groups, only one is supported:.  Please email the author to enable this. BC: {}", String::from_utf8_lossy(bcbytes).to_string());
                std::process::exit(1);
            }
            let bc_t = barcode_str_to_u64(bcbytes);
            maker[i] = bc_t;
        }
        maker
    };

    // Initialize map with barcodes of interest
    let mut metrics: HashMap<u64, barcode_counts> = HashMap::with_capacity(n_barcodes);

    for i in 0..bcs_i.len() {
        let bci = bcs_i[i] as usize;
        let bct = bc_translate[bci];
        let key = combine_bc_and_umi(bct as u64, umi_i[i] as u64);
        //let feature_index = feature_idx[i] as usize;
        let cnt = barcode_counts {
            //feature_mapped : String::from_utf8_lossy(feature_str[feature_index].as_bytes()).to_string(),
            total_reads: 0,
            barcode_corrected_reads: 0,
            conf_mapped: 0,
            nonconf_mapped_reads: 0,
            umi_corrected_reads: 0,
            unmapped_reads: 0,
        };
        metrics.insert(key, cnt);
    }


    parse_bam_for_metrics(bam_path.as_path(), &mut metrics, verbose);
    

    if (verbose) {
        println!("BAM parsing produced metrics for {} barcodes/UMIs.", metrics.len());
    }
    // Create vectors to output
    let mut non_conf_mapped = vec![0u32; n_barcodes];
    let mut conf_mapped = vec![0u32; n_barcodes];
    let mut bc_corrected = vec![0u32; n_barcodes];
    let mut umi_corrected = vec![0u32; n_barcodes];
    let mut unmapped = vec![0u32; n_barcodes];
    let mut reads= vec![0u32; n_barcodes];
    let mut barcodes = vec![0u64; n_barcodes];

    let mut found = 0;
    let mut notfound = 0;
    //for (i, (bc, umi)) in bcs_i.zip(umi_i).enumerate() {
    for i in 0..bcs_i.len() {
        let barcode_index = bcs_i[i] as usize;
        let bc = bc_translate[barcode_index];
        let umi = umi_i[i];
        let hash_key = combine_bc_and_umi(bc, umi);
        if !metrics.contains_key(&hash_key) {
            println!("Logic error! Previously entered BC no longer in hashmap.");
            std::process::exit(1);
        } else {
            let cnts = metrics.get(&hash_key).unwrap();
            if cnts.total_reads > 0 {
                found +=1;
            } else {
                notfound +=1;
            }
            non_conf_mapped[i] = cnts.nonconf_mapped_reads;
            conf_mapped[i] = cnts.conf_mapped;
            bc_corrected[i] = cnts.barcode_corrected_reads;
            umi_corrected[i] = cnts.umi_corrected_reads;
            unmapped[i] = cnts.unmapped_reads;
            reads[i] = cnts.total_reads;
            barcodes[i] = bc;
        }
    }
    if verbose {
        println!("Barcodes in original molecule_info.h5 found in BAM = {}, not found = {}", found, notfound);
    }
    let a = h5f.group("/").unwrap();
    write_array(&non_conf_mapped, "nonconf_mapped_reads", &a);
    write_array(&conf_mapped, "conf_mapped", &a); // This name is different
    write_array(&umi_corrected, "umi_corrected_reads", &a);
    write_array(&bc_corrected, "barcode_corrected_reads", &a);
    write_array(&unmapped, "unmapped_reads", &a);
    write_array(&reads, "reads", &a);
    write_array64(&barcodes, "barcode", &a);
}

#[cfg(test)]
mod tests {
    #[test]
    fn bc_to_umi() {
        let arr = "ACGTAAACTTTTCGAT".as_bytes();
        let result = super::barcode_str_to_u64(&arr);
        assert_eq!(453115747, result);
    }
}

/*  let in_bam_res = bam::Reader::from_path(&path);
let mut bam = match in_bam_res {
    Ok(f) => f,
    Err(error) => {
        match error {
            bam::ReaderPathError::InvalidPath => {
                println!("File {} was not found.", bam_name);
                std::process::exit(1);
            },
            bam::ReaderPathError::BGZFError {..}=> {
                println!("There was a problem reading the compressed file: {}", bam_name);
                std::process::exit(1);
            }
        };
    }
};
*/

