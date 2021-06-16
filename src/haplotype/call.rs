// use rust_htslib::bam;
// use rust_htslib::bam::record::{Record, Aux, Cigar};
// use bio::io::{fasta, bed};
//
// pub fn len(cig: &Cigar) -> u32 {
//     match *cig {
//         Cigar::Match(l) => l,
//         Cigar::Ins(l) => l,
//         Cigar::Del(l) => l,
//         Cigar::RefSkip(l) => l,
//         Cigar::SoftClip(l) => l,
//         Cigar::HardClip(l) => l,
//         Cigar::Pad(l) => l,
//         Cigar::Equal(l) => l,
//         Cigar::Diff(l) => l,
//         Cigar::Back(l) => l,
//     }
// }
//
// /// Find the position in the read of `ref_pos`.
// pub fn find_ref_pos_in_read(ref_pos: u32, rec: &Record) -> usize {
//
//     let mut read_pos = 0;
//     let mut tpl_pos = rec.pos() as u32;
//
//     if tpl_pos > ref_pos + 300 {
//         return 0;
//     }
//
//     // Avoid wrap-around at start of contig
//     //if tpl_pos < max(0, ref_pos as isize - 300) as u32 {
//     //    return rec.seq().len();
//     //}
//
//     for cig_elem in &rec.cigar() {
//
//         match *cig_elem {
//             Cigar::Match(l) | Cigar::Equal(l) | Cigar::Diff(l) => {
//                 read_pos += l;
//                 tpl_pos += l;
//             }
//
//             Cigar::Ins(l) => read_pos += l,
//
//             Cigar::Del(l) => tpl_pos += l,
//             Cigar::SoftClip(l) => read_pos += l,
//             _ => (),
//         }
//
//         if tpl_pos > ref_pos {
//             return max(0_i32, (read_pos as i32) - (tpl_pos - ref_pos) as i32) as usize;
//         }
//     }
//     rec.seq().len()
// }
//
//
// /// Return sequence of rec that aligns within the requested locus
// #[inline(never)]
// pub fn get_seq_bounded(rec: &Record, locus: &Locus) -> Vec<u8> {
//     let r_start = find_ref_pos_in_read(locus.start, rec);
//     let mut r_end = find_ref_pos_in_read(locus.end, rec);
//
//     // FIXME: this is a hack to prevent the interval from going negative
//     // there is an issue in find_ref_pos_in_read
//     r_end = max(r_end, r_start);
//
//     let rec_seq = rec.seq();
//     let mut seq = Vec::with_capacity(r_end - r_start);
//     for i in r_start..r_end {
//         seq.push(rec_seq[i]);
//     }
//
//     seq
// }
//
// #[inline(never)]
// pub fn get_qual_bounded(rec: &Record, locus: &Locus) -> Vec<u8> {
//     let r_start = find_ref_pos_in_read(locus.start, rec);
//     let r_end = find_ref_pos_in_read(locus.end, rec);
//
//     let rec_qual = rec.qual();
//     let mut qual = Vec::with_capacity(r_end - r_start);
//     for i in r_start..r_end {
//         qual.push(rec_qual[i]);
//     }
//
//     qual
// }
//
// pub fn print_seq(s: &Vec<u8>) -> String {
//     let mut bs = String::new();
//     for b in s.iter() {
//         bs.push(*b as char);
//     }
//     bs
// }
