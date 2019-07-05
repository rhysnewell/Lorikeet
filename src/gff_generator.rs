//use std;
//use std::io;
//use std::io::Read;
//use nix::unistd;
//use nix::sys::stat;
//use tempdir::TempDir;
//use tempfile;
//use bio::io::gff;
//
//
//pub trait NamedGffReader {
//    // Name of the stoit
//    fn name(&self) -> &str;
//
//    // Read a record into record parameter
//    fn read(&mut self, record: &mut gff::record::Record) -> ();
//
//    // Return the bam header of the final BAM file
//    fn header(&self) -> &bam::HeaderView;
//
//    fn finish(self);
//
//    // Number of reads that were detected
//    fn num_detected_primary_alignments(&self) -> u64;
//}
//
//pub trait NamedBamReaderGenerator<T> {
//    // For readers that map, start the process of mapping
//    fn start(self) -> T;
//}
//
//pub struct BamFileNamedReader {
//    stoit_name: String,
//    bam_reader: bam::Reader,
//    num_detected_primary_alignments: u64,
//}
//
//impl NamedBamReader for BamFileNamedReader {
//    fn name(&self) -> &str {
//        &(self.stoit_name)
//    }
//    fn read(&mut self, record: &mut bam::record::Record) -> () {
//        let res = self.bam_reader.read(record);
//        if res.is_ok() && !record.is_secondary() && !record.is_supplementary() {
//            self.num_detected_primary_alignments += 1;
//        }
//        return res;
//    }
//    fn header(&self) -> &bam::HeaderView {
//        self.bam_reader.header()
//    }
//    fn finish(self) {;}
//
//    fn num_detected_primary_alignments(&self) -> u64 {
//        return self.num_detected_primary_alignments
//    }
//}