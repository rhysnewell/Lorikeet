use bird_tool_utils::external_command_checker::*;
use std;
use std::io::Read;

pub fn check_for_bwa() {
    check_for_external_command_presence("BWA", "which bwa").expect("Failed to find installed BWA");
}

pub fn check_for_bwa_mem2() {
    check_for_external_command_presence("BWA", "which bwa-mem2")
        .expect("Failed to find installed BWA");
    default_version_check("bwa-mem2", "2.0", false, Some("bwa-mem2 version"))
        .expect("Failed to find sufficient version of bwa-mem2");
}

pub fn check_for_samtools() {
    check_for_external_command_presence("samtools", "which samtools")
        .expect("Failed to find installed samtools");
    default_version_check("samtools", "1.9", false, None)
        .expect("Failed to find sufficient version of samtools");
}

pub fn check_for_bcftools() {
    check_for_external_command_presence("bcftools", "which bcftools")
        .expect("Failed to find installed bcftools");
}

pub fn check_for_prodigal() {
    check_for_external_command_presence("prodigal", "which prodigal")
        .expect("Failed to find installed prodigal");
}

pub fn check_for_svim() {
    check_for_external_command_presence("svim", "which svim")
        .expect("Failed to find installed svim");
}

pub fn check_for_svim_asm() {
    check_for_external_command_presence("svim-asm", "which svim-asm")
        .expect("Failed to find installed svim-asm");
}

pub fn check_for_minimap2() {
    check_for_external_command_presence("minimap2", "which minimap2")
        .expect("Failed to find installed minimap2");
    default_version_check("minimap2", "2.24-r1122", false, None)
        .expect("Failed to find sufficient version of minimap2");
}

pub fn check_for_ngmlr() {
    check_for_external_command_presence("ngmlr", "which ngmlr")
        .expect("Failed to find ngmlr installed");
}

pub fn check_for_pggb() {
    check_for_external_command_presence("pggb", "which pggb")
        .expect("Failed to find pggb installed");
    info!("Valid pggb installation found")
}

pub fn check_for_fastani() {
    check_for_external_command_presence("fastaANI", "which fastANI")
        .expect("Failed to find installed fastANI");
    default_version_check("fastANI", "1.31", false, None)
        .expect("Failed to find sufficient version of fastANI");
}

pub fn check_for_dashing() {
    check_for_external_command_presence("dashing", "which dashing")
        .expect("Failed to find installed dashing. You may wish to use the finch precluster method if you are having problems with dashing.");
    default_version_check("dashing", "0.4.0", true, None)
        .expect("Failed to find sufficient version of dashing. You may wish to use the finch precluster method if you are having problems with dashing.");
}
