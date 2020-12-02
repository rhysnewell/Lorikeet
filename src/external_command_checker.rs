use std;
use std::io::Read;

pub fn check_for_bwa() {
    self::check_for_external_command_presence("BWA", "which bwa");
}

pub fn check_for_samtools() {
    self::check_for_external_command_presence("samtools", "which samtools");
}

pub fn check_for_samclip() {
    self::check_for_external_command_presence("samclip", "which samclip");
}

pub fn check_for_bcftools() {
    self::check_for_external_command_presence("bcftools", "which bcftools");
}

pub fn check_for_vt() {
    self::check_for_external_command_presence("vt", "which vt");
}

pub fn check_for_prokka() {
    self::check_for_external_command_presence("prokka", "which prokka");
}

pub fn check_for_prodigal() {
    self::check_for_external_command_presence("prodigal", "which prodigal");
}

pub fn check_for_pilon() {
    self::check_for_external_command_presence("pilon", "which pilon");
}

pub fn check_for_snippy() {
    self::check_for_external_command_presence("snippy", "which snippy");
}

pub fn check_for_svim() {
    self::check_for_external_command_presence("svim", "which svim");
}

pub fn check_for_svim_asm() {
    self::check_for_external_command_presence("svim-asm", "which svim-asm");
}

pub fn check_for_sniffles() {
    self::check_for_external_command_presence("sniffles", "which sniffles");
}

pub fn check_for_gatk() {
    self::check_for_external_command_presence("gatk", "which gatk");
}

pub fn check_for_freebayes() {
    self::check_for_external_command_presence("freebayes", "which freebayes");
}

pub fn check_for_freebayes_parallel() {
    self::check_for_external_command_presence("freebayes-parallel", "which freebayes-parallel");
}

pub fn check_for_fasta_generate_regions() {
    self::check_for_external_command_presence(
        "fasta_generate_regions.py",
        "which fasta_generate_regions.py",
    )
}

pub fn check_for_minimap2() {
    self::check_for_external_command_presence("minimap2", "which minimap2");
}

pub fn check_for_ngmlr() {
    self::check_for_external_command_presence("ngmlr", "which ngmlr");
}

fn check_for_external_command_presence(executable_name: &str, testing_cmd: &str) {
    debug!("Checking for {} ..", executable_name);
    let mut cmd = std::process::Command::new("bash");
    cmd.arg("-c")
        .arg(testing_cmd)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    let mut process = cmd.spawn().expect("Unable to execute bash");
    let es = process.wait().expect(&format!(
        "Failed to glean exitstatus while checking for presence of {}",
        executable_name
    ));
    if !es.success() {
        error!(
            "Could not find an available {} executable.",
            executable_name
        );
        let mut err = String::new();
        process
            .stderr
            .expect("Failed to grab stderr from failed executable finding process")
            .read_to_string(&mut err)
            .expect("Failed to read stderr into string");
        error!("The STDERR was: {:?}", err);
        panic!(
            "Cannot continue without {}. Testing for presence with `{}` failed",
            executable_name, testing_cmd
        );
    }
}
