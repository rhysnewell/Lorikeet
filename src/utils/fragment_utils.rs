use std::cmp::min;

use crate::reads::bird_tool_reads::BirdToolRead;
use crate::reads::cigar_utils::CigarUtils;
use crate::reads::read_utils::ReadUtils;
use crate::utils::quality_utils::QualityUtils;
use crate::utils::simple_interval::Locatable;

pub const DEFAULT_PCR_SNV_ERROR_RATE: f64 = 1e-4;
lazy_static! {
    static ref DEFAULT_PCR_SNV_ERROR_QUAL: u8 =
        QualityUtils::get_phred_score_from_obs_and_errors(DEFAULT_PCR_SNV_ERROR_RATE);
    static ref HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL: u8 = *DEFAULT_PCR_SNV_ERROR_QUAL / 2;
}

/**
 * Fix two overlapping reads from the same fragment by adjusting base qualities, if possible
 *
 *  Looks at the bases and alignment, and tries its best to create adjusted base qualities so that the observations
 * are not treated independently.  Sets the qualities of firstRead and secondRead to mimic a merged read or
 * nothing if the algorithm cannot create a meaningful one
 * @param pair two overlapping paired reads
 * @param setConflictingToZero if true, set base qualities to zero when mates have different base at overlapping position
 * @param halfOfPcrSnvQual half of phred-scaled quality of substitution errors from PCR. May not be negative.
 * @param halfOfPcrIndelQual half of phred-scaled quality of indel errors from PCR. May not be negative.
 */
pub fn adjust_quals_of_overlapping_paired_fragments(
    pair: (BirdToolRead, BirdToolRead),
    _set_conflicting_to_zero: bool,
    half_of_pcr_snv_qual: Option<u8>,
    _half_of_pcr_indel_qual: Option<u8>,
) -> (BirdToolRead, BirdToolRead) {
    let in_order = pair.0.get_soft_start().unwrap() < pair.1.get_soft_start().unwrap();
    let (mut first_read, mut second_read) = if in_order {
        (pair.0, pair.1)
    } else {
        (pair.1, pair.0)
    };

    // don't adjust fragments that do not overlap
    if first_read.get_end() < second_read.get_start()
        || first_read.get_contig() != second_read.get_contig()
    {
        return (first_read, second_read);
    }

    // the offset and cigar operator in the first read at the start of the left read
    let (offset, operator) = ReadUtils::get_read_index_for_reference_coordinate(
        first_read.get_soft_start().unwrap(),
        first_read.read.cigar(),
        second_read.get_start(),
    );
    // let operator = offset_and_operator.1;
    if offset.is_none()
        || match &operator {
            Some(op) => CigarUtils::is_clipping(op),
            None => false,
        }
    {
        return (first_read, second_read);
    }

    // Compute the final aligned base indexes for both since there might be right base softclips
    let first_read_end_base = ReadUtils::get_read_index_for_reference_coordinate(
        first_read.get_soft_start().unwrap(),
        first_read.read.cigar(),
        first_read.get_end(),
    )
    .0
    .unwrap();
    let second_read_end_base = ReadUtils::get_read_index_for_reference_coordinate(
        second_read.get_soft_start().unwrap(),
        second_read.read.cigar(),
        second_read.get_end(),
    )
    .0
    .unwrap();

    // TODO: we should be careful about the case where {@code operator} is a deletion; that is, when the second read start falls in a deletion of the first read
    // TODO: however, the issue is bigger than simply getting the start correctly, because the code below assumes that all bases of both reads are aligned in their overlap.
    // TODO: Any indel that occurs in one read and not the other will spoil things.  Really, the correct thing to do is a Smith-Waterman (or other) alignment of the reads
    // TODO: in their overlap and correct the double-counting for all aligned bases.
    // TODO: a cheaper solution would be to cap all quals in the overlap region by half of the PCR qual.
    let first_read_stop = offset.unwrap();
    let second_offset = ReadUtils::get_read_index_for_reference_coordinate(
        second_read.get_soft_start().unwrap(),
        second_read.read.cigar(),
        second_read.get_start(),
    )
    .0
    .unwrap(); //This operation handles softclipped bases in the qual/base array
    let num_overlapping_bases = min(
        first_read_end_base.saturating_sub(first_read_stop),
        second_read_end_base.saturating_sub(second_offset),
    ) + 1; // Add 1 here because if R1 ends on the same base that R2 starts then there is 1 base of overlap not 0

    let mut first_read_quals = first_read.read.qual().to_vec();
    let mut second_read_quals = second_read.read.qual().to_vec();

    let half_of_pcr_error_qual = match half_of_pcr_snv_qual {
        Some(val) => val,
        None => *HALF_OF_DEFAULT_PCR_SNV_ERROR_QUAL,
    };

    for i in 0..num_overlapping_bases {
        let first_read_index = first_read_stop + i;
        let second_read_index = second_offset + i;
        let first_read_base = first_read.bases[first_read_index];
        let second_read_base = second_read.bases[second_read_index];

        if first_read_base == second_read_base {
            first_read_quals[first_read_index] =
                min(first_read_quals[first_read_index], half_of_pcr_error_qual);
            second_read_quals[second_read_index] =
                min(second_read_quals[second_read_index], half_of_pcr_error_qual);
        } else {
            first_read_quals[first_read_index] = 0;
            second_read_quals[second_read_index] = 0;
        }
    }

    let first_read_cigar = first_read.read.cigar().take();
    let first_read_name = first_read.read.qname().to_vec();
    first_read.update(
        first_read_name.as_slice(),
        Some(&first_read_cigar),
        first_read.bases.clone(),
        first_read_quals.as_slice(),
    );

    let second_read_cigar = second_read.read.cigar().take();
    let second_read_name = second_read.read.qname().to_vec();
    second_read.update(
        second_read_name.as_slice(),
        Some(&second_read_cigar),
        second_read.bases.clone(),
        second_read_quals.as_slice(),
    );

    return (first_read, second_read);
    // match half_of_pcr_indel_qual {
    //     None => {
    //         // pass
    //     },
    //     Some(max_indel_qual) => {
    //
    //     }
    // }
}
