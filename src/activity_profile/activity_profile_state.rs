use num::traits::Float;
use utils::simple_interval::{Locatable, SimpleInterval};

/**
 * Captures the probability that a specific locus in the genome represents an "active" site containing
 * real variation.
 */
#[derive(Clone, Debug, PartialOrd, PartialEq)]
pub enum Type<T: Float + Copy> {
    None,
    HighQualitySoftClips(T),
}

impl<T: Float + Copy> Type<T> {
    pub fn new(high_quality_soft_clips_mean: T, threshold: T) -> Type<T> {
        if high_quality_soft_clips_mean >= threshold {
            Type::HighQualitySoftClips(high_quality_soft_clips_mean)
        } else {
            Type::None
        }
    }
}

#[derive(Clone, Debug, PartialOrd, PartialEq)]
pub struct ActivityProfileState {
    loc: SimpleInterval,
    pub active_prob: f64,
    result_state: Type<f64>,
}

impl ActivityProfileState {
    // When range-checking probabilities, we allow this much tolerance.
    pub const PROBABILITY_TOLERANCE: f64 = 0.01;

    /**
     * Create a new ActivityProfileState at loc with probability of being active of activeProb that maintains some
     * information about the result state and value
     *
     * The only state value in use is HighQualitySoftClips, and here the value is interpreted as the number
     * of bp affected by the soft clips.
     *
     * @param loc the position of the result profile (for debugging purposes)
     * @param activeProb the probability of being active (between 0 and 1)
     */
    pub fn new(
        loc: SimpleInterval,
        active_prob: f64,
        result_state: Type<f64>,
    ) -> ActivityProfileState {
        if loc.size() != 1 {
            panic!(
                "Location for an ActivityProfileState must have to size 1 bp but saw {:?}",
                loc
            )
        };

        ActivityProfileState {
            loc,
            active_prob,
            result_state,
        }
    }

    pub fn is_active_prob(&self) -> f64 {
        self.active_prob
    }

    /**
     * Set the probability that this site is active.
     *
     * Probabilities should be between 0.0 and 1.0, however this is not currently enforced
     * because the {@link BandPassActivityProfile} can sometimes generate probabilities that
     * slightly exceed 1.0 when moving probability mass around. We intend to fix this by
     * capping at 1.0, but first we must evaluate the effects of capping on the HaplotypeCaller.
     *
     * @param activeProb probability (should be between 0.0 and 1.0) that the site is active
     */
    pub fn set_is_active_prob(&mut self, active_prob: f64) {
        self.active_prob = active_prob
    }

    pub fn get_result_state(&self) -> &Type<f64> {
        &self.result_state
    }

    pub fn get_result_value(&self) -> f64 {
        match self.result_state {
            Type::None => 0.0,
            Type::HighQualitySoftClips(count) => count,
        }
    }

    /**
     * The offset of state w.r.t. our current region's start location
     * @param regionStartLoc the start of the region, as a Locatable
     * @return the position of this profile relative to the start of this region
     */
    pub fn get_offset(&self, region_start_loc: &SimpleInterval) -> i64 {
        (self.loc.get_start() as i64) - (region_start_loc.get_start() as i64)
    }

    pub fn get_loc(&self) -> &SimpleInterval {
        &self.loc
    }
}
