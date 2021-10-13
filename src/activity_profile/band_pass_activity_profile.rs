use activity_profile::activity_profile::{ActivityProfile, Profile};
use activity_profile::activity_profile_state::{ActivityProfileState, Type};
use assembly::assembly_region::AssemblyRegion;
use rayon::prelude::*;
use utils::math_utils::MathUtils;
use utils::simple_interval::{Locatable, SimpleInterval};

/**
 * A band pass filtering version of the activity profile
 *
 * Applies a band pass filter with a Gaussian kernel to the input state probabilities to smooth
 * them out of an interval
 */
#[derive(Debug)]
pub struct BandPassActivityProfile {
    filter_size: usize,
    sigma: f64,
    pub gaussian_kernel: Vec<f64>,
    pub activity_profile: ActivityProfile,
}

impl BandPassActivityProfile {
    pub const MAX_FILTER_SIZE: usize = 50;
    pub const MIN_PROB_TO_KEEP_IN_FILTER: f64 = 1e-5;
    pub const DEFAULT_SIGMA: f64 = 17.0;

    /**
     * Create an activity profile that implements a band pass filter on the states
     *  @param maxProbPropagationDistance region probability propagation distance beyond it's maximum size
     * @param activeProbThreshold  threshold for the probability of a profile state being active
     * @param maxFilterSize the maximum size of the band pass filter we are allowed to create, regardless of sigma
     * @param sigma the variance of the Gaussian kernel for this band pass filter
     * @param adaptiveFilterSize if true, use the kernel itself to determine the best filter size
     */
    pub fn new(
        max_prob_propagation_distance: usize,
        active_prob_threshold: f64,
        max_filter_size: usize,
        sigma: f64,
        adaptive_filter_size: bool,
        ref_idx: usize,
        tid: usize,
        contig_len: usize,
    ) -> BandPassActivityProfile {
        let activity_profile = ActivityProfile::new(
            max_prob_propagation_distance,
            active_prob_threshold,
            ref_idx,
            tid,
            contig_len,
        );
        let mut full_kernel = Self::make_kernel(max_filter_size, sigma);

        let filter_size = if adaptive_filter_size {
            Self::determine_filter_size(&full_kernel, Self::MIN_PROB_TO_KEEP_IN_FILTER)
        } else {
            max_filter_size
        };

        BandPassActivityProfile {
            activity_profile,
            sigma,
            filter_size,
            gaussian_kernel: Self::make_kernel(filter_size, sigma),
        }
    }

    fn make_kernel(filter_size: usize, sigma: f64) -> Vec<f64> {
        let band_size = 2 * filter_size + 1;
        // let mut kernel = vec![0.0; band_size];
        let mut kernel = (0..band_size)
            .into_par_iter()
            .map(|iii| MathUtils::normal_distribution(filter_size as f64, sigma, iii as f64))
            .collect::<Vec<f64>>();

        return MathUtils::normalize_sum_to_one(kernel);
    }

    fn determine_filter_size(kernel: &Vec<f64>, min_prob_to_keep_in_filter: f64) -> usize {
        let middle = (kernel.len() - 1) / 2;
        let mut filter_end = middle;

        while filter_end > 0 {
            if kernel[filter_end - 1] < min_prob_to_keep_in_filter {
                break;
            }
            filter_end -= 1
        }

        return middle - filter_end;
    }

    /**
     * Get the size (in bp) of the band pass filter
     * @return a positive integer
     */
    pub fn get_band_size(&self) -> usize {
        2 * self.filter_size + 1
    }

    /**
     * Get the filter size (which is the size of each wing of the band, minus the center point)
     * @return a positive integer
     */
    pub fn get_filtered_size(&self) -> usize {
        self.filter_size
    }

    /**
     * Get the Gaussian kernel sigma value
     * @return a positive double
     */
    pub fn get_sigma(&self) -> f64 {
        self.sigma
    }

    /**
     * Get the kernel of this band pass filter.  Do not modify returned result
     * @return the kernel used in this band pass filter
     */
    pub fn get_kernel(&self) -> &Vec<f64> {
        &self.gaussian_kernel
    }
}

impl Profile for BandPassActivityProfile {
    /**
     * Our maximize propagation distance is whatever our parent's is, plus our filter size
     *
     * Stops the profile from interpreting sites that aren't yet fully determined due to
     * propagation of the probabilities.
     *
     * @return the distance in bp we might move our probabilities around for some site i
     */
    fn get_max_prob_propagation_distance(&self) -> usize {
        self.activity_profile.get_max_prob_propagation_distance() + self.filter_size
    }

    /**
     * How many profile results are in this profile?
     * @return the number of profile results
     */
    fn size(&self) -> usize {
        self.activity_profile.get_state_list().len()
    }

    /**
     * Is this profile empty? (ie., does it contain no ActivityProfileStates?)
     * @return true if the profile is empty (ie., contains no ActivityProfileStates)
     */
    fn is_empty(&self) -> bool {
        self.activity_profile.get_state_list().is_empty()
    }

    /**
     * Get the span of this activity profile, which is from the start of the first state to the stop of the last
     * @return a potentially null SimpleInterval.  Will be null if this profile is empty
     */
    fn get_span(&self) -> Option<SimpleInterval> {
        if self.is_empty() {
            None
        } else {
            Some(
                self.activity_profile
                    .region_start_loc
                    .as_ref()
                    .unwrap()
                    .span_with(&self.activity_profile.region_stop_loc.as_ref().unwrap()),
            )
        }
    }

    fn get_contig(&self) -> usize {
        self.activity_profile
            .region_start_loc
            .as_ref()
            .unwrap()
            .get_contig()
    }

    fn get_end(&self) -> usize {
        self.activity_profile
            .region_stop_loc
            .as_ref()
            .unwrap()
            .get_end()
    }

    /**
     * Get the list of activity profile results in this object
     * @return a non-null, ordered list of activity profile results
     */
    fn get_state_list(&self) -> &Vec<ActivityProfileState> {
        &self.activity_profile.get_state_list()
    }

    // --------------------------------------------------------------------------------
    //
    // routines to add states to a profile
    //
    // --------------------------------------------------------------------------------

    /**
     * Add the next ActivityProfileState to this profile.
     *
     * Must be contiguous with the previously added result, or an IllegalArgumentException will be thrown
     *
     * @param state a well-formed ActivityProfileState result to incorporate into this profile
     */
    fn add(&mut self, state: ActivityProfileState) {
        let loc = state.get_loc();

        if self.is_empty() {
            self.activity_profile.region_start_loc = Some(loc.clone());
            self.activity_profile.region_stop_loc = Some(loc.clone());
        } else {
            if self
                .activity_profile
                .region_stop_loc
                .as_ref()
                .unwrap()
                .get_start()
                != loc.get_start() - 1
            {
                panic!(
                    "Bad add call to ActivityProfile: loc {:?} not immediately after last loc {:?}",
                    loc, self.activity_profile.region_stop_loc
                )
            }
            self.activity_profile.region_stop_loc = Some(loc.clone());
        }
        let processed_states = self.process_state(state);

        for processed_state in processed_states.into_iter() {
            self.incorporate_single_state(processed_state)
        }
    }

    /**
     * Band pass the probabilities in the ActivityProfile, producing a new profile that's band pass filtered
     * @return a new double[] that's the band-pass filtered version of this profile
     */
    fn process_state(&self, just_added_state: ActivityProfileState) -> Vec<ActivityProfileState> {
        let mut states = Vec::new();

        for super_state in self
            .activity_profile
            .process_state(just_added_state.clone())
        {
            if super_state.is_active_prob() > 0.0 {
                for i in (-(self.filter_size as i64)..=(self.filter_size as i64)).into_iter() {
                    let loc = self.get_loc_for_offset(&just_added_state.get_loc(), i);
                    match loc {
                        Some(loc) => {
                            let new_prob = super_state.is_active_prob()
                                * self.gaussian_kernel[(i + self.filter_size as i64) as usize];
                            states.push(ActivityProfileState::new(loc, new_prob, Type::None))
                        }
                        None => continue,
                    }
                }
            } else {
                states.push(just_added_state.clone())
            }
        }

        return states;
    }

    /**
     * Incorporate a single activity profile state into the current list of states
     *
     * If state's position occurs immediately after the last position in this profile, then
     * the state is appended to the state list.  If it's within the existing states list,
     * the prob of stateToAdd is added to its corresponding state in the list.  If the
     * position would be before the start of this profile, stateToAdd is simply ignored.
     *
     * @param stateToAdd the state we want to add to the states list
     */
    fn incorporate_single_state(&mut self, state_to_add: ActivityProfileState) {
        self.activity_profile.incorporate_single_state(state_to_add)
    }

    /**
     * Get the probabilities of the states as a single linear array of doubles
     * @return a non-null array
     */
    fn get_loc_for_offset(
        &self,
        relative_loc: &SimpleInterval,
        offset: i64,
    ) -> Option<SimpleInterval> {
        self.activity_profile
            .get_loc_for_offset(relative_loc, offset)
    }

    /**
     * Get the length of the current contig
     * @return the length in bp
     */
    fn get_current_contig_length(&self) -> usize {
        self.activity_profile.get_current_contig_length()
    }

    fn pop_ready_assembly_regions(
        &mut self,
        assembly_region_extension: usize,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Vec<AssemblyRegion> {
        self.activity_profile.pop_ready_assembly_regions(
            assembly_region_extension,
            min_region_size,
            max_region_size,
            force_conversion,
        )
    }

    fn pop_next_ready_assembly_region(
        &mut self,
        assembly_region_extension: usize,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Option<AssemblyRegion> {
        self.activity_profile.pop_next_ready_assembly_region(
            assembly_region_extension,
            min_region_size,
            max_region_size,
            force_conversion,
        )
    }

    fn find_end_of_region(
        &mut self,
        is_active_region: bool,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Option<usize> {
        self.activity_profile.find_end_of_region(
            is_active_region,
            min_region_size,
            max_region_size,
            force_conversion,
        )
    }

    fn find_best_cut_site(&self, end_of_active_region: usize, min_region_size: usize) -> usize {
        self.activity_profile
            .find_best_cut_site(end_of_active_region, min_region_size)
    }

    fn find_first_activity_boundary(
        &self,
        is_active_region: bool,
        max_region_size: usize,
    ) -> usize {
        self.activity_profile
            .find_first_activity_boundary(is_active_region, max_region_size)
    }

    fn get_prob(&self, index: usize) -> f64 {
        self.activity_profile.get_prob(index)
    }

    fn is_minimum(&self, index: usize) -> bool {
        self.activity_profile.is_minimum(index)
    }

    fn get_probabilities_as_array(&self) -> Vec<f64> {
        self.activity_profile.get_probabilities_as_array()
    }
}

/**
* Implement the extend method for ActivityProfile when
* given a parallel iterator of ActivityProfileState
*/
impl Extend<ActivityProfileState> for BandPassActivityProfile {
    fn extend<I>(&mut self, iter: I)
    where
        I: IntoIterator<Item = ActivityProfileState>,
    {
        let iter = iter.into_iter();
        iter.for_each(|state| self.add(state));
    }
}
