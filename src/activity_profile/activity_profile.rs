use activity_profile::activity_profile_state::ActivityProfileState;
use utils::simple_interval::SimpleInterval;

pub struct ActivityProfile {
    state_list: Vec<ActivityProfileState>,
    max_prob_propagation_distance: usize,
    active_prob_threshold: f64,
    region_start_loc: SimpleInterval,
    region_stop_loc: SimpleInterval,
    contig_len: usize,
}

impl ActivityProfile {

    /**
     * Create a empty ActivityProfile, restricting output to profiles overlapping intervals, if not null
     * @param maxProbPropagationDistance region probability propagation distance beyond its maximum size
     * @param activeProbThreshold threshold for the probability of a profile state being active
     */
    pub fn new(max_prob_propagation_distance: usize, active_prob_threshold: f64) -> ActivityProfile {
        ActivityProfile {
            state_list: Vec::new(),
            max_prob_propagation_distance,
            active_prob_threshold,
            region_start_loc: SimpleInterval::new(0, 0, 0),
            region_stop_loc: SimpleInterval::new(0, 0, 0),
            contig_len: 0,
        }
    }

    /**
     * How far away can probability mass be moved around in this profile?
     *
     * This distance puts an upper limit on how far, in bp, we will ever propagate probability mass around
     * when adding a new ActivityProfileState.  For example, if the value of this function is
     * 10, and you are looking at a state at bp 5, and we know that no states beyond 5 + 10 will have
     * their probability propagated back to that state.
     *
     * @return a positive integer distance in bp
     */
    pub fn get_max_prob_propagation_distance(&self) -> usize {
        self.max_prob_propagation_distance
    }
}