use activity_profile::activity_profile::ActivityProfile;

/**
 * A band pass filtering version of the activity profile
 *
 * Applies a band pass filter with a Gaussian kernel to the input state probabilities to smooth
 * them out of an interval
 */
pub struct BandPassActivityProfile {
    filter_size: usize,
    sigma: f64,
    gaussian_kernel: Vec<f64>,
    activity_profile: ActivityProfile
}

impl BandPassActivityProfile {
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
        tid: usize,
        contig_len: usize,
    ) -> BandPassActivityProfile {
        let activity_profile = ActivityProfile::new(max_prob_propagation_distance, active_prob_threshold, ref_idx, tid, contig_len);

        BandPassActivityProfile {
            activity_profile,
            sigma,

        }
    }

    fn make_kernel(filter_size: usize, sigma: f64) -> Vec<f64> {
        let band_size = 2 * filter_size + 1;
        let mut kernel = vec![0.0; band_size];

    }
}