use rand::prelude::*;

lazy_static! {
    static ref code_to_base: Vec<u8> = vec![b'A', b'C', b'G', b'T'];
    static ref DNA_SET: &'static [u8] = b"ATCG";
}

/**
 * Random DNA sequence generator.
 *
 * <p>
 *     Returned bases are always in upper case and one of the valid four nucleotides 'A', 'C', 'G' and 'T'.
 * </p>
 *
 * <p>
 *     This class is written in a way that ensure that the sequence of bases returned would only depend on the seed and
 *     not the sequence of calls to different base generating methods.
 * </p>
 * <p>
 *     For example a single call asking for 132 bases using {@link #nextBases(int)} would generate the same sequence
 *     as calling 132 times {@link #next_base()} if both random instances were initialized with the same seed.
 * </p>
 * Originally written as part of Broad Institute GATK tool kit and translated to Rust
 *
 * @author Rhys Newell &lt;rhys.newell@hdr.qut.edu.au&gt;
 */
pub struct RandomDNA {
    random: ThreadRng,
}

impl RandomDNA {
    /**
     * Constructs a new random DNA generator.
     * <p>
     * <p>
     * The seed would be the default which would depend on system properties and the current time as
     * described in {@link Random} documentation.
     * </p>
     */
    /**
     * Constructs a new random DNA generator.
     * <p>
     * <p>
     * The seed would be the default which would depend on system properties and the current time as
     * described in {@link Random} documentation.
     * </p>
     */
    pub fn new(rnd: ThreadRng) -> Self {
        Self { random: rnd }
    }

    /**
     * Puts a random sequence of base into an byte array.
     *
     * @param length how many based to put into the array.
     * @throws IllegalArgumentException if {@code dest} is {@code null} or the pair {@code offset} and {@code length}
     *                                  do not specify a valid index sub-interval in {@code dest}.
     */
    pub fn next_bases(&mut self, length: usize) -> String {
        let sequence: String = (0..length)
            .map(|_| {
                let idx = self.random.gen_range(0, DNA_SET.len());
                DNA_SET[idx] as char
            })
            .collect();

        return sequence;
    }

    // /**
    //  * Load new bases to the next bases queue to try to satisfy a required number of next-bases.
    //  *
    //  * <p>
    //  *     It might be that there is not enough space in the buffer to satisfy the target, in which case this will return
    //  * something close to the capacity of the buffer, expecting the caller to clear the buffer and repeat the request
    //  * for the remainder.
    //  * </p>
    //  *
    //  * @param targetSize the number of bases that the caller would like to have at its disposal.
    //  * @return the amount of bases available given the current content of the next-bases queue, its max capacity
    //  * and the requested target.
    //  */
    // fn load_next_bases(&mut self, target_size: usize) -> usize {
    //     let current_number_of_bases = self.next_bases_size;
    //
    //     if current_number_of_bases >= target_size {
    //         return target_size
    //     } else if NEXT_BASES_MAX_CAPACITY - current_number_of_bases < BASES_IN_AN_INT {
    //         // If we try to load more bases it would go over the max capacity.
    //         return current_number_of_bases
    //     } else {
    //         let result = min(NEXT_BASES_MAX_CAPACITY, target_size);
    //         let increase_in_ints = (result - current_number_of_bases + BASES_IN_AN_INT_MINUS_1) / BASES_IN_AN_INT;
    //         for i in 0..increase_in_ints {
    //             let random_int = self.random.gen();
    //         }
    //     }
    // }
}
