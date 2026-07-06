# Microstate Surrogate Toolbox

MATLAB functions for generating surrogate sequences from symbolic EEG microstate time series. The functions operate on **run-length encoded microstate segments**: each segment is treated as a consecutive run of identical labels with an associated duration.

The toolbox implements four elementary surrogate operations and three composite surrogates:

- **S1**: duration equalization
- **S2**: coverage equalization
- **S3**: coverage-preserved duration randomization
- **S4**: segment-order shuffling
- **S5**: S1 followed by S4
- **S6**: S2 followed by S4
- **S7**: S3 followed by S4

These surrogates are intended for testing which properties of a symbolic microstate sequence contribute to temporal statistics such as long-range temporal correlations, DFA/Hurst exponents, transition structure, or occupancy-related effects.

## Installation

Clone or download this repository and add the folder to your MATLAB path:

```matlab
addpath('microstate_surrogate_toolbox')
```

No MATLAB toolbox dependency is required.

## Input format

Each function expects a numeric vector of microstate labels, for example:

```matlab
seq = [1 1 1 2 2 3 3 3 1 1 4 4];
```

Both row and column vectors are accepted. The output preserves the input orientation.

## Quick start

```matlab
seq = [1 1 1 2 2 3 3 3 1 1 4 4 2 2 2 3 3];

S1 = generate_microstate_surrogate(seq, 'S1', 'Lm', 10);
S2 = generate_microstate_surrogate(seq, 'S2', 'states', 1:4, 'seed', 1);
S3 = generate_microstate_surrogate(seq, 'S3', 'states', 1:4, 'seed', 1);
S4 = generate_microstate_surrogate(seq, 'S4', 'seed', 1);
S5 = generate_microstate_surrogate(seq, 'S5', 'Lm', 10, 'seed', 1);
S6 = generate_microstate_surrogate(seq, 'S6', 'states', 1:4, 'seed', 1);
S7 = generate_microstate_surrogate(seq, 'S7', 'states', 1:4, 'seed', 1);
```

You can also call each elementary function directly:

```matlab
S1 = S1_durationEqualized(seq, 10);
S2 = S2_coverageEqualized(seq, 1:4, 1);
S3 = S3_coveragePreserved_durationRandomized(seq, 1:4, 1);
S4 = S4_orderShuffled(seq, 1, false);
```

## Surrogate definitions

| Surrogate | Function | Preserved | Modified / randomized | Purpose |
|---|---|---|---|---|
| Original | — | Segment order, durations, coverage, transition structure | — | Reference sequence |
| S1 | `S1_durationEqualized` | Symbolic segment order, segment labels, total sequence length | Segment durations are equalized to `Lm`; coverage may change | Tests whether heterogeneous dwell times contribute to the observed temporal structure |
| S2 | `S2_coverageEqualized` | Symbolic segment order, number of segments per state, total sequence length | Coverage is made approximately equal across states; durations are adjusted accordingly | Tests whether unequal state occupancy contributes to the observed temporal structure |
| S3 | `S3_coveragePreserved_durationRandomized` | Symbolic segment order, number of segments per state, total coverage per state, total sequence length | Durations of individual segments are randomized within each state | Tests whether the original dwell-time organization contributes beyond coverage |
| S4 | `S4_orderShuffled` | Individual segment labels and durations, total coverage per state, total sequence length | Order of complete label-duration segments is shuffled | Tests whether symbolic syntax / transition order contributes beyond dwell-time statistics |
| S5 | `generate_microstate_surrogate(...,'S5')` | Total sequence length | S1 followed by S4 | Removes both heterogeneous dwell times and original symbolic order |
| S6 | `generate_microstate_surrogate(...,'S6')` | Total sequence length; approximately equal coverage from S2 | S2 followed by S4 | Removes coverage imbalance and original symbolic order |
| S7 | `generate_microstate_surrogate(...,'S7')` | Total sequence length and original coverage from S3 | S3 followed by S4 | Preserves coverage while disrupting both original dwell-time organization and symbolic order |

## Function details

### S1: duration equalization

```matlab
seq_out = S1_durationEqualized(seq_in, Lm)
```

S1 preserves the symbolic segment order and sets all segment durations to a common duration `Lm` samples. If `Lm` is omitted, the median segment duration of the input sequence is used.

If the reconstructed sequence does not exactly match the original length, it is trimmed or the final segment is extended to preserve the original total length.

### S2: coverage equalization

```matlab
seq_out = S2_coverageEqualized(seq_in, states, seed)
```

S2 preserves the original symbolic segment order but adjusts the total occupancy of each microstate to obtain approximately equal coverage across states. If the sequence length `T` is not divisible by the number of states `K`, the remaining samples are assigned randomly to states.

Use `states` to enforce the complete state set, for example `1:7`. This is recommended for group analyses, because relying on `unique(seq_in)` can silently change `K` if a state is absent from a specific sequence.

### S3: coverage-preserved duration randomization

```matlab
seq_out = S3_coveragePreserved_durationRandomized(seq_in, states, seed)
```

S3 preserves the symbolic segment order and the original total coverage of each state, but randomizes individual segment durations within each state.

This surrogate disrupts the original dwell-time sequence without changing the overall occupancy of each microstate.

### S4: segment-order shuffling

```matlab
seq_out = S4_orderShuffled(seq_in, seed, avoidAdjacent)
```

S4 treats each segment as an intact label-duration unit and randomly permutes the order of these units. Thus, individual segment durations and total state coverage are preserved, while symbolic order and transition structure are disrupted.

The optional flag `avoidAdjacent` attempts to avoid placing two segments with the same label next to each other. This prevents post-shuffle merging of adjacent same-label segments. Because such a permutation may not always exist, the function falls back to the best random attempt and issues a warning.

## Composite surrogates

Composite surrogates are generated by cascading elementary operations:

```matlab
S5 = generate_microstate_surrogate(seq, 'S5', 'Lm', 10, 'seed', 1);
S6 = generate_microstate_surrogate(seq, 'S6', 'states', 1:7, 'seed', 1);
S7 = generate_microstate_surrogate(seq, 'S7', 'states', 1:7, 'seed', 1);
```

Definitions:

- **S5** = S1 duration equalization followed by S4 segment-order shuffling.
- **S6** = S2 coverage equalization followed by S4 segment-order shuffling.
- **S7** = S3 coverage-preserved duration randomization followed by S4 segment-order shuffling.

## Reproducibility

Functions involving randomization accept an optional `seed` argument. For reproducible results, set this argument explicitly:

```matlab
S3 = S3_coveragePreserved_durationRandomized(seq, 1:7, 42);
```

The wrapper also accepts `'seed'` as a name-value argument:

```matlab
S3 = generate_microstate_surrogate(seq, 'S3', 'states', 1:7, 'seed', 42);
```

## Recommended reporting in manuscripts

When using these functions, report:

1. The number of microstate classes `K`.
2. The sequence sampling rate.
3. The common duration `Lm` used for S1/S5.
4. Whether a fixed `states` vector was used.
5. The number of independent surrogate realizations.
6. Whether `avoidAdjacent` was enabled for order-shuffled surrogates.
7. The random seed policy, if exact reproducibility is required.

## Citation

If you use this toolbox in your research, please cite:

Zappasodi F, Croce P.
Scale-free Organization in Brain-state Dynamics via a Multivariater…. (to be updated upon publication)

## License

This project is released under the MIT License.
See the LICENSE file for details.
