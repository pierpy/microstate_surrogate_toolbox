function seq_out = generate_microstate_surrogate(seq_in, surrogateType, varargin)
% GENERATE_MICROSTATE_SURROGATE  General wrapper for S1-S7 microstate surrogates.
%
%   seq_out = generate_microstate_surrogate(seq_in, surrogateType)
%   seq_out = generate_microstate_surrogate(..., 'Name', Value)
%
%   surrogateType can be:
%       'S1' : duration equalization
%       'S2' : coverage equalization
%       'S3' : coverage-preserved duration randomization
%       'S4' : segment-order shuffling
%       'S5' : S1 followed by S4
%       'S6' : S2 followed by S4
%       'S7' : S3 followed by S4
%
%   Name-value options:
%       'Lm'            : common segment duration for S1/S5. Default: median duration.
%       'states'        : state labels for S2/S3/S6/S7. Default: unique(seq_in).
%       'seed'          : random seed. Default: [].
%       'avoidAdjacent' : logical flag for S4/S5/S6/S7. Default: false.
%
%   Examples:
%       seqS3 = generate_microstate_surrogate(seq, 'S3', 'states', 1:7, 'seed', 1);
%       seqS5 = generate_microstate_surrogate(seq, 'S5', 'Lm', 10, 'seed', 1);

    p = inputParser;
    addRequired(p, 'seq_in', @(x) isnumeric(x) && isvector(x) && ~isempty(x));
    addRequired(p, 'surrogateType', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Lm', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(p, 'states', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
    addParameter(p, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'avoidAdjacent', false, @(x) islogical(x) || isnumeric(x));
    parse(p, seq_in, surrogateType, varargin{:});

    Lm = p.Results.Lm;
    states = p.Results.states;
    seed = p.Results.seed;
    avoidAdjacent = logical(p.Results.avoidAdjacent);
    typ = upper(char(p.Results.surrogateType));

    switch typ
        case 'S1'
            seq_out = S1_durationEqualized(seq_in, Lm);

        case 'S2'
            seq_out = S2_coverageEqualized(seq_in, states, seed);

        case 'S3'
            seq_out = S3_coveragePreserved_durationRandomized(seq_in, states, seed);

        case 'S4'
            seq_out = S4_orderShuffled(seq_in, seed, avoidAdjacent);

        case 'S5'
            seq_tmp = S1_durationEqualized(seq_in, Lm);
            seq_out = S4_orderShuffled(seq_tmp, seed, avoidAdjacent);

        case 'S6'
            seq_tmp = S2_coverageEqualized(seq_in, states, seed);
            seq_out = S4_orderShuffled(seq_tmp, seed, avoidAdjacent);

        case 'S7'
            seq_tmp = S3_coveragePreserved_durationRandomized(seq_in, states, seed);
            seq_out = S4_orderShuffled(seq_tmp, seed, avoidAdjacent);

        otherwise
            error('Unknown surrogate type: %s. Use S1, S2, S3, S4, S5, S6, or S7.', typ);
    end
end
