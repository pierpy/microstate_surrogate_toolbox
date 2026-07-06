function seq_out = S2_coverageEqualized(seq_in, states, seed)
% S2_COVERAGEEQUALIZED  Coverage-equalized surrogate for microstate sequences.
%
%   seq_out = S2_coverageEqualized(seq_in)
%   seq_out = S2_coverageEqualized(seq_in, states)
%   seq_out = S2_coverageEqualized(seq_in, states, seed)
%
%   This surrogate preserves the original symbolic segment order while
%   adjusting microstate occupancies to obtain approximately equal coverage
%   across states. Segment durations are modified accordingly.
%
%   Preserves:
%       - symbolic segment order
%       - segment labels/order of runs
%       - number of segments per microstate
%       - total sequence length
%
%   Modifies:
%       - microstate coverage, which is made approximately equal across states
%       - segment durations, as required to satisfy the coverage constraint
%
%   Notes:
%       If T is the sequence length and K is the number of states, each state
%       is assigned floor(T/K) samples. The remaining samples are assigned to
%       randomly selected states, so coverage is exactly equal only when T is
%       divisible by K.
%
%   Inputs:
%       seq_in : vector of integer microstate labels
%       states : optional vector of allowed states. Default: unique(seq_in).
%                Use this to enforce a fixed state set, e.g. 1:7.
%       seed   : optional random seed for reproducibility
%
%   Output:
%       seq_out: surrogate sequence with the same length as seq_in

    validateattributes(seq_in, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'seq_in', 1);
    wasRow = isrow(seq_in);
    seq = seq_in(:).';
    T = numel(seq);

    if nargin < 2 || isempty(states)
        states = unique(seq, 'stable');
    else
        states = states(:).';
    end
    if nargin >= 3 && ~isempty(seed)
        rng(seed);
    end

    [labels, ~] = local_run_length_encoding(seq);
    nSeg = numel(labels);
    new_lengths = zeros(1, nSeg);

    base = floor(T / numel(states));
    target = base * ones(1, numel(states));
    rem = T - sum(target);
    if rem > 0
        p = randperm(numel(states));
        target(p(1:rem)) = target(p(1:rem)) + 1;
    end

    for k = 1:numel(states)
        st = states(k);
        idx = find(labels == st);
        nk = numel(idx);
        if nk == 0
            error('State %g is not present in the sequence; coverage cannot be equalized for absent states.', st);
        end
        if target(k) < nk
            error('State %g target coverage (%d samples) is smaller than its number of segments (%d).', st, target(k), nk);
        end
        new_lengths(idx) = local_random_positive_composition(target(k), nk);
    end

    out = local_rebuild_sequence(labels, new_lengths);
    if numel(out) ~= T
        error('Internal error: output length does not match input length.');
    end

    seq_out = out;
    if ~wasRow
        seq_out = seq_out(:);
    end
end

function [labels, lengths] = local_run_length_encoding(seq)
    changeIdx = [1, find(diff(seq) ~= 0) + 1, numel(seq) + 1];
    starts = changeIdx(1:end-1);
    ends = changeIdx(2:end) - 1;
    labels = seq(starts);
    lengths = ends - starts + 1;
end

function parts = local_random_positive_composition(total, nParts)
    parts = ones(1, nParts);
    remaining = total - nParts;
    if remaining > 0
        cuts = randi(nParts, 1, remaining);
        for c = cuts
            parts(c) = parts(c) + 1;
        end
    end
end

function out = local_rebuild_sequence(labels, lengths)
    out = zeros(1, sum(lengths));
    pos = 1;
    for i = 1:numel(labels)
        L = lengths(i);
        out(pos:pos+L-1) = labels(i);
        pos = pos + L;
    end
end
