function seq_out = S3_coveragePreserved_durationRandomized(seq_in, states, seed)
% S3_COVERAGEPRESERVED_DURATIONRANDOMIZED
% Coverage-preserved duration-randomized surrogate for microstate sequences.
%
%   seq_out = S3_coveragePreserved_durationRandomized(seq_in)
%   seq_out = S3_coveragePreserved_durationRandomized(seq_in, states)
%   seq_out = S3_coveragePreserved_durationRandomized(seq_in, states, seed)
%
%   This surrogate preserves the original symbolic segment order and the
%   original coverage of each microstate, while randomizing the durations of
%   individual segments within each state.
%
%   Preserves:
%       - symbolic segment order
%       - segment labels/order of runs
%       - number of segments per microstate
%       - total coverage of each microstate
%       - total sequence length
%
%   Randomizes:
%       - individual segment durations within each microstate
%
%   Inputs:
%       seq_in : vector of integer microstate labels
%       states : optional vector of allowed states. Default: unique(seq_in).
%       seed   : optional random seed for reproducibility
%
%   Output:
%       seq_out: surrogate sequence with the same length as seq_in and the
%                same per-state coverage as seq_in

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

    for k = 1:numel(states)
        st = states(k);
        idx = find(labels == st);
        nk = numel(idx);
        Tk = sum(seq == st);

        if nk == 0
            if Tk == 0
                continue;
            else
                error('State %g has coverage but no segments.', st);
            end
        end
        if Tk < nk
            error('State %g coverage (%d samples) is smaller than its number of segments (%d).', st, Tk, nk);
        end
        new_lengths(idx) = local_random_positive_composition(Tk, nk);
    end

    out = local_rebuild_sequence(labels, new_lengths);
    if numel(out) ~= T
        error('Internal error: output length does not match input length.');
    end

    for k = 1:numel(states)
        st = states(k);
        if sum(out == st) ~= sum(seq == st)
            error('Internal error: coverage mismatch for state %g.', st);
        end
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
        bins = randi(nParts, 1, remaining);
        for b = bins
            parts(b) = parts(b) + 1;
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
