function seq_out = S4_orderShuffled(seq_in, seed, avoidAdjacent)
% S4_ORDERSHUFFLED  Segment-order shuffled surrogate for microstate sequences.
%
%   seq_out = S4_orderShuffled(seq_in)
%   seq_out = S4_orderShuffled(seq_in, seed)
%   seq_out = S4_orderShuffled(seq_in, seed, avoidAdjacent)
%
%   This surrogate treats each microstate segment as an intact label-duration
%   unit and randomly permutes the order of these units.
%
%   Preserves:
%       - segment labels
%       - individual segment durations
%       - number of segments before reconstruction
%       - total coverage of each microstate
%       - total sequence length
%
%   Randomizes:
%       - symbolic segment order
%       - transition structure / symbolic syntax
%
%   Notes:
%       If avoidAdjacent is true, the function attempts to find a permutation
%       in which no two consecutive segments have the same label. This avoids
%       merging adjacent same-label segments after reconstruction. If no such
%       permutation is found after repeated attempts, the best available random
%       permutation is used and a warning is issued.
%
%   Inputs:
%       seq_in        : vector of integer microstate labels
%       seed          : optional random seed for reproducibility
%       avoidAdjacent : optional logical flag. Default: false
%
%   Output:
%       seq_out: surrogate sequence with the same length and coverage as seq_in

    validateattributes(seq_in, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'seq_in', 1);
    wasRow = isrow(seq_in);
    seq = seq_in(:).';

    if nargin >= 2 && ~isempty(seed)
        rng(seed);
    end
    if nargin < 3 || isempty(avoidAdjacent)
        avoidAdjacent = false;
    end

    [labels, lengths] = local_run_length_encoding(seq);
    nSeg = numel(labels);

    if nSeg <= 1
        seq_out = seq_in;
        return;
    end

    if avoidAdjacent
        maxAttempts = 10000;
        bestPerm = randperm(nSeg);
        bestConflicts = local_count_adjacent_conflicts(labels(bestPerm));
        perm = bestPerm;
        found = bestConflicts == 0;
        attempt = 1;
        while ~found && attempt < maxAttempts
            candidate = randperm(nSeg);
            conflicts = local_count_adjacent_conflicts(labels(candidate));
            if conflicts < bestConflicts
                bestConflicts = conflicts;
                bestPerm = candidate;
            end
            if conflicts == 0
                perm = candidate;
                found = true;
                break;
            end
            attempt = attempt + 1;
        end
        if ~found
            perm = bestPerm;
            warning('S4_orderShuffled:AdjacentLabels', ...
                'Could not avoid all adjacent equal labels; using permutation with %d adjacent conflict(s).', bestConflicts);
        end
    else
        perm = randperm(nSeg);
    end

    out = local_rebuild_sequence(labels(perm), lengths(perm));
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

function n = local_count_adjacent_conflicts(labels)
    n = sum(diff(labels) == 0);
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
