function seq_out = S1_durationEqualized(seq_in, Lm)
% S1_DURATIONEQUALIZED  Duration-equalized surrogate for microstate sequences.
%
%   seq_out = S1_durationEqualized(seq_in)
%   seq_out = S1_durationEqualized(seq_in, Lm)
%
%   This surrogate preserves the symbolic order of microstate segments/runs,
%   while setting all segment durations to a common value Lm. If Lm is not
%   provided, it is set to the median segment duration of seq_in, rounded to
%   the nearest integer.
%
%   Preserves:
%       - symbolic segment order
%       - segment labels
%       - total sequence length
%
%   Modifies:
%       - segment durations, which are equalized to Lm samples
%       - microstate coverage may change as a consequence of duration equalization
%
%   Notes:
%       If the reconstructed sequence is longer than the original sequence,
%       it is trimmed to the original length. If it is shorter, the last
%       segment is extended to match the original length.
%
%   Inputs:
%       seq_in : vector of integer microstate labels
%       Lm     : optional positive integer segment duration in samples
%
%   Output:
%       seq_out: surrogate sequence with the same length as seq_in

    validateattributes(seq_in, {'numeric'}, {'vector', 'nonempty'}, mfilename, 'seq_in', 1);
    wasRow = isrow(seq_in);
    seq = seq_in(:);
    N = numel(seq);

    [labels, lengths] = local_run_length_encoding(seq);

    if nargin < 2 || isempty(Lm)
        Lm = round(median(lengths));
    end
    validateattributes(Lm, {'numeric'}, {'scalar', 'integer', 'positive'}, mfilename, 'Lm', 2);

    out = zeros(0, 1);
    for i = 1:numel(labels)
        out = [out; repmat(labels(i), Lm, 1)]; %#ok<AGROW>
    end

    if numel(out) > N
        out = out(1:N);
    elseif numel(out) < N
        out = [out; repmat(labels(end), N - numel(out), 1)]; %#ok<AGROW>
    end

    seq_out = out;
    if wasRow
        seq_out = seq_out.';
    end
end

function [labels, lengths] = local_run_length_encoding(seq)
    changeIdx = [1; find(diff(seq) ~= 0) + 1; numel(seq) + 1];
    starts = changeIdx(1:end-1);
    ends = changeIdx(2:end) - 1;
    labels = seq(starts);
    lengths = ends - starts + 1;
end
