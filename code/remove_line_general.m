function [Yout, pairs_s, Zser_s, Bend_s, tap_s] = remove_line_general(Yin, C, a, b)
% Remove the branch (a,b) from a Ybus built with buildY(),
% and return the filtered line lists for this scenario.
% Outputs:
%   Yout    : Ybus after removing the branch
%   pairs_s : remaining line endpoints
%   Zser_s  : remaining series impedances
%   Bend_s  : remaining per-end shunts (j*B/2)
%   tap_s   : remaining off-nominal taps

    % find the line
    idx = find( (C.pairs(:,1)==a & C.pairs(:,2)==b) | ...
                (C.pairs(:,1)==b & C.pairs(:,2)==a), 1);
    if isempty(idx)
        error('remove_line_general: line %d-%d not found', a, b);
    end

    % reconstruct this line's own Y-stamp to subtract it safely
    i = C.pairs(idx,1); j = C.pairs(idx,2);
    z = C.Zser(idx);   y = 1/z;
    bEnd = C.Bend(idx);
    t = 1;
    if isfield(C,'tap') && ~isempty(C.tap)
        t = C.tap(idx); if t==0, t=1; end
    end

    Yout = Yin;

    % subtract the shunt halves
    Yout(i,i) = Yout(i,i) - bEnd;
    Yout(j,j) = Yout(j,j) - bEnd;

    % subtract the series/tap contributions
    Yout(i,i) = Yout(i,i) - y/(t^2);
    Yout(j,j) = Yout(j,j) - y;
    Yout(i,j) = Yout(i,j) + y/t;   % minus( -y/t )
    Yout(j,i) = Yout(j,i) + y/t;

    % filter the line lists
    keep = true(size(C.pairs,1),1);
    keep(idx) = false;

    pairs_s = C.pairs(keep,:);
    Zser_s  = C.Zser(keep,:);
    Bend_s  = C.Bend(keep,:);
    if isfield(C,'tap') && ~isempty(C.tap)
        tap_s = C.tap(keep,:);
    else
        tap_s = ones(sum(keep),1);
    end
end
