function Y = buildY(C)
% Build Ybus using series Z, per-end shunt C.Bend (j*B/2), and real off-nominal taps C.tap.
    N = max(C.pairs(:));
    Y = zeros(N,N);
    for e = 1:size(C.pairs,1)
        i = C.pairs(e,1); j = C.pairs(e,2);
        z = C.Zser(e); y = 1/z;
        bEnd = C.Bend(e);                 % j*B/2 at each end
        t = 1;
        if isfield(C,'tap') && ~isempty(C.tap)
            t = C.tap(e); if t==0, t=1; end
        end

        % add line charging (per end)
        Y(i,i) = Y(i,i) + bEnd;
        Y(j,j) = Y(j,j) + bEnd;

        % off-nominal transformer model (tap on i-side)
        Y(i,i) = Y(i,i) + y/(t^2);
        Y(j,j) = Y(j,j) + y;
        Y(i,j) = Y(i,j) - y/t;
        Y(j,i) = Y(j,i) - y/t;
    end
end
