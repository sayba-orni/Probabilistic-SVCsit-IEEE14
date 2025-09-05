function T = rank_candidates(C, scenarios, candidates, lambda, alpha)
M = numel(candidates);
meanLoss  = nan(M,1);
meanPhi   = nan(M,1);
violProb  = nan(M,1);
cvarScore = nan(M,1);
meanScore = nan(M,1);

for m = 1:M
    if candidates(m).type == "bus"
        nm = sprintf('BUS-%d', candidates(m).k);
    else
        nm = sprintf('MID-%d-%d', candidates(m).a, candidates(m).b);
    end
    R = evaluate_candidate(candidates(m), C, scenarios);
    valid = R.converged & isfinite(R.lossMW) & isfinite(R.phiV);
    nSc = numel(valid); nBad = nSc - nnz(valid);
    fprintf('%-10s  infeasible/nonconv: %3d / %3d (%.1f%%)\n', nm, nBad, nSc, 100*nBad/max(1,nSc));

    if ~any(valid)
        bigM = 5e3;
        meanLoss(m)=bigM; meanPhi(m)=bigM; violProb(m)=1; cvarScore(m)=bigM; meanScore(m)=bigM; 
        continue;
    end
    loss = R.lossMW(valid); phi = R.phiV(valid); viol = R.vviol(valid);
    J = loss + lambda*phi;

    meanLoss(m)  = mean(loss);
    meanPhi(m)   = mean(phi);
    violProb(m)  = mean(double(viol));
    meanScore(m) = mean(J);
    cvarScore(m) = cvar(J, alpha);
end

name = strings(M,1);
for m=1:M
    if candidates(m).type=="bus"
        name(m) = sprintf('BUS-%d',candidates(m).k);
    else
        name(m) = sprintf('MID-%d-%d',candidates(m).a,candidates(m).b);
    end
end

T = table(name, meanLoss, meanPhi, violProb, meanScore, cvarScore);
T = sortrows(T, 'cvarScore');
end
