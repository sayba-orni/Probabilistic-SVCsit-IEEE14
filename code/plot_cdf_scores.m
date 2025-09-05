function plot_cdf_scores(lambda, alpha, Nmc, seed)
% PLOT_CDF_SCORES  Empirical CDF of scenario scores J for:
%   - Baseline (no SVC)
%   - Best candidate BUS-9
%
% J = lossMW + lambda * phiV
%
% Usage:
%   plot_cdf_scores                % uses defaults (lambda=100, alpha=0.90, Nmc=500, seed=42)
%   plot_cdf_scores(150,0.95,1000,123)

    if nargin < 1 || isempty(lambda), lambda = 100; end
    if nargin < 2 || isempty(alpha),  alpha  = 0.90; end
    if nargin < 3 || isempty(Nmc),    Nmc    = 500;  end
    if nargin < 4 || isempty(seed),   seed   = 42;   end

    clc; fprintf('>>> plot_cdf_scores: start (lambda=%g, alpha=%.2f, Nmc=%d, seed=%d)\n', ...
                 lambda, alpha, Nmc, seed);

    % setup & paths
    here     = fileparts(mfilename('fullpath'));
    figDir   = fullfile(here, 'figures');
    if ~exist(figDir,'dir'), mkdir(figDir); end

    %  rebuild same scenarios for a fair head-to-head
    rng(seed,'twister');
    C = config14();
    S = build_scenarios(C, Nmc, 'p_outage', 0.05, 'timeMode', "random");

    %  baseline: per-scenario J
    [J_base, ok_base] = compute_J_plain(C, S, lambda);
    J_base = J_base(ok_base);

    % BUS-9 candidate: per-scenario J 
    cand = struct('type',"bus",'k',9,'a',[],'b',[]);
    R = evaluate_candidate(cand, C, S);
    ok_cand = R.converged & isfinite(R.lossMW) & isfinite(R.phiV);
    J_bus9  = R.lossMW(ok_cand) + lambda * R.phiV(ok_cand);

    %  stats
    base.mean  = mean(J_base);
    base.cvar  = cvar_local(J_base, alpha);
    cand.mean  = mean(J_bus9);
    cand.cvar  = cvar_local(J_bus9, alpha);

    q_base = quantile(J_base, alpha);
    q_cand = quantile(J_bus9,  alpha);

    fprintf('Baseline:   n=%d, mean=%.3f, CVaR=%.3f (q_%.2f=%.3f)\n', numel(J_base), base.mean, base.cvar, alpha, q_base);
    fprintf('BUS-9:      n=%d, mean=%.3f, CVaR=%.3f (q_%.2f=%.3f)\n', numel(J_bus9), cand.mean, cand.cvar, alpha, q_cand);

    %  plot empirical CDFs (no toolbox dependency)
    [xB, FB] = ecdf_manual(J_base);
    [x9, F9] = ecdf_manual(J_bus9);

    f = figure('Color','w','Name','CDF of scenario scores (Baseline vs BUS-9)');
    plot(xB, FB, '-', 'LineWidth', 1.8); hold on; grid on;
    plot(x9, F9, '-', 'LineWidth', 1.8);

    % mark alpha-quantiles (VaR) and CVaR means
    yl = ylim;
    plot([q_base q_base], yl, ':', 'LineWidth', 1.0);
    plot([q_cand q_cand], yl, ':', 'LineWidth', 1.0);

    % small markers at CVaR (draw as points on top axis)
    plot(base.cvar, 1.0, 'o', 'MarkerSize', 6, 'LineWidth', 1.2);
    plot(cand.cvar, 1.0, 'o', 'MarkerSize', 6, 'LineWidth', 1.2);

    legend({ ...
        sprintf('Baseline  (n=%d)', numel(J_base)), ...
        sprintf('BUS-9     (n=%d)', numel(J_bus9)), ...
        sprintf('Baseline VaR_{%.2f}=%.2f', alpha, q_base), ...
        sprintf('BUS-9 VaR_{%.2f}=%.2f',    alpha, q_cand), ...
        sprintf('Baseline CVaR=%.2f', base.cvar), ...
        sprintf('BUS-9  CVaR=%.2f',   cand.cvar) }, ...
        'Location','southeast');

    xlabel(sprintf('Scenario score  J = Loss + \\lambda\\cdotVoltPen   (\\lambda = %g)', lambda));
    ylabel('Empirical CDF  F_J(x)');
    title(sprintf('CDF of scenario scores (N=%d, seed=%d) — \\alpha=%.2f', Nmc, seed, alpha));

    % save
    outPNG = fullfile(figDir, sprintf('cdf_scores_lambda%d_alpha%02d_seed%d.png', ...
                    round(lambda), round(100*alpha), seed));
    exportgraphics(f, outPNG, 'Resolution', 200);
    fprintf('>>> wrote %s\n', outPNG);
    fprintf('>>> plot_cdf_scores: done\n');
end



function [J, ok] = compute_J_plain(C, S, lambda)
    % Per-scenario baseline score J using solver_plain
    Y0    = buildY(C);
    Sbase = C.Sbase;
    Nsc   = numel(S);
    J     = nan(Nsc,1);
    ok    = false(Nsc,1);

    for s = 1:Nsc
        % outage-aware network
        if S(s).outage == "none"
            Y = Y0; pairs_s = C.pairs; Zser_s = C.Zser;
        else
            ab = sscanf(S(s).outage,'%d-%d');
            [Y, pairs_s, Zser_s] = remove_line_general(Y0, C, ab(1), ab(2));
        end
        if is_islanded(pairs_s, max(C.pairs(:))), continue; end

        % loads/gens (pu)
        Pd = (C.Pd_MW   .* (1 + S(s).epsLoad)) / Sbase;
        Qd = (C.Qd_MVAr .* (1 + S(s).epsLoad)) / Sbase;
        Pg =  C.Pg_MW   / Sbase;
        Qg =  C.Qg_MVAr / Sbase;

        % renewables (pu)
        Ppv_MW   = pvPlantMW_datasheet(S(s).G,  S(s).TA,  C.PV);
        Pwind_MW = windPlantMW_formula(S(s).Vw, C.WF);
        if isfield(C,'bus_solar') && ~isempty(C.bus_solar)
            Pg(C.bus_solar) = Pg(C.bus_solar) + Ppv_MW / Sbase;
        end
        if isfield(C,'bus_wind') && ~isempty(C.bus_wind)
            Pg(C.bus_wind)  = Pg(C.bus_wind)  + Pwind_MW / Sbase;
        end

        % solve
        res = solver_plain(Y, pairs_s, Zser_s, Sbase, Pg, Qg, Pd, Qd);
        if ~(isfield(res,'converged') && res.converged), continue; end

        Vabs = res.V_abs(:);
        if any(~isfinite(Vabs)) || any(Vabs < 0.8 | Vabs > 1.2), continue; end

        [phiV,~] = volt_penalty(Vabs(2:end), C.Vband);
        J(s)     = res.totPloss + lambda * phiV;
        ok(s)    = isfinite(J(s));
    end
end

function [x, F] = ecdf_manual(v)
    v = sort(v(:));
    n = numel(v);
    x = v;
    F = (1:n)'/n;
end

function c = cvar_local(x, alpha)
    x  = sort(x(:),'ascend');
    k  = max(1, ceil(alpha*numel(x)));
    c  = mean(x(k:end));
end
