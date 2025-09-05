function plot_tail_voltage(lambda, alpha, Nmc, seed, whichQuant)
% PLOT_TAIL_VOLTAGE  Compare bus voltage magnitudes for a tail scenario
% picked from the Baseline (no SVC) distribution and re-solved with BUS-9.
%
% The scenario chosen is the ceil(whichQuant*N)th smallest J for the
% Baseline (no SVC), where J = Loss + lambda * VoltPenalty.
%
% Usage:
%   plot_tail_voltage                       % defaults (100, 0.90, 500, 42, alpha)
%   plot_tail_voltage(150, 0.95, 1000, 7)   % custom (quantile defaults to alpha)

    if nargin < 1 || isempty(lambda),    lambda = 100;  end
    if nargin < 2 || isempty(alpha),     alpha  = 0.90; end
    if nargin < 3 || isempty(Nmc),       Nmc    = 500;  end
    if nargin < 4 || isempty(seed),      seed   = 42;   end
    if nargin < 5 || isempty(whichQuant),whichQuant = alpha; end

    clc;
    fprintf('>>> plot_tail_voltage: start (lambda=%g, alpha=%.2f, Nmc=%d, seed=%d, q=%.2f)\n', ...
            lambda, alpha, Nmc, seed, whichQuant);

    %  folders
    here   = fileparts(mfilename('fullpath'));
    figDir = fullfile(here, 'figures');
    if ~exist(figDir,'dir'), mkdir(figDir); end

    %  recreate scenarios (so it aligns with other plots)
    rng(seed,'twister');
    C = config14();
    S = build_scenarios(C, Nmc, 'p_outage', 0.05, 'timeMode', "random");

    % compute per-scenario J for the baseline and pick tail index
    [J_base, okb] = compute_J_plain(C, S, lambda);
    idxValid = find(okb);
    if isempty(idxValid)
        error('No valid baseline scenarios to select a tail case.');
    end
    Jv = J_base(idxValid);
    [Jv_sorted, order] = sort(Jv, 'ascend');
    k = max(1, min(numel(Jv_sorted), ceil(whichQuant * numel(Jv_sorted))));
    idxTail = idxValid(order(k));
    fprintf('Selected scenario s=%d as tail case: J_base=%.3f (quantile q=%.2f)\n', ...
            idxTail, J_base(idxTail), whichQuant);

    % solve that scenario: baseline vs BUS-9
    % Build outage-aware network for *this* scenario
    [net, ren] = scenario_network_and_inj(C, S(idxTail));

    % baseline solve
    resB = solver_plain(net.Y, net.pairs, net.Zser, C.Sbase, ren.Pg_pu, ren.Qg_pu, ren.Pd_pu, ren.Qd_pu);

    % BUS-9 solve (same scenario)
    res9 = solver_bus(net.Y, net.pairs, net.Zser, C.Sbase, 9, C.Vref, C.SVC_Q_MAX, ...
                      ren.Pg_pu, ren.Qg_pu, ren.Pd_pu, ren.Qd_pu);

    assert(resB.converged && res9.converged, 'One of the solves did not converge.');

    % --- check losses & volt penalty for info
    [phiB, badB] = volt_penalty(resB.V_abs(2:end), C.Vband);
    [phi9, bad9] = volt_penalty(res9.V_abs(2:end), C.Vband);
    JB = resB.totPloss + lambda*phiB;
    J9 = res9.totPloss + lambda*phi9;

    % plot
    f = figure('Color','w','Name','Tail scenario: voltage profile (Baseline vs BUS-9)');
    hold on; grid on;
    nb = numel(resB.V_abs);
    bidx = 1:nb;
    plot(bidx, resB.V_abs, '-o', 'LineWidth', 1.8, 'MarkerSize', 5);
    plot(bidx, res9.V_abs, '-s', 'LineWidth', 1.8, 'MarkerSize', 5);

    yline(C.Vband(1), ':', 'LineWidth', 1.0);
    yline(C.Vband(2), ':', 'LineWidth', 1.0);

    xlabel('Bus index');
    ylabel('|V| (pu)');
    title(sprintf('Tail scenario (seed=%d, N=%d): J_{base}=%.2f, J_{BUS-9}=%.2f', seed, Nmc, JB, J9));
    legend({'Baseline (no SVC)', 'BUS-9', 'V_{min}', 'V_{max}'}, 'Location','best');

    % annotate SVC info
    txt = sprintf('BUS-9 SVC: Q=%.2f MVAr  (mode=%s)', res9.Qsvc_MVAr, res9.svcType);
    xl = xlim; yl = ylim;
    text(xl(1)+0.02*(xl(2)-xl(1)), yl(2)-0.05*(yl(2)-yl(1)), txt, 'FontWeight','bold');

    % save
    outPNG = fullfile(figDir, ...
        sprintf('voltage_tail_seed%d_N%d_q%02d.png', seed, Nmc, round(100*whichQuant)));
    exportgraphics(f, outPNG, 'Resolution', 200);
    fprintf('>>> wrote %s\n', outPNG);
    fprintf('>>> plot_tail_voltage: done\n');
end



function [J, ok] = compute_J_plain(C, S, lambda)
    Y0    = buildY(C);
    Sbase = C.Sbase;
    Nsc   = numel(S);
    J     = nan(Nsc,1);
    ok    = false(Nsc,1);

    for s = 1:Nsc
        % outage-aware
        if S(s).outage == "none"
            Y = Y0; pairs_s = C.pairs; Zser_s = C.Zser;
        else
            ab = sscanf(S(s).outage,'%d-%d');
            [Y, pairs_s, Zser_s] = remove_line_general(Y0, C, ab(1), ab(2));
        end
        if is_islanded(pairs_s, max(C.pairs(:)))
            continue;
        end

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

        res = solver_plain(Y, pairs_s, Zser_s, Sbase, Pg, Qg, Pd, Qd);
        if ~(isfield(res,'converged') && res.converged), continue; end
        Vabs = res.V_abs(:);
        if any(~isfinite(Vabs)) || any(Vabs < 0.8 | Vabs > 1.2), continue; end

        [phiV,~] = volt_penalty(Vabs(2:end), C.Vband);
        J(s)     = res.totPloss + lambda * phiV;
        ok(s)    = isfinite(J(s));
    end
end

function [net, ren] = scenario_network_and_inj(C, S1)
    % Build Y and active lists for *one* scenario + PU injections with renewables.
    Y0 = buildY(C);
    if S1.outage == "none"
        net.Y = Y0; net.pairs = C.pairs; net.Zser = C.Zser;
    else
        ab = sscanf(S1.outage,'%d-%d');
        [net.Y, net.pairs, net.Zser] = remove_line_general(Y0, C, ab(1), ab(2));
    end
    if is_islanded(net.pairs, max(C.pairs(:)))
        error('Selected tail scenario is islanded.');
    end

    Sbase = C.Sbase;
    ren.Pd_pu = (C.Pd_MW   .* (1 + S1.epsLoad)) / Sbase;
    ren.Qd_pu = (C.Qd_MVAr .* (1 + S1.epsLoad)) / Sbase;
    ren.Pg_pu =  C.Pg_MW / Sbase;
    ren.Qg_pu =  C.Qg_MVAr / Sbase;

    Ppv_MW   = pvPlantMW_datasheet(S1.G,  S1.TA,  C.PV);
    Pwind_MW = windPlantMW_formula(S1.Vw, C.WF);
    if isfield(C,'bus_solar') && ~isempty(C.bus_solar)
        ren.Pg_pu(C.bus_solar) = ren.Pg_pu(C.bus_solar) + Ppv_MW / Sbase;
    end
    if isfield(C,'bus_wind') && ~isempty(C.bus_wind)
        ren.Pg_pu(C.bus_wind)  = ren.Pg_pu(C.bus_wind)  + Pwind_MW / Sbase;
    end
end
