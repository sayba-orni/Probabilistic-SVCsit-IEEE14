function R = evaluate_candidate(cand, C, S)
% Evaluate one SVC candidate (bus or mid-line) across all scenarios S.
% Relies on remove_line_general() returning:
%   [Y_new, pairs_s, Zser_s, Bend_s] for the post-outage network.

    Y_base = buildY(C);
    Sbase  = C.Sbase;
    Vref   = C.Vref;

    Nsc = numel(S);
    lossMW = nan(Nsc,1);
    phiV   = nan(Nsc,1);
    vviol  = false(Nsc,1);
    qsvc   = nan(Nsc,1);
    conv   = false(Nsc,1);

    for s = 1:Nsc
        % Outage-aware Y and active line set 
        if S(s).outage == "none"
            Y       = Y_base;
            pairs_s = C.pairs;
            Zser_s  = C.Zser;
            Bend_s  = C.Bend;
            out_ab  = [];
        else
            ab     = sscanf(S(s).outage, '%d-%d');
            out_ab = ab(:).';
            % IMPORTANT: use the 4-output remover so pairs/Z are consistent
            [Y, pairs_s, Zser_s, Bend_s, tap_s] = remove_line_general(Y_base, C, out_ab(1), out_ab(2));
        end

        % If the outage islands the grid, skip this scenario (N/A)
        if is_islanded(pairs_s, max(C.pairs(:)))
            continue;
        end

        % Loads & gens (pu) 
        Pd = (C.Pd_MW   .* (1 + S(s).epsLoad)) / Sbase;
        Qd = (C.Qd_MVAr .* (1 + S(s).epsLoad)) / Sbase;
        Pg =  C.Pg_MW   / Sbase;
        Qg =  C.Qg_MVAr / Sbase;

        % Renewables (pu) 
        Ppv_MW   = pvPlantMW_datasheet(S(s).G,  S(s).TA,  C.PV);
        Pwind_MW = windPlantMW_formula(S(s).Vw, C.WF);
        if isfield(C,'bus_solar') && ~isempty(C.bus_solar)
            Pg(C.bus_solar) = Pg(C.bus_solar) + Ppv_MW / Sbase;
        end
        if isfield(C,'bus_wind') && ~isempty(C.bus_wind)
            Pg(C.bus_wind)  = Pg(C.bus_wind)  + Pwind_MW / Sbase;
        end

        %  Solve PF with SVC model 
        if cand.type == "bus"
            % Pass the scenario's active line set so losses are computed on it
            res = solver_bus(Y, pairs_s, Zser_s, Sbase, cand.k, Vref, C.SVC_Q_MAX, Pg, Qg, Pd, Qd);

        else % cand.type == "mid"
            % If the candidate line is outaged or not present after filtering, skip
            if ~isempty(out_ab) && isequal(sort(out_ab), sort([cand.a cand.b]))
                continue;
            end
            eCase = find( (pairs_s(:,1)==cand.a & pairs_s(:,2)==cand.b) | ...
                          (pairs_s(:,1)==cand.b & pairs_s(:,2)==cand.a), 1);
            if isempty(eCase)
                continue;
            end
            res = solver_mid(pairs_s, Zser_s, Bend_s, Sbase, eCase, Vref, C.Bmin, C.Bmax, Pg, Qg, Pd, Qd);
        end

        % Validate solution 
        conv(s) = isfield(res,'converged') && res.converged;
        if ~conv(s), continue; end

        Vabs = res.V_abs(:);
        if any(~isfinite(Vabs)) || any(Vabs < 0.8 | Vabs > 1.2)
            conv(s) = false; continue;
        end

        % Metrics 
        lossMW(s) = res.totPloss;

        % Voltage penalty on non-slack buses (exclude bus 1)
        if numel(Vabs) >= 2
            Vcheck = Vabs(2:end);
        else
            Vcheck = Vabs;
        end
        [phiV(s), vviol(s)] = volt_penalty(Vcheck, C.Vband);

        if isfield(res,'Qsvc_MVAr')
            qsvc(s) = res.Qsvc_MVAr;
        end
    end

    R = struct('lossMW',lossMW, 'phiV',phiV, 'vviol',vviol, ...
               'qsvc',qsvc,     'converged',conv);
end
