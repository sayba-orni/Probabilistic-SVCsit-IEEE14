function res = solver_bus(Y_base, pairs, Zser, Sbase, svc_bus, Vref_svc, SVC_Q_MAX, Pg, Qg, Pd, Qd)
% SOLVER_BUS  Robust NR power flow with a single SVC at an existing bus.
% Active-set NR + Armijo backtracking, variable-susceptance SVC model.

    N = size(Y_base,1);

    % MVAr limit -> susceptance bounds (pu)
    Bmax =  (SVC_Q_MAX) / (max(Vref_svc,1e-6)^2 * Sbase);
    Bmin = -Bmax;

    % net schedules (pu)
    Psch = Pg - Pd;
    Qsch = Qg - Qd;

    % flat start
    Vmag = ones(1,N); Vmag(1) = 1.06;
    Vang = zeros(1,N);
    Bsvc = 0.0; 
    mode = "REG";      % REG = regulate |Vsvc| ; SAT = B at limits

    tol = 1e-6; max_iter = 80; iter = 0; err = 1;

    while err>tol && iter<max_iter
        iter = iter + 1;

        % Y with current SVC susceptance
        Yit = Y_base;  
        Yit(svc_bus,svc_bus) = Yit(svc_bus,svc_bus) + 1i*Bsvc;

        if mode=="REG", Vmag(svc_bus) = Vref_svc; end

        % mismatch (buses 2..N)
        [P,Q] = pq_injections(Vmag,Vang,Yit);
        dP = Psch - P;  dQ = Qsch - Q;
        M  = [dP(2:N) dQ(2:N)]';
        err = max(abs(M)); if err<=tol, break; end

        % Jacobian
        [J1,J2,J3,J4] = jac_blocks(Vmag,Vang,Yit);
        if mode=="REG"
            idx = svc_bus-1;     % column in the 2..N block
            J2(:,idx) = 0;
            J4(:,idx) = 0;  
            J4(idx,idx) = -Vmag(svc_bus)^2;   % dQ/dB at SVC bus
        end
        J  = [J1 J2; J3 J4];

        % Newton step (tiny Tikhonov regularization)
        dx = (J + 1e-12*eye(size(J))) \ M; 
        dx = dx.';
        PQ = N-1; 
        dth = dx(1:PQ); 
        dVm = dx(PQ+1:end);
        dB  = (mode=="REG") * dVm(svc_bus-1);

        % Armijo backtracking + bounds projection
        alpha = 1.0; accepted=false; merit0 = 0.5*(M.'*M);
        while alpha > 1/1024
            Vang_try = Vang;  
            Vang_try(2:N) = Vang(2:N) + alpha*dth;

            dV_full = zeros(1,N); dV_full(2:N) = alpha*dVm;
            Vmag_try = Vmag + dV_full;
            if mode=="REG", Vmag_try(svc_bus) = Vref_svc; end

            B_try   = Bsvc + alpha*dB;
            mode_try = mode;
            if mode=="REG" && (B_try<Bmin || B_try>Bmax)
                B_try    = min(max(B_try,Bmin),Bmax);
                mode_try = "SAT";
            end

            Ytr = Y_base; 
            Ytr(svc_bus,svc_bus) = Ytr(svc_bus,svc_bus) + 1i*B_try;
            [Pt,Qt] = pq_injections(Vmag_try,Vang_try,Ytr);
            Mt = [ (Psch(2:N)-Pt(2:N))' ; (Qsch(2:N)-Qt(2:N))' ];
            merit_try = 0.5*(Mt.'*Mt);

            if merit_try < merit0
                Vang = Vang_try; Vmag = Vmag_try; Bsvc = B_try; mode = mode_try;
                accepted = true; 
                break;
            else
                alpha = alpha/2;
            end
        end

        if ~accepted
            % minimal damped step
            Vang(2:N) = Vang(2:N) + 1e-3*dth;
            if mode=="REG"
                Vmag(2:N) = Vmag(2:N) + 1e-3*dVm; 
                Vmag(svc_bus) = Vref_svc;
                Bsvc = min(max(Bsvc + 1e-3*dB, Bmin), Bmax);
            else
                Vmag(2:N) = Vmag(2:N) + 1e-3*dVm;
            end
        end

        % re-enter regulation if interior & near setpoint
        if mode=="SAT"
            if Bsvc>Bmin+1e-6 && Bsvc<Bmax-1e-6 && abs(Vmag(svc_bus)-Vref_svc)<5e-4
                mode="REG"; Vmag(svc_bus)=Vref_svc;
            end
        end
    end

    % final quantities
    Yf = Y_base; 
    Yf(svc_bus,svc_bus) = Yf(svc_bus,svc_bus) + 1i*Bsvc;
    Vdeg = Vang*180/pi;

    Pg1_MW   = Sbase * calc_injection_P(1, Vmag, Vang, Yf);
    Qg1_MVAr = Sbase * calc_injection_Q(1, Vmag, Vang, Yf);

    Qsvc_MVAr = -(Vmag(svc_bus)^2) * Bsvc * Sbase;
    svcType   = 'Inductive'; 
    if Qsvc_MVAr >= 0, svcType = 'Capacitive'; end

    % currents & losses for the actual network
    lines = pairs; 
    Zlist = Zser;
    L = size(lines,1); 
    I_mag = zeros(L,1); I_ang = zeros(L,1); 
    Pl = zeros(L,1);    Ql = zeros(L,1);
    for ee = 1:L
        a = lines(ee,1); b = lines(ee,2);
        [I_mag(ee), I_ang(ee), Pl(ee), Ql(ee)] = current_and_lineloss( ...
            Vmag(a), Vdeg(a), Vmag(b), Vdeg(b), Zlist(ee), Sbase);
    end

    res.iter      = iter;
    res.err       = err;
    res.converged = (err<=tol);
    res.mode      = char(mode);
    res.svc_bus   = svc_bus;
    res.B_svc     = Bsvc;
    res.V_abs     = Vmag;
    res.V_deg     = Vdeg;
    res.Qsvc_MVAr = Qsvc_MVAr;
    res.svcType   = svcType;

    res.genMW     = [Pg1_MW,   Pg(2:end)*Sbase];
    res.genMVAr   = [Qg1_MVAr, Qg(2:end)*Sbase];
    res.loadMW    = Pd * Sbase;
    res.loadMVAr  = Qd * Sbase;

    res.lines     = lines;  
    res.I_mag     = I_mag; 
    res.I_ang     = I_ang;
    res.Lmw       = Pl;     
    res.Lmvar     = Ql;
    res.totPloss  = sum(Pl); 
    res.totQloss  = sum(Ql);
end

% LOCAL HELPERS 

function [P,Q] = pq_injections(Vmag, Vang, Y)
    N = numel(Vmag); 
    P = zeros(1,N); 
    Q = zeros(1,N);
    for i = 1:N
        for k = 1:N
            th = angle(Y(i,k)) + Vang(k) - Vang(i);
            g  = abs(Y(i,k));
            P(i) = P(i) + Vmag(i)*Vmag(k)*g*cos(th);
            Q(i) = Q(i) - Vmag(i)*Vmag(k)*g*sin(th);
        end
    end
end

function [J1,J2,J3,J4] = jac_blocks(Vmag, Vang, Y)
    N = numel(Vmag); 
    n = N-1;
    J1 = zeros(n,n); J2 = zeros(n,n); 
    J3 = zeros(n,n); J4 = zeros(n,n);
    [P,Q] = pq_injections(Vmag, Vang, Y);
    for i = 2:N
        for k = 2:N
            if i == k
                Yii = Y(i,i); gii = abs(Yii); thi = angle(Yii);
                J1(i-1,k-1) = -Q(i) - Vmag(i)^2*gii*sin(thi);
                J2(i-1,k-1) =  P(i) + Vmag(i)^2*gii*cos(thi);
                J3(i-1,k-1) =  P(i) - Vmag(i)^2*gii*cos(thi);
                J4(i-1,k-1) =  Q(i) - Vmag(i)^2*gii*sin(thi);
            else
                th = angle(Y(i,k)) + Vang(k) - Vang(i);
                g  = abs(Y(i,k));
                J1(i-1,k-1) = -Vmag(i)*Vmag(k)*g*sin(th);
                J2(i-1,k-1) =  Vmag(i)*Vmag(k)*g*cos(th);
                J3(i-1,k-1) = -Vmag(i)*Vmag(k)*g*cos(th);
                J4(i-1,k-1) = -Vmag(i)*Vmag(k)*g*sin(th);
            end
        end
    end
end

function [Iabs,Iang,Ploss,Qloss] = current_and_lineloss(Vm_i, ang_i_deg, Vm_j, ang_j_deg, Z, Sbase)
    Vi = Vm_i * exp(1i*ang_i_deg*pi/180);
    Vj = Vm_j * exp(1i*ang_j_deg*pi/180);
    I  = (Vi - Vj) / Z;
    Iabs = abs(I); 
    Iang = angle(I)*180/pi;
    S_loss_pu = (Iabs^2) * Z;            % I^2 * Z
    Ploss = Sbase * real(S_loss_pu);
    Qloss = Sbase * imag(S_loss_pu);
end

function P = calc_injection_P(i, Vmag, Vang, Y)
    N = numel(Vmag); 
    P = 0;
    for k = 1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        P  = P + Vmag(i)*Vmag(k)*abs(Y(i,k))*cos(th);
    end
end

function Q = calc_injection_Q(i, Vmag, Vang, Y)
    N = numel(Vmag); 
    Q = 0;
    for k = 1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        Q  = Q - Vmag(i)*Vmag(k)*abs(Y(i,k))*sin(th);
    end
end
