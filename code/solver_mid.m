function res = solver_mid(pairs, Zser, Bend, Sbase, eCase, Vref, Bmin, Bmax, Pg, Qg, Pd, Qd)
% SOLVER_MID  Robust NR power flow with an SVC at the MIDPOINT of line eCase.
% - Splits the selected line into two half-π sections with a new bus m=N0+1.
% - Bend is the per-end shunt susceptance (j*B_end), i.e., j*(B_total/2).

    %  Build base Y for the current scenario network 
    N0 = max(pairs(:));
    Y0 = zeros(N0,N0);
    for e = 1:size(pairs,1)
        a = pairs(e,1); b = pairs(e,2);
        y = 1/Zser(e);   B = Bend(e);              % per-end shunt jB_end
        Y0(a,a)=Y0(a,a)+y+B;  Y0(b,b)=Y0(b,b)+y+B;
        Y0(a,b)=Y0(a,b)-y;    Y0(b,a)=Y0(b,a)-y;
    end

    %  Select & split the target line (a,b) 
    a = pairs(eCase,1); 
    b = pairs(eCase,2);
    Zab = Zser(eCase);             
    Bab = Bend(eCase);             % j*(B_total/2) per end
    y_ab = 1/Zab;

    m = N0 + 1;                    % new midpoint bus index
    Y  = zeros(m,m); 
    Y(1:N0,1:N0) = Y0;

    % remove original (a,b) branch (series + per-end shunts)
    Y(a,a) = Y(a,a) - y_ab - Bab;
    Y(b,b) = Y(b,b) - y_ab - Bab;
    Y(a,b) = 0; 
    Y(b,a) = 0;

    % add two half-π sections a–m and m–b
    z_half = Zab/2;                % series Z/2 on each half
    y_half = 1/z_half;
    B_half = Bab/2;                % per-end jB becomes jB/2 on each new half-end
    % a–m
    Y(a,a) = Y(a,a) + y_half + B_half;
    Y(m,m) = Y(m,m) + y_half + B_half;
    Y(a,m) = Y(a,m) - y_half;  
    Y(m,a) = Y(m,a) - y_half;
    % m–b
    Y(b,b) = Y(b,b) + y_half + B_half;
    Y(m,m) = Y(m,m) + y_half + B_half;
    Y(b,m) = Y(b,m) - y_half;  
    Y(m,b) = Y(m,b) - y_half;

    % ---------- NR setup ----------
    N    = m;
    Vmag = ones(1,N); Vmag(1)=1.06; Vmag(2:N0)=1.00; Vmag(m)=Vref;
    Vang = zeros(1,N);
    Bsvc = 0.0; 
    mode = "REG";

    % schedules (append zeros for the new bus m)
    Psch = [Pg 0] - [Pd 0];
    Qsch = [Qg 0] - [Qd 0];

    tol = 1e-6; max_iter = 60; iter = 0; err = 1;

    % NR with active-set + backtracking 
    while err>tol && iter<max_iter
        iter = iter + 1;

        Yit = Y; 
        Yit(m,m) = Yit(m,m) + 1i*Bsvc;         % SVC susceptance at m

        [P,Q] = pq_injections(Vmag,Vang,Yit);
        dP = Psch - P; 
        dQ = Qsch - Q;
        M  = [dP(2:N) dQ(2:N)]';
        err = max(abs(M)); 
        if err<=tol, break; end

        [J1,J2,J3,J4] = jac_blocks(Vmag,Vang,Yit);
        if mode=="REG"
            idx = m-1;                 % column matching bus m in 2..N
            J2(:,idx) = 0; 
            J4(:,idx) = 0; 
            J4(idx,idx) = -Vmag(m)^2;  % dQ/dB at m
        end
        J = [J1 J2; J3 J4];

        dx = (J + 1e-12*eye(size(J))) \ M; 
        dx = dx.';
        PQ = N-1; 
        dth = dx(1:PQ); 
        dVm = dx(PQ+1:end);
        dB  = (mode=="REG") * dVm(m-1);

        % Armijo backtracking + bounds projection
        alpha=1.0; accepted=false; M0=0.5*(M.'*M);
        while alpha>1/1024
            Vang_try=Vang;  Vang_try(2:N)=Vang(2:N)+alpha*dth;
            dV=zeros(1,N);  dV(2:N)=alpha*dVm;
            Vmag_try=Vmag + dV; 
            if mode=="REG", Vmag_try(m)=Vref; end

            B_try = Bsvc + alpha*dB; mode_try=mode;
            if mode=="REG" && (B_try<Bmin || B_try>Bmax)
                B_try = min(max(B_try,Bmin),Bmax); 
                mode_try="SAT";
            end

            Ytr=Y; Ytr(m,m)=Ytr(m,m)+1i*B_try;
            [Pt,Qt]=pq_injections(Vmag_try,Vang_try,Ytr);
            Mt=[ (Psch(2:N)-Pt(2:N))' ; (Qsch(2:N)-Qt(2:N))' ];
            if 0.5*(Mt.'*Mt) < M0
                Vang=Vang_try; Vmag=Vmag_try; Bsvc=B_try; mode=mode_try; 
                accepted=true; 
                break;
            else
                alpha = alpha/2;
            end
        end

        if ~accepted
            % minimal damped step
            Vang(2:N)=Vang(2:N)+1e-3*dth;
            if mode=="REG"
                Vmag(2:N)=Vmag(2:N)+1e-3*dVm; Vmag(m)=Vref;
                Bsvc=min(max(Bsvc + 1e-3*dB, Bmin), Bmax);
            else
                Vmag(2:N)=Vmag(2:N)+1e-3*dVm;
            end
        end

        % re-enter regulation if back inside bounds near setpoint
        if mode=="SAT"
            if Bsvc>Bmin+1e-6 && Bsvc<Bmax-1e-6 && abs(Vmag(m)-Vref)<5e-4
                mode="REG"; Vmag(m)=Vref;
            end
        end
    end

    %  Post-processing 
    Yfinal = Y; 
    Yfinal(m,m) = Yfinal(m,m) + 1i*Bsvc;
    Vdeg = Vang*180/pi;

    % SVC MVAr (positive = capacitive)
    Qsvc_MVAr = -(Vmag(m)^2)*Bsvc*Sbase;
    svcType   = 'Inductive'; 
    if Qsvc_MVAr >= 0, svcType = 'Capacitive'; end

    % line currents & losses (replace a-b by a-m and m-b)
    basePairs = pairs; baseZ = Zser;
    keep = ~( (basePairs(:,1)==a & basePairs(:,2)==b) | ...
              (basePairs(:,1)==b & basePairs(:,2)==a) );
    basePairs = basePairs(keep,:); 
    baseZ     = baseZ(keep);

    lines = [basePairs; a m; m b];
    Zs    = [baseZ;     z_half; z_half];

    L = size(lines,1);
    I_mag=zeros(L,1); I_ang=zeros(L,1); Pl=zeros(L,1); Ql=zeros(L,1);
    for ee=1:L
        ii = lines(ee,1); 
        jj = lines(ee,2);
        [I_mag(ee), I_ang(ee), Pl(ee), Ql(ee)] = current_and_lineloss( ...
            Vmag(ii), Vdeg(ii), Vmag(jj), Vdeg(jj), Zs(ee), Sbase);
    end

    % slack injections (bus 1)
    Pg1_MW   = Sbase * calc_injection_P(1, Vmag, Vang, Yfinal);
    Qg1_MVAr = Sbase * calc_injection_Q(1, Vmag, Vang, Yfinal);

    % assemble output
    res.iter=iter; res.err=err; res.converged=(err<=tol); res.mode=char(mode);
    res.svc_bus=m; res.B_svc=Bsvc; res.Qsvc_MVAr=Qsvc_MVAr; res.svcType=svcType;
    res.V_abs=Vmag; res.V_deg=Vdeg;

    res.genMW    = [Pg1_MW,   Pg(2:end)*Sbase, 0];
    res.genMVAr  = [Qg1_MVAr, Qg(2:end)*Sbase, 0];
    res.loadMW   = [Pd 0] * Sbase;
    res.loadMVAr = [Qd 0] * Sbase;

    res.lines=lines; res.I_mag=I_mag; res.I_ang=I_ang;
    res.Lmw=Pl; res.Lmvar=Ql; res.totPloss=sum(Pl); res.totQloss=sum(Ql);
end



function [P,Q] = pq_injections(Vmag, Vang, Y)
    N=numel(Vmag); P=zeros(1,N); Q=zeros(1,N);
    for i=1:N
        for k=1:N
            th = angle(Y(i,k)) + Vang(k) - Vang(i);
            g  = abs(Y(i,k));
            P(i) = P(i) + Vmag(i)*Vmag(k)*g*cos(th);
            Q(i) = Q(i) - Vmag(i)*Vmag(k)*g*sin(th);
        end
    end
end

function [J1,J2,J3,J4] = jac_blocks(Vmag,Vang,Y)
    N=numel(Vmag); n=N-1;
    J1=zeros(n,n); J2=zeros(n,n); J3=zeros(n,n); J4=zeros(n,n);
    [P,Q] = pq_injections(Vmag,Vang,Y);
    for i=2:N
        for k=2:N
            if i==k
                Yii = Y(i,i); gii=abs(Yii); thi=angle(Yii);
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

function [Iabs,Iang,Ploss,Qloss] = current_and_lineloss(Vm_i,ang_i_deg,Vm_j,ang_j_deg,Z,Sbase)
    Vi = Vm_i * exp(1i*ang_i_deg*pi/180);
    Vj = Vm_j * exp(1i*ang_j_deg*pi/180);
    I  = (Vi - Vj) / Z;
    Iabs = abs(I); Iang = angle(I)*180/pi;
    S_loss_pu = (Iabs^2) * Z;      % I^2 * Z
    Ploss = Sbase * real(S_loss_pu);
    Qloss = Sbase * imag(S_loss_pu);
end

function P = calc_injection_P(i,Vmag,Vang,Y)
    N=numel(Vmag); P=0;
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        P = P + Vmag(i)*Vmag(k)*abs(Y(i,k))*cos(th);
    end
end

function Q = calc_injection_Q(i,Vmag,Vang,Y)
    N=numel(Vmag); Q=0;
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        Q  = Q - Vmag(i)*Vmag(k)*abs(Y(i,k))*sin(th);
    end
end
