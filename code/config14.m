function C = config14()
% CONFIG14 — IEEE14 bus with all synchronous generators turned off except small reactive support on bus 2 . user can modify this as necessary for their case

% system / SVC settings 
C.Sbase     = 100;                % MVA
C.Vband     = [0.95 1.05];        % acceptable |V| band
C.Vref      = 1.00;               % SVC voltage setpoint
C.SVC_Q_MAX = 300;                % MVAr rating (symmetric)
C.Bmax      =  (C.SVC_Q_MAX) / (C.Vref^2 * C.Sbase);
C.Bmin      = -C.Bmax;

% network topology (IEEE-14, per-unit on 100 MVA) 
% Columns: from  to   R(pu)     X(pu)     B_total(pu)
raw = [ ...
   1   2   0.01938   0.05917   0.0528
   1   5   0.05403   0.22304   0.0492
   2   3   0.04699   0.19797   0.0438
   2   4   0.05811   0.17632   0.0340
   2   5   0.05695   0.17388   0.0346
   3   4   0.06701   0.17103   0.0128
   4   5   0.01335   0.04211   0.0000
   4   7   0.00000   0.20912   0.0000
   4   9   0.00000   0.55618   0.0000
   5   6   0.00000   0.25202   0.0000
   6  11   0.09498   0.19890   0.0000
   6  12   0.12291   0.25581   0.0000
   6  13   0.06615   0.13027   0.0000
   7   8   0.00000   0.17615   0.0000
   7   9   0.00000   0.11001   0.0000
   9  10   0.03181   0.08450   0.0000
   9  14   0.12711   0.27038   0.0000
  10  11   0.08205   0.19207   0.0000
  12  13   0.22092   0.19988   0.0000
  13  14   0.17093   0.34802   0.0000 ];

C.pairs = raw(:,1:2);
R       = raw(:,3);
X       = raw(:,4);
Btot    = raw(:,5);

% series impedance and per-end shunt susceptance
C.Zser  = R + 1i*X;          % length = nLines
C.Bend  = 1i*(Btot/2);       % per end: j*B_total/2

% buses for SVC search 
C.svc_buses = 2:14;          % all non-slack buses
C.mid_list  = C.pairs;       % all lines as mid-line candidates

% base injections (MW / MVAr) 
Pd = zeros(1,14); Qd = zeros(1,14);
Pd( 2)=21.7;  Qd( 2)=12.7;
Pd( 3)=94.2;  Qd( 3)=19.0;
Pd( 4)=47.8;  Qd( 4)= 7.6;
Pd( 5)= 7.6;  Qd( 5)= 1.6;
Pd( 9)=29.5;  Qd( 9)=16.6;
Pd(10)= 9.0;  Qd(10)= 5.8;
Pd(11)= 3.5;  Qd(11)= 1.8;
Pd(12)= 6.1;  Qd(12)= 1.6;
Pd(13)=13.5;  Qd(13)= 5.8;
Pd(14)=14.9;  Qd(14)= 5.0;

C.Pd_MW   = Pd;
C.Qd_MVAr = Qd;

% Scheduled conventional generators (slack balances)
C.Pg_MW   = zeros(1,14);
C.Qg_MVAr = zeros(1,14);
C.Qg_MVAr(2) = 30;     % small local MVAr support (optional)

%  where PV / wind are injected 
C.bus_wind  = 9;       % wind at bus 9
C.bus_solar = 14;      % PV at bus 14

% PV module / plant (located Netrokona) 
C.PV.N_modules = round(30e6 / 280);  % ~30 MW plant, 280 W/module
C.PV.VMPP  = 31.28;
C.PV.IMPP  = 8.95;
C.PV.VOC   = 37.82;
C.PV.ISC   = 9.47;
C.PV.NOT   = 43;                      % °C
C.PV.Kv    = -0.0032 * C.PV.VOC;      % V/°C  (-0.32% Voc/°C)
C.PV.Ki    = -0.0006 * C.PV.ISC;      % A/°C  (-0.06% Isc/°C)

% Wind farm in Cox's Bazar rated 15 MW 
C.WF.PrMW  = 15;
C.WF.v_in  = 3;        % m/s (cut-in)
C.WF.v_r   = 12.5;     % m/s (rated)
C.WF.v_out = 25;       % m/s (cut-out)

% data files 
here       = fileparts(mfilename('fullpath'));   % .../code
projRoot   = fileparts(here);                    % project root
C.data.pv  = fullfile(projRoot,'data','pvdata2.csv');
C.data.wind= fullfile(projRoot,'data','windspeed.csv');

end
