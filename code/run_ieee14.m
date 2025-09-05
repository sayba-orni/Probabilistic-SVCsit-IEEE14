function run_ieee14
clc; close all;
fprintf('>>> run_ieee14: start\n');
dbstop if error

% reproducible randomness for the main run
rng(42,'twister');

fprintf('>>> loading config14() ... ');
C  = config14();                 % no args
assert(~isempty(C) && isstruct(C), 'config14 returned empty/not a struct');
fprintf('ok\n');
need = {'Sbase','Vband','Vref','pairs','Zser','Bend','svc_buses','mid_list','data','PV','WF','bus_wind','bus_solar'};
for k = 1:numel(need)
    assert(isfield(C,need{k}), 'config14 missing field: %s', need{k});
end

Nmc = 500;
fprintf('>>> building %d scenarios ... ', Nmc);
S = build_scenarios(C, Nmc, 'p_outage', 0.05, 'timeMode', "random");
assert(numel(S)==Nmc, 'build_scenarios did not return Nmc scenarios');
fprintf('ok\n');

fprintf('>>> building candidate set ... ');
k = 0; candidates = struct('type',{},'k',{},'a',{},'b',{});
for b = C.svc_buses
    k=k+1; candidates(k) = struct('type',"bus",'k',b,'a',[],'b',[]);
end
for e = 1:size(C.mid_list,1)
    k=k+1; candidates(k) = struct('type',"mid",'k',[],'a',C.mid_list(e,1),'b',C.mid_list(e,2));
end
assert(~isempty(candidates), 'no candidates constructed');
fprintf('ok (%d candidates)\n', numel(candidates));

lambda = 100; alpha = 0.90;
fprintf('>>> ranking (lambda=%g, alpha=%g) ... ', lambda, alpha);
T = rank_candidates(C, S, candidates, lambda, alpha);
assert(istable(T) && height(T)>0, 'rank_candidates returned empty');
fprintf('ok (%d rows)\n', height(T));

if exist('baseline_no_svc','file')==2
    B = baseline_no_svc(C, S, lambda, alpha);
    fprintf(['Baseline (no SVC): meanLoss=%.3f MW, meanPhi=%.5f pu-sum, ', ...
             'violProb=%.3f, meanScore=%.3f, CVaR=%.3f\n'], ...
            B.meanLoss, B.meanPhi, B.violProb, B.meanScore, B.cvarScore);
else
    warning('baseline_no_svc.m not found — skipping baseline print.');
    B = struct();
end

set(0,'DefaultFigureVisible','on');

thisFile = mfilename('fullpath');
scriptDir = fileparts(thisFile);    % folder containing run_ieee14.m
projDir   = scriptDir;              
figDir    = fullfile(projDir, 'figures');
resDir    = fullfile(projDir, 'results');

if ~exist(resDir, 'dir'), mkdir(resDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end

outCSV = fullfile(resDir, 'ranking14.csv');
writetable(T, outCSV);
fprintf('>>> wrote %s\n', outCSV);


if ~isempty(fieldnames(B))
    outBASE = fullfile(resDir, 'baseline_no_svc.txt');
    fid = fopen(outBASE,'w');
    fprintf(fid, 'meanLoss=%.6f\nmeanPhi=%.6f\nviolProb=%.6f\nmeanScore=%.6f\ncvarScore=%.6f\n', ...
        B.meanLoss, B.meanPhi, B.violProb, B.meanScore, B.cvarScore);
    fclose(fid);
    fprintf('>>> wrote %s\n', outBASE);
end

% bar chart of CVaR
f = figure('Color','w','Visible','on','Name','IEEE-14 SVC ranking');
bar(T.cvarScore); grid on
xticks(1:height(T));
xticklabels(T.name);
xtickangle(45);
ylabel(sprintf('CVaR_{%.2f}(Loss + \\lambda·VoltPen)', alpha));
title('IEEE 14-bus — SVC candidate ranking (lower is better)');
drawnow;

outPNG = fullfile(figDir, 'ranking14.png');
exportgraphics(f, outPNG, 'Resolution', 200);
fprintf('>>> wrote %s\n', outPNG);


DO_REPEAT = true;               % <— set false to skip
if DO_REPEAT
    if exist('repeatability_check','file')==2
        seeds = [21, 42, 84];
        repeatability_check(C, seeds, Nmc, lambda, alpha);
    else
        warning('repeatability_check.m not found — skipping repeatability.');
    end
end

fprintf('>>> run_ieee14: done\n');
end
