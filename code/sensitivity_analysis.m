function sensitivity_analysis
% Sensitivity sweeps for lambda and SVC rating.
% - Uses the SAME scenario set for all runs (seeded) to isolate settings.
% - Produces CSV + figures showing whether BUS-9 (and 9-corridor mids)
%   remain on top.

clc; fprintf('>>> sensitivity_analysis: start\n');
dbstop if error

% fixed experiment knobs 
seed     = 42;          % repeatability
Nmc      = 500;         % scenarios
alpha    = 0.90;        % CVaR tail
lambdas  = [50 100 150];
ratings  = [150 300 450];   % MVAr (symmetric +/-)

%  build once: config + scenarios (seeded), and candidate set ---
rng(seed,'twister');
C0 = config14();
assert(all(isfield(C0,{'Sbase','Vband','Vref','pairs','Zser','Bend','svc_buses','mid_list','PV','WF','bus_wind','bus_solar'})), ...
    'config14 check failed');

fprintf('>>> building scenarios (N=%d, seed=%d) ... ', Nmc, seed);
S  = build_scenarios(C0, Nmc, 'p_outage',0.05, 'timeMode',"random");
fprintf('ok\n');

% candidate list
cand = mkcands(C0);

% dirs
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
figDir = fullfile(here, 'figures');
if ~exist(resDir,'dir'), mkdir(resDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end

% containers
Rows = [];

% A) lambda sweep (keep rating = config's default) 
fprintf('\n--- Lambda sweep ---\n');
for L = lambdas
    C = C0;  % unchanged rating
    T = rank_candidates(C, S, cand, L, alpha);
    T = sortrows(T,'cvarScore','ascend');

    % record top-5 names
    topNames = string(T.name(1:min(5,height(T))));
    Rows = [Rows; pack_row("lambda", L, NaN, topNames)]; %#ok<AGROW>

    % print short report
    fprintf('lambda=%-4g  top-5: %s\n', L, strjoin(cellstr(topNames)', ', '));

    % quick figure: CVaR of top-10 for each lambda
    f = figure('Color','w','Visible','off','Name',sprintf('CVaR by candidate (lambda=%g)',L));
    k = min(10, height(T));
    bar(T.cvarScore(1:k)); grid on
    xticks(1:k); xticklabels(T.name(1:k)); xtickangle(40);
    ylabel(sprintf('CVaR_{%.2f}(Loss + \\lambda·VoltPen)', alpha));
    title(sprintf('Top-%d CVaR — \\lambda=%g', k, L));
    drawnow;
    outPNG = fullfile(figDir, sprintf('sensitivity_lambda_%g.png', L));
    exportgraphics(f, outPNG, 'Resolution', 200);
    close(f);
    fprintf('>>> wrote %s\n', outPNG);
end

%  B) rating sweep (keep lambda = 100)
fprintf('\n--- SVC rating sweep ---\n');
lambda = 100;
for Rmvar = ratings
    C = C0;
    C.SVC_Q_MAX = Rmvar;                       % MVAr
    C.Bmax      =  (C.SVC_Q_MAX) / (C.Vref^2 * C.Sbase);
    C.Bmin      = -C.Bmax;

    T = rank_candidates(C, S, cand, lambda, alpha);
    T = sortrows(T,'cvarScore','ascend');

    % record top-5 names
    topNames = string(T.name(1:min(5,height(T))));
    Rows = [Rows; pack_row("rating", NaN, Rmvar, topNames)]; %#ok<AGROW>

    % print short report
    fprintf('rating=%-3g MVAr  top-5: %s\n', Rmvar, strjoin(cellstr(topNames)', ', '));

    % quick figure: CVaR of top-10 for each rating
    f = figure('Color','w','Visible','off','Name',sprintf('CVaR by candidate (rating=%g MVAr)',Rmvar));
    k = min(10, height(T));
    bar(T.cvarScore(1:k)); grid on
    xticks(1:k); xticklabels(T.name(1:k)); xtickangle(40);
    ylabel(sprintf('CVaR_{%.2f}(Loss + \\lambda·VoltPen)', alpha));
    title(sprintf('Top-%d CVaR — SVC rating=%g MVAr', k, Rmvar));
    drawnow;
    outPNG = fullfile(figDir, sprintf('sensitivity_rating_%gMVAr.png', Rmvar));
    exportgraphics(f, outPNG, 'Resolution', 200);
    close(f);
    fprintf('>>> wrote %s\n', outPNG);
end

%  write summary CSV 
Tsum = struct2table(Rows);
outCSV = fullfile(resDir, sprintf('sensitivity_summary_seed%d_N%d.csv', seed, Nmc));
writetable(Tsum, outCSV);
fprintf('\n>>> wrote %s\n', outCSV);

%  quick sanity statements (BUS-9 + corridor check) 
isBus9 = any(strcmp(Tsum.top1, "BUS-9"));
corridorTags = ["MID-9-10","MID-7-9","MID-9-14"];
hasCorridorTop3 = arrayfun(@(i) any(ismember(string(Tsum{i,["top1","top2","top3"]}),corridorTags)), 1:height(Tsum));

fprintf('\nSummary checks:\n');
fprintf(' - BUS-9 appears as top-1 in %d/%d cases.\n', isBus9 + 0*sum(~isBus9), height(Tsum)); % simple boolean print
fprintf(' - 9-corridor (MID-9-10 / MID-7-9 / MID-9-14) appears within top-3 in %d/%d cases.\n', sum(hasCorridorTop3), height(Tsum));

fprintf('>>> sensitivity_analysis: done\n');

function cand = mkcands(C)
    k=0; cand = struct('type',{},'k',{},'a',{},'b',{});
    for b = C.svc_buses, k=k+1; cand(k)=struct('type',"bus",'k',b,'a',[],'b',[]); end
    for e=1:size(C.mid_list,1), k=k+1; cand(k)=struct('type',"mid",'k',[],'a',C.mid_list(e,1),'b',C.mid_list(e,2)); end
end

function row = pack_row(kind, lam, rat, topNames)
    % pad to 5 entries
    tpad = strings(1,5); tpad(1:numel(topNames)) = topNames;
    row = struct('kind',string(kind), ...
                 'lambda',lam, ...
                 'ratingMVAr',rat, ...
                 'top1',tpad(1),'top2',tpad(2),'top3',tpad(3),'top4',tpad(4),'top5',tpad(5));
end
end
