function make_topk_table(k, Nmc, seed, lambda, alpha)
% MAKE_TOPK_TABLE  Build a compact table: Baseline vs top-k candidates.
% Outputs:
%   results/topk_table.csv
%   figures/topk_table.png  (simple figure rendering of the table)
%
% Usage:
%   make_topk_table                % defaults: k=10, Nmc=500, seed=42, lambda=100, alpha=0.90
%   make_topk_table(5,500,42,100,0.95)

    if nargin<1 || isempty(k),       k = 10;   end
    if nargin<2 || isempty(Nmc),     Nmc = 500; end
    if nargin<3 || isempty(seed),    seed = 42; end
    if nargin<4 || isempty(lambda),  lambda = 100; end
    if nargin<5 || isempty(alpha),   alpha = 0.90; end

    clc; fprintf('>>> make_topk_table: start (k=%d, Nmc=%d, seed=%d, lambda=%g, alpha=%.2f)\n',...
                 k, Nmc, seed, lambda, alpha);

    here   = fileparts(mfilename('fullpath'));
    resDir = fullfile(here, 'results');
    figDir = fullfile(here, 'figures');
    if ~exist(resDir,'dir'), mkdir(resDir); end
    if ~exist(figDir,'dir'), mkdir(figDir); end

    csvPath  = fullfile(resDir, 'ranking14.csv');
    basePath = fullfile(resDir, 'baseline_no_svc.txt');
    assert(exist(csvPath,'file')==2, 'Missing %s. Run run_ieee14 first.', csvPath);

    %  read ranking & pick top-k by CVaR
    T = readtable(csvPath, 'TextType','string');
    req = ["name","meanLoss","meanPhi","violProb","meanScore","cvarScore"];
    for j = 1:numel(req)
        assert(any(strcmpi(T.Properties.VariableNames, req(j))), ...
            'ranking14.csv missing column: %s', req(j));
    end
    T = sortrows(T, 'cvarScore', 'ascend');
    top = T(1:min(k,height(T)),:);


    haveBaseline = false;
    base = struct('name',"Baseline (no SVC)",'meanLoss',NaN,'meanPhi',NaN,'violProb',NaN,'meanScore',NaN,'cvarScore',NaN);
    if exist(basePath,'file')==2
        try
            txt = fileread(basePath);
            baseParsed = try_parse_baseline(txt);
            % normalize cvar/cvarScore
            if isfield(baseParsed,'cvarScore') && isfinite(baseParsed.cvarScore)
                base.cvarScore = baseParsed.cvarScore;
            elseif isfield(baseParsed,'cvar') && isfinite(baseParsed.cvar)
                base.cvarScore = baseParsed.cvar;
            end
            base.meanLoss  = baseParsed.meanLoss;
            base.meanPhi   = baseParsed.meanPhi;
            base.violProb  = baseParsed.violProb;
            base.meanScore = baseParsed.meanScore;
            haveBaseline = all(isfinite([base.meanLoss base.meanPhi base.violProb base.meanScore base.cvarScore]));
        catch
            haveBaseline = false;
        end
    end
    if ~haveBaseline
        fprintf('>>> recomputing baseline ... ');
        rng(seed,'twister');
        C  = config14();
        S  = build_scenarios(C, Nmc, 'p_outage',0.05, 'timeMode',"random");
        B  = baseline_no_svc(C, S, lambda, alpha);
        base.meanLoss  = B.meanLoss;
        base.meanPhi   = B.meanPhi;
        base.violProb  = B.violProb;
        base.meanScore = B.meanScore;
        base.cvarScore = B.cvarScore;
        fprintf('ok\n');
    end

    % assemble table (baseline row + top-k)
    name = ["BASELINE"; top.name];
    meanLoss  = [base.meanLoss;  top.meanLoss];
    meanPhi   = [base.meanPhi;   top.meanPhi];
    violProb  = [base.violProb;  top.violProb];
    meanScore = [base.meanScore; top.meanScore];
    cvarScore = [base.cvarScore; top.cvarScore];

    TT = table(name, meanLoss, meanPhi, violProb, meanScore, cvarScore);

    % save CSV
    outCSV = fullfile(resDir, sprintf('topk_table_k%d.csv', height(top)));
    writetable(TT, outCSV);
    fprintf('>>> wrote %s\n', outCSV);

    % quick figure of the table (lightweight for the paper's appendix)
    f = figure('Color','w','Name','Baseline vs Top-k (by CVaR)','Position',[100 100 900 420]);
    uit = uitable(f, 'Data', TT{:,:}, 'ColumnName', TT.Properties.VariableNames, ...
                      'RowName', {}, 'Units','normalized', 'Position',[0 0 1 1]);
    % widen the first col
    cw = uit.ColumnWidth;
    if iscell(cw)
        cw{1} = 200;
        uit.ColumnWidth = cw;
    else
        uit.ColumnWidth = {200, 'auto', 'auto', 'auto', 'auto', 'auto'};
    end
    drawnow;
    outPNG = fullfile(figDir, sprintf('topk_table_k%d.png', height(top)));
    exportgraphics(f, outPNG, 'Resolution', 200);
    fprintf('>>> wrote %s\n', outPNG);

    fprintf('>>> make_topk_table: done\n');
end



function base = try_parse_baseline(txt)
% Pull numbers by labels; supports units; accepts CVaR/cvarScore label variants.
base = struct('meanLoss',NaN,'meanPhi',NaN,'violProb',NaN,'meanScore',NaN,'cvarScore',NaN,'cvar',NaN);
lines = regexp(txt, '\r\n|\r|\n', 'split');
bln = '';
for i=1:numel(lines)
    if ~isempty(regexp(lines{i}, '(?i)baseline', 'once'))
        bln = strtrim(lines{i}); break;
    end
end
if isempty(bln), bln = strrep(txt, sprintf('\n'), ' '); end
base.meanLoss   = find_key_number(bln, {'meanLoss'});
base.meanPhi    = find_key_number(bln, {'meanPhi'});
base.violProb   = find_key_number(bln, {'violProb'});
base.meanScore  = find_key_number(bln, {'meanScore'});
cv_guess        = find_key_number(bln, {'CVaR','CVAR','cvar','cvarScore'}, true);
base.cvarScore  = cv_guess;  % keep both names for flexibility
base.cvar       = cv_guess;
end

function num = find_key_number(str, keys, use_last_if_missing)
if nargin<3, use_last_if_missing=false; end
num = NaN; found=false;
for k=1:numel(keys)
    key = regexptranslate('escape', keys{k});
    pat = ['(?i)\b' key '\b\s*[:=]\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)'];
    tok = regexp(str, pat, 'tokens', 'once');
    if ~isempty(tok)
        num = str2double(tok{1}); found=true; break;
    end
end
if ~found && use_last_if_missing
    nums = regexp(str, '([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)', 'tokens');
    if ~isempty(nums)
        num = str2double(nums{end}{1});
    end
end
end
