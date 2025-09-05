function summarize_for_paper
% Robust summary generator (paper-ready):
% - Reads results/ranking14.csv
% - Tries parsing results/baseline_no_svc.txt
% - If parsing fails, recomputes baseline with same seed/Nmc
% - Prints improvements and writes results/summary.txt

clc; fprintf('>>> summarize_for_paper: start\n');

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
csvPath  = fullfile(resDir, 'ranking14.csv');
basePath = fullfile(resDir, 'baseline_no_svc.txt');

assert(exist(csvPath,'file')==2, 'Missing %s. Run run_ieee14 first.', csvPath);

% read candidate ranking & pick winner by CVaR
T = readtable(csvPath, 'TextType','string');
req = ["name","meanLoss","meanPhi","violProb","meanScore","cvarScore"];
for j = 1:numel(req)
    assert(any(strcmpi(T.Properties.VariableNames, req(j))), ...
        'ranking14.csv missing column: %s', req(j));
end
[~,ix] = sort(T.cvarScore, 'ascend');
T = T(ix,:);
W = T(1,:);

% --- baseline: parse if possible, else recompute
haveBaseline = false;
base = struct('meanLoss',NaN,'meanPhi',NaN,'violProb',NaN,'meanScore',NaN,'cvar',NaN);

if exist(basePath,'file')==2
    try
        txt = fileread(basePath);
        baseParsed = try_parse_baseline(txt);
        % normalize field name cvar/cvarScore
        if isfield(baseParsed,'cvar'); base.cvar = baseParsed.cvar; end
        if isfield(baseParsed,'cvarScore'); base.cvar = baseParsed.cvarScore; end
        base.meanLoss  = baseParsed.meanLoss;
        base.meanPhi   = baseParsed.meanPhi;
        base.violProb  = baseParsed.violProb;
        base.meanScore = baseParsed.meanScore;
        haveBaseline = all(isfinite(struct2array(base)));
        if ~haveBaseline
            warning('Parser could not extract all baseline fields from %s. Will recompute.', basePath);
        end
    catch ME
        warning('Failed reading/parsing %s (%s). Will recompute baseline.', basePath, ME.message);
    end
else
    warning('Baseline file not found: %s. Will recompute baseline.', basePath);
end

if ~haveBaseline
    fprintf('>>> recomputing baseline (no SVC) ... ');
    rng(42,'twister');
    C  = config14();
    Nmc = 500;
    S  = build_scenarios(C, Nmc, 'p_outage', 0.05, 'timeMode', "random");
    B  = baseline_no_svc(C, S, 100, 0.90);   % expected fields may be cvar or cvarScore
    base.meanLoss  = getf(B, ["meanLoss"]);
    base.meanPhi   = getf(B, ["meanPhi"]);
    base.violProb  = getf(B, ["violProb"]);
    base.meanScore = getf(B, ["meanScore"]);
    base.cvar      = getf(B, ["cvar","cvarScore"]);
    fprintf('ok\n');
end

%  improvements (winner vs baseline)
imp.cvarAbs   = base.cvar - W.cvarScore;
imp.cvarRel   = 100 * imp.cvarAbs / base.cvar;
imp.meanAbs   = base.meanScore - W.meanScore;
imp.meanRel   = 100 * imp.meanAbs / base.meanScore;
imp.violAbs   = base.violProb - W.violProb;
imp.violRel   = 100 * imp.violAbs / max(base.violProb, eps);

% print to console
fprintf('\n=== Baseline (no SVC) ===\n');
fprintf('meanLoss = %.3f MW\n', base.meanLoss);
fprintf('meanPhi  = %.5f pu-sum\n', base.meanPhi);
fprintf('violProb = %.3f\n', base.violProb);
fprintf('meanScore= %.3f\n', base.meanScore);
fprintf('CVaR     = %.3f\n', base.cvar);

fprintf('\n=== Winner (by CVaR) ===\n');
fprintf('name     = %s\n', W.name);
fprintf('meanLoss = %.3f MW\n', W.meanLoss);
fprintf('meanPhi  = %.5f pu-sum\n', W.meanPhi);
fprintf('violProb = %.3f\n', W.violProb);
fprintf('meanScore= %.3f\n', W.meanScore);
fprintf('CVaR     = %.3f\n', W.cvarScore);

fprintf('\n=== Improvements (Winner vs Baseline) ===\n');
fprintf('CVaR     : %+ .3f  (%+.1f%%)\n', imp.cvarAbs, imp.cvarRel);
fprintf('MeanScore: %+ .3f  (%+.1f%%)\n', imp.meanAbs, imp.meanRel);
fprintf('ViolProb : %+ .3f  (%+.1f%%)\n', imp.violAbs, imp.violRel);

% --- write summary
outTxt = fullfile(resDir, 'summary.txt');
fid = fopen(outTxt, 'w');  assert(fid>0, 'Cannot write %s', outTxt);
fprintf(fid, 'Baseline (no SVC): meanLoss=%.3f MW, meanPhi=%.5f, violProb=%.3f, meanScore=%.3f, CVaR=%.3f\n', ...
    base.meanLoss, base.meanPhi, base.violProb, base.meanScore, base.cvar);
fprintf(fid, 'Winner (by CVaR): %s | meanLoss=%.3f MW, meanPhi=%.5f, violProb=%.3f, meanScore=%.3f, CVaR=%.3f\n', ...
    W.name, W.meanLoss, W.meanPhi, W.violProb, W.meanScore, W.cvarScore);
fprintf(fid, 'Improvements (Winner vs Baseline): CVaR=%+.3f (%+.1f%%), MeanScore=%+.3f (%+.1f%%), ViolProb=%+.3f (%+.1f%%)\n', ...
    imp.cvarAbs, imp.cvarRel, imp.meanAbs, imp.meanRel, imp.violAbs, imp.violRel);
fclose(fid);

fprintf('\n>>> wrote %s\n', outTxt);
fprintf('>>> summarize_for_paper: done\n');
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


function v = getf(s, names)
% Return first existing finite field s.(names{k})
v = NaN;
for k=1:numel(names)
    if isfield(s, names{k}) && isfinite(s.(names{k}))
        v = s.(names{k}); return;
    end
end
end
