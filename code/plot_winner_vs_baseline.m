function plot_winner_vs_baseline()
% Compare baseline vs winner visually
fprintf('>>> plot_winner_vs_baseline: start\n');

T = readtable(fullfile('results','ranking14.csv'));
winner = T(1,:);

% Load scenario results (assume saved .mat from run_ieee14)
load(fullfile('results','scenario_data.mat'), 'baselineRes', 'winnerRes');

figure('Color','w');
subplot(1,2,1);
boxplot([baselineRes.lossMW, winnerRes.lossMW], ...
        'Labels', {'Baseline','Winnfunction plot_winner_vs_baseline
% PLOT_WINNER_VS_BASELINE – simple bar plot of baseline vs winner CVaR.

here     = fileparts(mfilename('fullpath'));
resDir   = fullfile(here, 'results');
figDir   = fullfile(here, 'figures');
if ~exist(figDir,'dir'), mkdir(figDir); end

csvPath  = fullfile(resDir, 'ranking14.csv');
basePath = fullfile(resDir, 'baseline_no_svc.txt');
assert(exist(csvPath, 'file')==2 && exist(basePath,'file')==2, ...
    'Missing results. Run run_ieee14 first.');

T   = readtable(csvPath, 'TextType','string');
[~,ix] = sort(T.cvarScore, 'ascend'); T = T(ix,:);
winner = T(1,:);

txt = fileread(basePath);
pat = 'CVaR\s*=\s*([0-9.+-Ee]+)';
tok = regexp(txt, pat, 'tokens', 'once');
assert(~isempty(tok), 'Could not parse baseline CVaR.');
baseCvar = str2double(tok{1});

vals = [baseCvar, winner.cvarScore];
labels = {'Baseline', char(winner.name)};

figure('Color','w','Name','CVaR: Baseline vs Winner');
bar(vals); grid on; xticklabels(labels); ylabel('CVaR (MW-equivalent)');
title('Tail Risk (CVaR_{0.9}) Comparison');
out = fullfile(figDir, 'winner_vs_baseline_cvar.png');
exportgraphics(gcf, out, 'Resolution', 200);
fprintf('>>> wrote %s\n', out);
end
er'});
ylabel('MW Loss');
title('Distribution of active power losses');

subplot(1,2,2);
boxplot([baselineRes.phiV, winnerRes.phiV], ...
        'Labels', {'Baseline','Winner'});
ylabel('Voltage penalty');
title('Distribution of voltage deviations');

outPNG = fullfile('figures','winner_vs_baseline.png');
if ~exist('figures','dir'), mkdir figures; end
exportgraphics(gcf, outPNG, 'Resolution', 200);
fprintf('>>> wrote %s\n', outPNG);
end
