function inspect_tail_scenario()
% Inspect one extreme scenario for intuition
fprintf('>>> inspect_tail_scenario: start\n');

load(fullfile('results','scenario_data.mat'), 'scenarios', 'winnerRes');

% Find worst-case scenario index for winner (max loss)
[~, idx] = max(winnerRes.lossMW);

fprintf('Worst-case scenario #%d:\n', idx);
disp(scenarios(idx));

fprintf('Winner results:\n');
disp(struct('V_abs', winnerRes.V_abs(idx,:), ...
            'totPloss', winnerRes.lossMW(idx), ...
            'phiV', winnerRes.phiV(idx)));
end
