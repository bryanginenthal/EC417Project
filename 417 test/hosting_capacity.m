% hosting_capacity.m
clc; clear; close all;

formatdata;          % bussys, bus, branch, nb, etc.
define_constants;

load('ren_profiles_case39.mat');   % Pwind_bus_MW, Psolar_bus_MW

% nb x T renewable matrix and peak-renewable hour
if size(Pwind_bus_MW,1) == nb
    P_ren_all = Pwind_bus_MW + Psolar_bus_MW;
else
    P_ren_all = (Pwind_bus_MW + Psolar_bus_MW).';
end

total_ren   = sum(P_ren_all,1);
[~, t_peak] = max(total_ren);
P_ren_base  = P_ren_all(:, t_peak);

fprintf('t_peak = %d, P_ren_peak = %.2f MW\n', t_peak, sum(P_ren_base));

% Line limits and PTDF
Fmax = branch(:, RATE_A);
Fmax(Fmax == 0) = 1e4;
PTDF = makePTDF(bussys);

% Base DC-OPF (reference flows)
base_res = dcopf(bussys);
if ~isfield(base_res,'success') || ~base_res.success
    error('Base DC-OPF did not converge.');
end
Pf0 = base_res.branch(:, PF);

% Injection pattern for scaling (net sum = 0)
slack_bus = 1;
delta_inj_pattern = P_ren_base;
delta_inj_pattern(slack_bus) = delta_inj_pattern(slack_bus) - sum(P_ren_base);

% Sweep renewable scale factor
alpha_list = linspace(0, 3, 13);
MaxLoading = nan(size(alpha_list));

for k = 1:length(alpha_list)
    delta_inj = alpha_list(k) * delta_inj_pattern;
    Pf_alpha  = Pf0 + PTDF * delta_inj;
    MaxLoading(k) = max(abs(Pf_alpha) ./ Fmax);
end

% Approximate alpha where max loading hits 1
idx_over = find(MaxLoading >= 1, 1, 'first');
if isempty(idx_over)
    alpha_host = alpha_list(end);
elseif idx_over == 1
    alpha_host = alpha_list(1);
else
    a1 = alpha_list(idx_over-1);
    a2 = alpha_list(idx_over);
    L1 = MaxLoading(idx_over-1);
    L2 = MaxLoading(idx_over);
    alpha_host = a1 + (1 - L1) * (a2 - a1) / (L2 - L1);
end

fprintf('alpha_host ≈ %.2f\n', alpha_host);

% Plot hosting curve
figure;
plot(alpha_list, MaxLoading, '-o', ...
    'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', 'auto');
hold on;
yline(1, '--', 'LineWidth', 1.2);
plot(alpha_host, 1, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
grid on;

xlabel('\alpha');
ylabel('Max line loading');
title('Hosting capacity (peak renewable)');
legend('Max loading','Limit','Host cap','Location','northwest');


