clc; clear; close all;

formatdata;
load('ren_profiles_case39.mat');

% Choose one hour in the 24-hour window
t = 1;

% Base renewable injection at that hour (MW)
P_ren_base = Pwind_bus_MW(:, t) + Psolar_bus_MW(:, t);

% Line limits
Fmax = branch(:, RATE_A);
Fmax(Fmax == 0) = 1e4;

% Renewable penetration sweep
alpha_list = linspace(0, 3, 13);

Feasible   = zeros(size(alpha_list));
MaxLoading = zeros(size(alpha_list));

for k = 1:length(alpha_list)

    alpha = alpha_list(k);

    % Scaled renewables
    P_ren = alpha * P_ren_base;

    % Modify case (treat renewables as negative load)
    mpc = bussys;
    mpc.bus(:, PD) = mpc.bus(:, PD) - P_ren;

    % Run DC-OPF
    try
        results = dcopf(mpc);
    catch
        results.success = 0;
    end

    if isfield(results, 'success') && results.success
        Feasible(k) = 1;

        % Branch flows (MW)
        Pf = results.branch(:, PF);
        Pt = results.branch(:, PT);

        flow = max(abs(Pf), abs(Pt));
        loading = flow ./ Fmax;

        MaxLoading(k) = max(loading);
    else
        Feasible(k) = 0;
        MaxLoading(k) = nan;
    end

end

% Hosting capacity estimate
idx_fail = find(Feasible == 0, 1, 'first');

if isempty(idx_fail)
    alpha_host = alpha_list(end);
else
    alpha_host = alpha_list(max(1, idx_fail-1));
end

fprintf('Estimated hosting capacity alpha ~ %.2f (DC-OPF feasibility).\n', alpha_host);

% Plot feasibility curve
figure;
plot(alpha_list, Feasible, '-o', 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'auto');
xlabel('Renewable Penetration Scale \alpha');
ylabel('Feasibility (1 = Feasible)');
title('Task 5: Hosting Capacity Feasibility (DC-OPF)');
grid on;
ylim([-0.1 1.1]);

% Plot max loading curve (only where feasible)
figure;
plot(alpha_list, MaxLoading, '-o', 'LineWidth', 1.8, ...
    'MarkerSize', 7, 'MarkerFaceColor', 'auto');
yline(1, '--', 'LineWidth', 1.2);
xlabel('Renewable Penetration Scale \alpha');
ylabel('Max Line Loading (|Flow| / Limit)');
title('Task 5: Hosting Capacity Curve (DC-OPF Results)');
grid on;
