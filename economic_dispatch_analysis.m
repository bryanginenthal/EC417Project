clc; clear; close all;

formatdata;
load('peak_placement_result.mat');   % windbus, solarbus

% Number of generators
N = ng;

% Linearized cost function coefficients (approximated)
c = linspace(1, 2, N);  % Simple merit-order proxy

% Generator limits (Pmin and Pmax)
Pmin = gen(:, PMIN).';
Pmax = gen(:, PMAX).';

% Total load demand
Pd0 = bus(:, PD);
Pd_total = sum(Pd0);

% PTDF model data
Fmax = branch(:, RATE_A);
Fmax(Fmax == 0) = 1e4;

PTDF = makePTDF(bussys);

gen_bus = gen(:, GEN_BUS);
Cg = zeros(nb, N);
for i = 1:N
    Cg(gen_bus(i), i) = 1;
end

H = PTDF * Cg;

% Evaluation hour
if ~exist('startIdx24','var')
    startIdx24 = 1;
end
t_eval = startIdx24;

% Renewable penetration sweep
alpha_list = linspace(0, 2, 9);

% Storage
Cost_simple = zeros(size(alpha_list));
Cost_ptdf   = zeros(size(alpha_list));
CongCost    = zeros(size(alpha_list));
NumCong     = zeros(size(alpha_list));

Pg_simple_store = zeros(N, length(alpha_list));
Pg_ptdf_store   = zeros(N, length(alpha_list));

opts = optimoptions('linprog','Display','none');

for k = 1:length(alpha_list)

    alpha = alpha_list(k);

    % Renewable injections (MW)
    load('ren_profiles_case39.mat');
    t = 1;  % or your chosen hour in the 24h window

    P_ren_base = Pwind_bus_MW(:, t) + Psolar_bus_MW(:, t);
    P_ren = alpha * P_ren_base;


    inj_ren = zeros(nb,1);
    inj_ren(solarbus) = Ps;
    inj_ren(windbus)  = Pw;

    % Total load demand
    Pd = Pd_total - sum(inj_ren);

    % Decision variable: Power output P = [P1, P2, ..., PN]
    f = c';

    % Constraints
    Aeq = ones(1, N);
    beq = Pd;

    lb = Pmin';
    ub = Pmax';

    % Simple Dispatch (No line constraints)
    A = [];
    b_ineq = [];

    [P_opt_simple, ~, exitflag_s] = linprog(f, A, b_ineq, Aeq, beq, lb, ub, opts);

    if exitflag_s == 1
        Pg_simple_store(:,k) = P_opt_simple;
        Cost_simple(k) = sum(c .* P_opt_simple');
    else
        Pg_simple_store(:,k) = nan(N,1);
        Cost_simple(k) = nan;
    end

    % PTDF-based DC-OPF (With line constraints)
    rhs_shift = PTDF * (inj_ren - Pd0);

    A = [ H;
         -H ];

    b_ineq = [ Fmax - rhs_shift;
               Fmax + rhs_shift ];

    [P_opt_ptdf, ~, exitflag_c] = linprog(f, A, b_ineq, Aeq, beq, lb, ub, opts);

    if exitflag_c == 1
        Pg_ptdf_store(:,k) = P_opt_ptdf;
        Cost_ptdf(k) = sum(c .* P_opt_ptdf');

        f_line = H * P_opt_ptdf + rhs_shift;
        loading = abs(f_line) ./ Fmax;
        NumCong(k) = sum(loading >= 1 - 1e-6);
    else
        Pg_ptdf_store(:,k) = nan(N,1);
        Cost_ptdf(k) = nan;
        NumCong(k) = nan;
    end

    % Congestion cost proxy
    CongCost(k) = Cost_ptdf(k) - Cost_simple(k);

end

% Display results
Result = table(alpha_list(:), Cost_simple(:), Cost_ptdf(:), CongCost(:), NumCong(:), ...
    'VariableNames', {'Alpha','Cost_Simple','Cost_PTDF','Congestion_Cost','Num_Congested_Lines'});

disp(Result);

%% ---------------- Plot Styling Helpers ----------------
lw = 1.8;
ms = 7;

%% Task 8: Simple vs PTDF Cost vs Renewable Penetration
figure;
plot(alpha_list, Cost_simple, '-o', ...
    'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', 'auto'); 
hold on;
plot(alpha_list, Cost_ptdf, '-o', ...
    'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', 'auto');

xlabel('Renewable Penetration Scale \alpha');
ylabel('Total Linear Cost (Arb. Units)');
title('Task 8: Simple vs PTDF Cost vs Renewable Penetration');
legend('Simple Dispatch','PTDF-based DC-OPF','Location','northeast');
grid on;

%% Task 7: Congestion Cost vs Renewable Penetration
figure;
plot(alpha_list, CongCost, '-o', ...
    'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', 'auto');

xlabel('Renewable Penetration Scale \alpha');
ylabel('Congestion Cost (Arb. Units)');
title('Task 7: Congestion Cost vs Renewable Penetration');
grid on;

%% Task 7: Congested Lines vs Renewable Penetration
figure;
plot(alpha_list, NumCong, '-o', ...
    'LineWidth', lw, 'MarkerSize', ms, 'MarkerFaceColor', 'auto');

xlabel('Renewable Penetration Scale \alpha');
ylabel('Number Of Congested Lines');
title('Task 7: Congested Lines vs Renewable Penetration');
grid on;

%% Optional: Dispatch comparison at highest alpha
[~, kmax] = max(alpha_list);

figure;
bar([Pg_simple_store(:,kmax), Pg_ptdf_store(:,kmax)]);
xlabel('Generator Index');
ylabel('MW');
title(sprintf('Dispatch Comparison at \\alpha = %.2f', alpha_list(kmax)));
legend('Simple','PTDF','Location','northeast');
grid on;

