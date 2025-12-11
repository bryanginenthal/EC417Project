clc; clear; close all;

% Using a reduced dataset here purely to speed up run time.
% Only the essential comments are kept.

formatdata;        % defines bussys, bus, branch, gen, nb, nl, ng
define_constants;

Pd0      = bus(:, PD);
Pd_total = sum(Pd0);

Fmax = branch(:, RATE_A);
Fmax(Fmax == 0) = 1e4;

PTDF = makePTDF(bussys);

% Generator data
N       = ng;
gen_bus = gen(:, GEN_BUS);
Pmin    = gen(:, PMIN);
Pmax    = gen(:, PMAX);

% Simple linear cost
c = linspace(10, 20, N);

% Bus–generator incidence
Cg = zeros(nb, N);
for i = 1:N
    Cg(gen_bus(i), i) = 1;
end
H = PTDF * Cg;

% Simple renewable model
ren_bus      = 30;
base_ren_MW  = 0.3 * Pd_total;
alpha_list   = linspace(0, 2, 9);

Cost_simple = zeros(size(alpha_list));
Cost_ptdf   = zeros(size(alpha_list));
CongCost    = zeros(size(alpha_list));
NumCong     = zeros(size(alpha_list));

Pg_simple_store = zeros(N, length(alpha_list));
Pg_ptdf_store   = zeros(N, length(alpha_list));

opts = optimoptions('linprog','Display','none');

% ---------------------------------------------------------
% Main sweep over renewable penetration
% ---------------------------------------------------------
for k = 1:length(alpha_list)

    alpha = alpha_list(k);
    inj_ren         = zeros(nb,1);
    inj_ren(ren_bus)= alpha * base_ren_MW;

    Pd_eff = Pd_total - sum(inj_ren);

    % ---------------- Simple dispatch (no line limits)
    f   = c';
    Aeq = ones(1, N);
    beq = Pd_eff;

    [Pg_s, ~, flag_s] = linprog(f, [], [], Aeq, beq, Pmin, Pmax, opts);

    if flag_s == 1
        Pg_simple_store(:,k) = Pg_s;
        Cost_simple(k)       = sum(c .* Pg_s');
    else
        Pg_simple_store(:,k) = nan;
        Cost_simple(k)       = nan;
    end

    % ---------------- PTDF-based DC-OPF
    rhs_shift = PTDF * (inj_ren - Pd0);

    Aineq = [ H;
             -H ];
    bineq = [ Fmax - rhs_shift;
              Fmax + rhs_shift ];

    [Pg_c, ~, flag_c] = linprog(f, Aineq, bineq, Aeq, beq, Pmin, Pmax, opts);

    if flag_c == 1
        Pg_ptdf_store(:,k) = Pg_c;
        Cost_ptdf(k)       = sum(c .* Pg_c');

        flows    = H*Pg_c + rhs_shift;
        loading  = abs(flows)./Fmax;
        NumCong(k)= sum(loading >= 1 - 1e-6);
    else
        Pg_ptdf_store(:,k) = nan;
        Cost_ptdf(k)       = nan;
        NumCong(k)         = nan;
    end

    CongCost(k) = Cost_ptdf(k) - Cost_simple(k);
end

% ---------------------------------------------------------
% Plots (Tasks 7 & 8)
% ---------------------------------------------------------
lw = 1.8; ms = 7;

% Total cost comparison
figure;
plot(alpha_list, Cost_simple,'-o','LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','auto'); hold on;
plot(alpha_list, Cost_ptdf,  '-o','LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','auto');
grid on;
xlabel('Renewable Penetration \alpha');
ylabel('Total Cost');
title('Total Cost vs Renewable Penetration');
legend('Simple','PTDF-OPF','Location','northwest');

% Congestion cost
figure;
plot(alpha_list, CongCost,'-o','LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','auto');
grid on;
xlabel('\alpha');
ylabel('Congestion Cost');
title('Congestion Cost vs Renewable Penetration');

% Number of congested lines
figure;
plot(alpha_list, NumCong,'-o','LineWidth',lw,'MarkerSize',ms,'MarkerFaceColor','auto');
grid on;
xlabel('\alpha');
ylabel('Congested Lines');
title('Line Congestion vs Renewable Penetration');

% Generator dispatch comparison at max alpha
[~, idx] = max(alpha_list);
figure;
bar([Pg_simple_store(:,idx), Pg_ptdf_store(:,idx)]);
grid on;
xlabel('Generator');
ylabel('MW');
title(sprintf('Dispatch Comparison at \\alpha = %.2f', alpha_list(idx)));
legend('Simple','PTDF','Location','northwest');

