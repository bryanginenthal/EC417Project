%% DC OPF with PTDF for IEEE 39-bus system + LMPs
% Requirements:
%   - MATPOWER on the MATLAB path 
%   - Optimization Toolbox 

clear; clc; 


%% -------------------------------------------------
% Load case and basic dimensions
% --------------------------------------------------
define_constants;          % MATPOWER named indices (BUS_I, PD, etc.)
mpc = loadcase('case39');

bus    = mpc.bus;
branch = mpc.branch;
gen    = mpc.gen;
gencost= mpc.gencost;

nb = size(bus,    1);      % number of buses
nl = size(branch, 1);      % number of lines
ng = size(gen,    1);      % number of generators

%% -------------------------------------------------
% Loads and generators
% --------------------------------------------------

% Base loads (MW)
Pd0 = bus(:, PD);

% Optionally increase the load on the system to make congestion more likely
scale_vec = ones(nb,1);
heavy_buses = [15 16 18 21 23 24 25];  % pick some nodes of the system
scale_vec(heavy_buses) = 1;          % 
Pd = Pd0 .* scale_vec;                 % modified loads

% Generator data
gen_bus = gen(:, GEN_BUS);             % bus index for each generator
Pg_min  = gen(:, PMIN);                % MW
Pg_max  = gen(:, PMAX);                % MW

% Linear cost
c = [1 2 1.1 2 3 3.2 1 1 1 5].'; %c = ones(10,1);

%% -------------------------------------------------
% Line limits (Fmax)
% --------------------------------------------------

Fmax = branch(:, RATE_A);            % MW line limits from the case
Fmax(Fmax == 0) = 1e4;               % treat 0 as "very large" initially

% Tighten to force congestion 
Fmax(1) = 600;


%% -------------------------------------------------
% Build PTDF with MATPOWER
% --------------------------------------------------
% H = makePTDF(mpc) returns nl x nb matrix such that:
%      f (MW) = H * P_inj (MW),   sum(P_inj) = 0.
PTDF = makePTDF(mpc);                % nl x nb

%% -------------------------------------------------
% Bus-generator incidence matrix
% --------------------------------------------------

Cg = zeros(nb, ng);                  % nb x ng
for i = 1:ng
    Cg(gen_bus(i), i) = 1;
end

%% -------------------------------------------------
% DC-OPF constraints with PTDF formulation
%
% Cost: c*Pg
%
% inj = Cg*Pg - Pd       (bus injections, MW)
% f   = PTDF * inj       (line flows, MW)
%
% Constraints:
%   sum_b inj(b) = 0
%   -Fmax <= f <= Fmax
%   Pg_min <= Pg <= Pg_max
% --------------------------------------------------

% Power balance: sum_b (Cg*Pg - Pd) = 0  => (1^T * Cg) Pg = sum(Pd)
Aeq = (ones(1, nb) * Cg);            % 1 x ng
beq = sum(Pd);                       % scalar

% Line constraints:
% f = PTDF * (Cg*Pg - Pd)
%
% Upper: f <= Fmax
%   PTDF*Cg*Pg - PTDF*Pd <= Fmax
%   (PTDF*Cg) Pg <= Fmax + PTDF*Pd
A_line_pos = PTDF * Cg;              % nl x ng
b_line_pos = Fmax + PTDF * Pd;       % nl x 1

% Lower: -Fmax <= f  <=>  -f <= Fmax
%   -PTDF*(Cg*Pg - Pd) <= Fmax
%   -PTDF*Cg*Pg + PTDF*Pd <= Fmax
A_line_neg = -PTDF * Cg;             % nl x ng
b_line_neg = Fmax - PTDF * Pd;       % nl x 1

Aineq = [A_line_pos; A_line_neg];    % 2*nl x ng
bineq = [b_line_pos; b_line_neg];    % 2*nl x 1

% Generator bounds
lb = Pg_min;
ub = Pg_max;

%% -------------------------------------------------
% Solve LP with linprog
% --------------------------------------------------


options = optimoptions('linprog', ...
                       'Algorithm','interior-point', ...
                       'Display','iter');

[Pg_opt, cost_opt, exitflag, output, lambda] = ...
    linprog(c, Aineq, bineq, Aeq, beq, lb, ub, options);

if exitflag ~= 1
    warning('linprog did not converge to an optimal solution.');
end

%% -------------------------------------------------
% Recover flows, injections, LMPs
% --------------------------------------------------

inj_opt = Cg * Pg_opt - Pd;          % nb x 1
f_opt   = PTDF * inj_opt;            % nl x 1

% Duals from linprog
lambda_eq = lambda.eqlin;            % scalar (power balance)
mu        = lambda.ineqlin;          % 2*nl x 1 (line constraints)

mu_pos = mu(1:nl);                   % upper flow limits (f <= Fmax)
mu_neg = mu(nl+1:end);               % lower flow limits (-f <= Fmax)

% Correct LMP formula:
% LMP_k = d(opt cost) / d(Pd_k)
%       = -lambda_eq - sum_ell PTDF(ell,k)*(mu_pos(ell) - mu_neg(ell))
LMP = -lambda_eq * ones(nb,1) - PTDF.' * (mu_pos - mu_neg);

%% -------------------------------------------------
% Diagnostics: binding constraints and unique LMPs
% --------------------------------------------------

tol_bind = 1e-5;
Ax = Aineq * Pg_opt;

is_bind = abs(Ax - bineq) <= tol_bind;
bind_idx = find(is_bind);

fprintf('\n===== Binding inequality constraints (within %.1e) =====\n', tol_bind);
if isempty(bind_idx)
    fprintf('  None.\n');
else
    for k = 1:length(bind_idx)
        j = bind_idx(k);
        if j <= nl
            ell = j; which = 'upper';
        else
            ell = j - nl; which = 'lower';
        end
        fprintf('  Aineq(%3d,:)  line %3d  %-5s  dual = %.6f,  f = %8.2f,  limit = %8.2f\n', ...
            j, ell, which, lambda.ineqlin(j), f_opt(ell), Fmax(ell));
    end
end

fprintf('\nUnique LMP values (rounded to 4 decimals):\n');
disp(unique(round(LMP,4)));

fprintf('\nTotal cost: %.2f $/h\n', cost_opt);

fprintf('\nOptimal generation:\n');
for i = 1:ng
    fprintf('  Gen %2d at bus %2d: %8.2f MW\n', i, gen_bus(i), Pg_opt(i));
end

fprintf('\nLine flows (first 10):\n');
for ell = 1:nl
    fprintf('  Line %2d (%2d -> %2d): %8.2f MW (limit +/- %.1f MW)\n', ...
        ell, branch(ell, F_BUS), branch(ell, T_BUS), f_opt(ell), Fmax(ell));
end



%% -------------------------------------------------
% Plot LMPs vs bus index
% -------------------------------------------------
figure; 
bar(1:nb, LMP);
grid on;
xlabel('Bus index');
ylabel('LMP [$/MWh]');
title('Locational Marginal Prices (IEEE 39-bus)');
set(gca, 'XTick', 1:nb);

%% -------------------------------------------------
% Plot network: buses (LMPs), flows, and congestion
% -------------------------------------------------

figure; clf; hold on; box on;

% Build graph using from/to buses
line_from = branch(:, F_BUS);
line_to   = branch(:, T_BUS);

G = graph(line_from, line_to);   % undirected graph

hG = plot(G, 'Layout', 'force', 'NodeLabel', []);
bus_x = hG.XData;
bus_y = hG.YData;
delete(hG);      % remove default graph plot, keep coordinates
hold on;

% Parameters for congestion visualization
cong_threshold = 1;   % consider congested if |f| >= 99.9% of limit

% Plot lines with flows and congestion highlighting
for ell = 1:nl
    b_from = line_from(ell);
    b_to   = line_to(ell);

    x1 = bus_x(b_from);  y1 = bus_y(b_from);
    x2 = bus_x(b_to);    y2 = bus_y(b_to);

    flow  = f_opt(ell);      % MW
    limit = Fmax(ell);       % MW

    if limit <= 0
        loading = 0;         % treat as unconstrained
    else
        loading = abs(flow) / limit;
    end

    % Style: red & thick if congested, gray otherwise
    if loading >= cong_threshold
        lineColor = 'r';           % congested
        lineWidth = 3;
    else
        lineColor = [0.4 0.4 0.4]; % non-congested
        lineWidth = 1.5;
    end

    % Draw the line
    plot([x1, x2], [y1, y2], '-', ...
        'Color', lineColor, 'LineWidth', lineWidth);

    % Midpoint for annotation
    xm = 0.5*(x1 + x2);
    ym = 0.5*(y1 + y2);

    % Annotate flow
    text(xm, ym, sprintf('%.1f MW', flow), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 7);

    % Direction arrow (positive from -> to)
    dx = x2 - x1;
    dy = y2 - y1;
    if abs(dx) + abs(dy) > 1e-6
        net_scale   = max(max(bus_x) - min(bus_x), max(bus_y) - min(bus_y));
        arrow_scale = 0.04 * net_scale;

        if flow >= 0
            xa = xm - 0.3*dx;
            ya = ym - 0.3*dy;
            quiver(xa, ya, arrow_scale*dx, arrow_scale*dy, 0, ...
                   'MaxHeadSize', 1.5, 'LineWidth', 0.8);
        else
            xa = xm + 0.3*dx;
            ya = ym + 0.3*dy;
            quiver(xa, ya, -arrow_scale*dx, -arrow_scale*dy, 0, ...
                   'MaxHeadSize', 1.5, 'LineWidth', 0.8);
        end
    end
end

% Plot buses colored by LMP
scatter(bus_x, bus_y, 80, LMP, 'filled');
colormap(jet);
cb = colorbar;
ylabel(cb, 'LMP [$/MWh]');

% Color limits for better contrast
cmin = min(LMP);
cmax = max(LMP);
if abs(cmax - cmin) < 1e-3
    caxis([cmin - 0.01, cmax + 0.01]);
else
    caxis([cmin, cmax]);
end

% Annotate buses with ID and LMP
for k = 1:nb
    text(bus_x(k) + 0.01, bus_y(k) + 0.01, ...
        sprintf('%d\n%.2f', k, LMP(k)), ...
        'FontSize', 7, 'FontWeight', 'bold');
end

% Legend for congested vs non-congested lines
h1 = plot(nan, nan, 'r-',  'LineWidth', 3);
h2 = plot(nan, nan, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
legend([h1, h2], {'Congested line', 'Non-congested line'}, ...
       'Location', 'bestoutside');

title('IEEE 39-bus DC-OPF: LMPs, Power Flows, and Congestion (graph layout)');
axis equal;
axis off;
hold off;
