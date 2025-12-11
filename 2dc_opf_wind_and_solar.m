clearvars
clc
formatdata;

%% -------------------------------------------------
% Loads and generators
% --------------------------------------------------

% Base loads (MW)
Pd0 = bus(:, PD);


% Generator data
gen_bus = gen(:, GEN_BUS);             % bus index for each generator
Pg_min  = gen(:, PMIN);                % MW
Pg_max  = gen(:, PMAX);                % MW


% Costs from IEEE-39 Bus system - polynomical costs requiring 3 vars
%                                 [0.01 0.3 0.2] for each generator
c = ones(ng,3);
c(:,1) = 0.01;
c(:,2) = .3;
c(:,3) = .2;

% Treating solar/wind as "free" to run
%% -------------------------------------------------
% Line limits (Fmax)
% --------------------------------------------------

Fmax = branch(:, RATE_A);            % MW line limits from the case
Fmax(Fmax == 0) = 1e4;               % treat 0 as "very large" initially

%% -------------------------------------------------
% Build PTDF with MATPOWER
% --------------------------------------------------
% H = makePTDF(bussys) returns nl x nb matrix such that:
%      f (MW) = H * P_inj (MW),   sum(P_inj) = 0.
PTDF = makePTDF(bussys);                % nl x nb

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
% Cost: c2*Pg^2 + c1*Pg + c0
%
% inj = Cg*Pg + Ps + Pw - Pd - addload   (bus injections, MW)
%
% f   = PTDF * inj       (line flows, MW)
%
% Constraints:
%   sum_b inj(b) = 0
%   -Fmax <= f <= Fmax
%   Pg_min <= Pg <= Pg_max
%   0 <= zs <= 1
%   0 <= zw <= 1
%   0 <= Ps <= solargen * zs
%   0 <= Pw <= windgen * zw
% --------------------------------------------------

Cons = [];
% vars - Pg1-Pg10, Psolar1-Psolar39, Pwind1-Pwind39, zsolar1-zsolar39,
% zwind1-zwind39
pg = sdpvar(ng,1);
psolar = sdpvar(nb,1);
pwind = sdpvar(nb,1);
zsolar = binvar(nb,1);
zwind = binvar(nb,1);
x = [pg;psolar;pwind;zsolar;zwind];

% Power balance: sum_b (Cg*Pg+Ps+Pw- Pd) = 0
% => (1^T * Cg) Pg + (Ps + Pw) = sum(Pd)
A_peq = [(ones(1, nb) * Cg), ones(1,nb), ones(1,nb), zeros(1,nb), zeros(1,nb)];            % 1 x ng + 4 * nb
b_peq = sum(Pd0);                % scalar

Cons = [Cons; A_peq*x == b_peq];

% zs + zw = 2
%A_zs = [zeros(1,ng), zeros(1,nb), zeros(1,nb), ones(1,nb), ones(1,nb)]; % 1 x ng + 4*nb
%b_zs = 2;   % scalar


% zw = 1
%A_zw = [zeros(1,ng), zeros(1,nb), zeros(1,nb), zeros(1,nb), ones(1,nb)]; % 1 x ng + 4*nb
%b_zw = 1;   % scalar

Aeq = [A_peq]; %;A_zs]; %;A_zw];
beq = [b_peq]; %;b_zs]; %;b_zw];

% Line constraints:
% f = PTDF * (Cg*Pg + Ps + Pw - Pd)
%
% Upper: f <= Fmax
%   PTDF*(Cg*Pg+Ps+Pw - Pd0) <= Fmax
%   (PTDF*Cg) Pg + (PTDF*Ps) + (PTDF*Pw)  <= Fmax + PTDF*(Pd0)

A_line_pos = [PTDF * Cg, PTDF, PTDF, zeros(nl,nb), zeros(nl,nb)] ;      % nl x ng + 4*nb
b_line_pos = Fmax + PTDF * (Pd0);       % nl x 1
                                                                        
% Lower: -Fmax <= f  <=>  -f <= Fmax
%   -PTDF*(Cg*Pg - Pd) <= Fmax
%   -PTDF*Cg*Pg + PTDF*Pd <= Fmax
A_line_neg = [-PTDF * Cg, -PTDF, -PTDF, zeros(nl,nb), zeros(nl,nb)];             % nl x ng + 4*nb
b_line_neg = Fmax - PTDF * Pd0;       % nl x 1

Cons = [Cons; A_line_pos*x <= b_line_pos];
Cons = [Cons; A_line_neg*x <= b_line_neg];

%   0 <= Ps <= solargen * zs
%   Ps - solargen*zs <= 0
A_solar = [zeros(length(solargen),ng), ones(length(solargen),nb), zeros(length(solargen),nb), -solargen.*ones(length(solargen),nb), zeros(length(solargen),nb)];          % #scenarios(length(solargen)) x ng + 4*nb
%   0 <= Pw <= windgen * zw
%   Pw - windgen*zw <= 0
A_wind = [zeros(length(windgen),ng), zeros(length(windgen),nb), ones(length(windgen),nb), zeros(length(windgen),nb), -windgen.*ones(length(windgen),nb)];            % #scenarios(length(windgen)) x ng + 4*nb

Aineq = [A_line_pos; A_line_neg; A_solar; A_wind];    % 2*nl + 2 x ng+4*nb
bineq = [b_line_pos; b_line_neg; zeros(length(solargen),1); zeros(length(windgen),1)];    % 2*nl + 2*#scenarios x 1

%Cons = [Cons; A_solar*x <= zeros(length(solargen),1)];
%Cons = [Cons; A_wind*x <= zeros(length(windgen),1)];

% Generator bounds
%   Pg_min <= Pg <= Pg_max
Cons = [Cons; Pg_min <= pg];
Cons = [Cons; pg <= Pg_max];
Cons = [Cons; 0 <= psolar];
for i = 1:length(solargen)
    Cons = [Cons; psolar <= solargen(i).*zsolar];
end

Cons = [Cons; 0 <= pwind];
for i = 1:length(windgen)
    Cons = [Cons; pwind <= windgen(i).*zwind];
end

Cons = [Cons; sum(zsolar) == 1];
Cons = [Cons; sum(zwind) == 1];


% ---- YALMIP variables and constraints
%x1 = sdpvar(1,1); x2 = sdpvar(1,1);
zs  = sdpvar(nb,1);
zw = sdpvar(nb,1);
%% -- Solve cost with intlinprog

c = [ones(1,ng), zeros(1,4*nb)]; % linear cost for generators, no other costs considered
obj = c * x;
opts = sdpsettings('solver','intlinprog','verbose',0);
diagsln = optimize(Cons, obj, opts);

% ---- Report
code = diagsln.problem; % 0=OK, 1=infeasible, 2=unbounded, 3=inaccurate, 4=numerical
if code ~= 0
    fprintf('YALMIP/solver status: problem=%d, info="%s"\n', code, diagsln.info);
else
    fprintf('Optimal solution found.\n');
end

x_opt = value(x);
fval  = value(obj);

gensol = x_opt(1:ng);
solarsol = x_opt(ng+1:ng+nb);
windsol = x_opt(ng+nb+1:ng+2*nb);
solarsel = x_opt(ng+2*nb+1:ng+3*nb);
windsel = x_opt(ng+3*nb+1:end);

windbus = find(windsel~=0);
solarbus = find(solarsel~=0);
fprintf('Power at each generator:\n')
disp(gensol)
fprintf('Wind generator placed at bus %1d producing %.5f MW.\n', windbus, windsol(windbus)*windcap)

fprintf('Solar generator placed at bus %1d producing %.5f MW.\n', solarbus, solarsol(solarbus)*solarcap)


%%
bussys.bus(solarbus,3) = bussys.bus(solarbus,3) - solarsol(solarbus)*solarcap;
bussys.bus(windbus,3) = bussys.bus(windbus,3) - windsol(windbus)*windcap;
dcopf(bussys)


%%
inj_opt = Cg * gensol - bussys.bus(:,3);          % nb x 1
f_opt   = PTDF * inj_opt;   
LMP = ones(nb,1) * 10.947;
%% -------------------------------------------------
% Plot network: buses (LMPs), flows, and congestion
% -------------------------------------------------

% Taken from classwork example code - 'DC_opf_with_lmps_39bus_system.m'


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

title('Adjusted IEEE 39-bus DC-OPF: LMPs, Power Flows, and Congestion (graph layout)');
axis equal;
axis off;
hold off;