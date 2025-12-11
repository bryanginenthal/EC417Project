clearvars
clc
formatdata;

Pd0 = bus(:, PD);

gen_bus = gen(:, GEN_BUS);
Pg_min  = gen(:, PMIN);
Pg_max  = gen(:, PMAX);

c = ones(ng,3);
c(:,1) = 0.01;
c(:,2) = .3;
c(:,3) = .2;

Fmax = branch(:, RATE_A);
Fmax(Fmax == 0) = 1e4;

PTDF = makePTDF(bussys);

Cg = zeros(nb, ng);
for i = 1:ng
    Cg(gen_bus(i), i) = 1;
end

Cons = [];

pg = sdpvar(ng,1);
psolar = sdpvar(nb,1);
pwind = sdpvar(nb,1);
zsolar = binvar(nb,1);
zwind = binvar(nb,1);
x = [pg;psolar;pwind;zsolar;zwind];

A_peq = [(ones(1, nb) * Cg), ones(1,nb), ones(1,nb), zeros(1,nb), zeros(1,nb)];
b_peq = sum(Pd0);
Cons = [Cons; A_peq*x == b_peq];

A_line_pos = [PTDF * Cg, PTDF, PTDF, zeros(nl,nb), zeros(nl,nb)] ;
b_line_pos = Fmax + PTDF * (Pd0);

A_line_neg = [-PTDF * Cg, -PTDF, -PTDF, zeros(nl,nb), zeros(nl,nb)];
b_line_neg = Fmax - PTDF * Pd0;

Cons = [Cons; A_line_pos*x <= b_line_pos];
Cons = [Cons; A_line_neg*x <= b_line_neg];

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

c = [ones(1,ng), zeros(1,4*nb)];
obj = c * x;

opts = sdpsettings('solver','intlinprog','verbose',0);
diagsln = optimize(Cons, obj, opts);

code = diagsln.problem;
if code ~= 0
    fprintf('YALMIP/solver status: problem=%d, info="%s"\n', code, diagsln.info);
else
    fprintf('Optimal solution found.\n');
end

x_opt = value(x);

gensol   = x_opt(1:ng);
solarsol = x_opt(ng+1:ng+nb);
windsol  = x_opt(ng+nb+1:ng+2*nb);
solarsel = x_opt(ng+2*nb+1:ng+3*nb);
windsel  = x_opt(ng+3*nb+1:end);

windbus  = find(windsel~=0);
solarbus = find(solarsel~=0);

fprintf('Power at each generator:\n')
disp(gensol)

fprintf('Wind generator placed at bus %1d producing %.5f MW.\n', ...
    windbus, windsol(windbus)*windcap)

fprintf('Solar generator placed at bus %1d producing %.5f MW.\n', ...
    solarbus, solarsol(solarbus)*solarcap)

%%
bussys.bus(solarbus,3) = bussys.bus(solarbus,3) - solarsol(solarbus)*solarcap;
bussys.bus(windbus,3)  = bussys.bus(windbus,3)  - windsol(windbus)*windcap;
dcopf(bussys)

%%
inj_opt = Cg * gensol - bussys.bus(:,3);
f_opt   = PTDF * inj_opt;
LMP = ones(nb,1) * 10.947;

%% -------------------------------------------------
% LMP + congestion graph (as in your original)
%% -------------------------------------------------
figure; clf; hold on; box on;

line_from = branch(:, F_BUS);
line_to   = branch(:, T_BUS);

G = graph(line_from, line_to);

hG = plot(G, 'Layout', 'force', 'NodeLabel', []);
bus_x = hG.XData;
bus_y = hG.YData;
delete(hG);
hold on;

cong_threshold = 1;

for ell = 1:nl
    b_from = line_from(ell);
    b_to   = line_to(ell);

    x1 = bus_x(b_from);  y1 = bus_y(b_from);
    x2 = bus_x(b_to);    y2 = bus_y(b_to);

    flow  = f_opt(ell);
    limit = Fmax(ell);

    if limit <= 0
        loading = 0;
    else
        loading = abs(flow) / limit;
    end

    if loading >= cong_threshold
        lineColor = 'r';
        lineWidth = 3;
    else
        lineColor = [0.4 0.4 0.4];
        lineWidth = 1.5;
    end

    plot([x1, x2], [y1, y2], '-', ...
        'Color', lineColor, 'LineWidth', lineWidth);

    xm = 0.5*(x1 + x2);
    ym = 0.5*(y1 + y2);

    text(xm, ym, sprintf('%.1f MW', flow), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 7);

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

scatter(bus_x, bus_y, 80, LMP, 'filled');
colormap(jet);
cb = colorbar;
ylabel(cb, 'LMP [$/MWh]');

cmin = min(LMP);
cmax = max(LMP);
if abs(cmax - cmin) < 1e-3
    caxis([cmin - 0.01, cmax + 0.01]);
else
    caxis([cmin, cmax]);
end

for k = 1:nb
    text(bus_x(k) + 0.01, bus_y(k) + 0.01, ...
        sprintf('%d\n%.2f', k, LMP(k)), ...
        'FontSize', 7, 'FontWeight', 'bold');
end

h1 = plot(nan, nan, 'r-',  'LineWidth', 3);
h2 = plot(nan, nan, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);
legend([h1, h2], {'Congested line', 'Non-congested line'}, ...
       'Location', 'bestoutside');

title('Adjusted IEEE 39-bus DC-OPF: LMPs, Power Flows, and Congestion (graph layout)');
axis equal;
axis off;
hold off;

