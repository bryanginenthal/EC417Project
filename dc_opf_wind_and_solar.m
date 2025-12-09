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
