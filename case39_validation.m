clc; clear; close all;

define_constants;

% Load IEEE-39
mpc = loadcase('case39');

bus    = mpc.bus;
branch = mpc.branch;
gen    = mpc.gen;

% Baseline AC power flow
results = runpf(mpc);

% Convergence check
if results.success
    disp('Baseline AC Power Flow Converged.');
else
    disp('Baseline AC Power Flow Did Not Converge.');
end

% Basic system summary
nb = size(bus,1);
nl = size(branch,1);
ng = size(gen,1);

fprintf('Buses: %d\n', nb);
fprintf('Branches: %d\n', nl);
fprintf('Generators: %d\n', ng);

% Bus type counts
bus_types = bus(:, BUS_TYPE);
fprintf('Slack buses: %d\n', sum(bus_types == 3));
fprintf('PV buses: %d\n', sum(bus_types == 2));
fprintf('PQ buses: %d\n', sum(bus_types == 1));

% Generator limits check (simple snapshot table)
GenTable = table( ...
    gen(:, GEN_BUS), gen(:, PMIN), gen(:, PMAX), ...
    'VariableNames', {'GenBus','Pmin','Pmax'});

disp('Generator Capability Limits (first 10):');
disp(GenTable(1:min(10,height(GenTable)),:));

% Branch ratings check (RATE_A)
RateA = branch(:, RATE_A);
fprintf('Branches with RATE_A = 0: %d\n', sum(RateA == 0));

% Save a small report struct
baseline.success = results.success;
baseline.nb = nb; baseline.nl = nl; baseline.ng = ng;
baseline.bus_type_counts = [sum(bus_types==3), sum(bus_types==2), sum(bus_types==1)];
baseline.GenTable = GenTable;
baseline.RateA_zero = sum(RateA == 0);

save('task2_baseline_report.mat','baseline');

disp('Saved: task2_baseline_report.mat');
