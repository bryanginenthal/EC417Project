clc; clear; close all;

define_constants;

% Load the IEEE-39 case
mpc = loadcase('case39');

bus    = mpc.bus;
branch = mpc.branch;
gen    = mpc.gen;

% Run a baseline AC power flow
results = runpf(mpc);

if results.success
    disp('Baseline AC PF converged.');
else
    disp('AC PF did not converge.');
end

% Basic counts
nb = size(bus,1);
nl = size(branch,1);
ng = size(gen,1);

fprintf('Buses: %d\n', nb);
fprintf('Branches: %d\n', nl);
fprintf('Generators: %d\n', ng);

% Bus type breakdown
types = bus(:, BUS_TYPE);
fprintf('Slack: %d\n', sum(types == 3));
fprintf('PV: %d\n',    sum(types == 2));
fprintf('PQ: %d\n',    sum(types == 1));

% Quick look at generator limits
GenTable = table(gen(:,GEN_BUS), gen(:,PMIN), gen(:,PMAX), ...
    'VariableNames', {'GenBus','Pmin','Pmax'});

disp('Generator limits (first 10):');
disp(GenTable(1:min(10,height(GenTable)),:));

% Branch rating check
RateA = branch(:, RATE_A);
fprintf('Branches with RATE_A = 0: %d\n', sum(RateA == 0));

% Save basic info
baseline.success = results.success;
baseline.nb = nb;
baseline.nl = nl;
baseline.ng = ng;
baseline.bus_type_counts = [sum(types==3), sum(types==2), sum(types==1)];
baseline.GenTable = GenTable;
baseline.RateA_zero = sum(RateA == 0);

save('task2_baseline_report.mat','baseline');
disp('Saved task2_baseline_report.mat');

