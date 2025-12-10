clear
clc
formatdata;

Pd0 = bus(:, PD);

gen_bus = gen(:, GEN_BUS);
Pg_min  = gen(:, PMIN);
Pg_max  = gen(:, PMAX);

Fmax = branch(:, RATE_A);
Fmax(Fmax == 0) = 1e4;

PTDF = makePTDF(bussys);

Cg = zeros(nb, ng);
for i = 1:ng
    Cg(gen_bus(i), i) = 1;
end

%% -------------------------------------------------
% (added step 1) Define 24-hour window start from 30-day peaks
%% -------------------------------------------------
if exist('bestSolarIdx','var') && exist('bestWindIdx','var')
    startIdx24 = min(bestSolarIdx, bestWindIdx);        % (added step 1)
elseif exist('bestSolarIdx','var')
    startIdx24 = bestSolarIdx;                          % (added step 1)
elseif exist('bestWindIdx','var')
    startIdx24 = bestWindIdx;                           % (added step 1)
else
    startIdx24 = 1;                                     % (added step 1)
end

T = 24;                                                  % (added step 1)
Nseries = min(length(solargen), length(windgen));        % (added step 1)
if startIdx24 + T - 1 > Nseries
    startIdx24 = max(1, Nseries - T + 1);               % (added step 1)
end
hoursIdx24 = startIdx24:(startIdx24+T-1);               % (added step 1)

%% =================================================
% WIND block
%% =================================================

% (added step 1) Peak hour within the 24-hour window
[~, iW] = max(windgen(hoursIdx24));                     % (added step 1)
peakWindIdx24 = hoursIdx24(iW);                         % (added step 1)
peakWindCF24  = windgen(peakWindIdx24);                 % (added step 1)

Pw_MW_peak = peakWindCF24 * windcap;                    % (added step 1)

pg    = sdpvar(ng,1);
pwind = sdpvar(nb,1);
zwind = binvar(nb,1);

ConsW = [];
ConsW = [ConsW; sum(Cg*pg) + sum(pwind) == sum(Pd0)];
fW = PTDF*(Cg*pg + pwind - Pd0);
ConsW = [ConsW; -Fmax <= fW <= Fmax];
ConsW = [ConsW; Pg_min <= pg <= Pg_max];

ConsW = [ConsW; 0 <= pwind <= Pw_MW_peak .* zwind];     % (added step 1)
ConsW = [ConsW; sum(zwind) == 1];

objW = sum(pg);
optsW = sdpsettings('solver','intlinprog','verbose',0);
optimize(ConsW, objW, optsW);

zwind_opt = value(zwind);
windbus = find(zwind_opt > 0.5);

fprintf('(added step 1) WIND peak hour idx: %d\n', peakWindIdx24);
fprintf('Wind placed at bus %d (CF %.4f, peak %.2f MW).\n', ...
    windbus, peakWindCF24, Pw_MW_peak);

%% =================================================
% SOLAR block
%% =================================================

% (added step 1) Peak hour within the 24-hour window
[~, iS] = max(solargen(hoursIdx24));                    % (added step 1)
peakSolarIdx24 = hoursIdx24(iS);                        % (added step 1)
peakSolarCF24  = solargen(peakSolarIdx24);              % (added step 1)

Ps_MW_peak = peakSolarCF24 * solarcap;                  % (added step 1)

pg     = sdpvar(ng,1);
psolar = sdpvar(nb,1);
zsolar = binvar(nb,1);

ConsS = [];
ConsS = [ConsS; sum(Cg*pg) + sum(psolar) == sum(Pd0)];
fS = PTDF*(Cg*pg + psolar - Pd0);
ConsS = [ConsS; -Fmax <= fS <= Fmax];
ConsS = [ConsS; Pg_min <= pg <= Pg_max];

ConsS = [ConsS; 0 <= psolar <= Ps_MW_peak .* zsolar];   % (added step 1)
ConsS = [ConsS; sum(zsolar) == 1];

objS = sum(pg);
optsS = sdpsettings('solver','intlinprog','verbose',0);
optimize(ConsS, objS, optsS);

zsolar_opt = value(zsolar);
solarbus = find(zsolar_opt > 0.5);

fprintf('(added step 1) SOLAR peak hour idx: %d\n', peakSolarIdx24);
fprintf('Solar placed at bus %d (CF %.4f, peak %.2f MW).\n', ...
    solarbus, peakSolarCF24, Ps_MW_peak);

%% -------------------------------------------------
% (added step 1) Save outputs for summary script
%% -------------------------------------------------
save('peak_placement_result.mat', ...
    'windbus','solarbus', ...
    'startIdx24','hoursIdx24', ...
    'peakWindIdx24','peakSolarIdx24', ...
    'peakWindCF24','peakSolarCF24', ...
    'windcap','solarcap');

