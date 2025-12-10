define_constants;

bussys  = loadcase('case39');
bus     = bussys.bus;
branch  = bussys.branch;
gen     = bussys.gen;
gencost = bussys.gencost;

nb = size(bus,    1);
nl = size(branch, 1);
ng = size(gen,    1);

%% -------------------------------
% Solar (capacity fraction series)
%% -------------------------------
solardata = readtable("ISONE_solarcapacity_modeled_24yr_EPT.csv");

solarcap = solardata.Capacity(1);           % MW total capacity
solargen = solardata.ISONE_pv_pwr;          % capacity fraction (0–1)

% Align range (consistent with your earlier trimming intent)
solargen = solargen(6:end-1);

% Per-unit optional
solargen_pu = solargen * solarcap / bussys.baseMVA;

%% -------------------------------
% Wind (capacity fraction series)
%% -------------------------------
windgentot = readtable("2024_isone_all_windplant_timeseries_locations.csv");
windcap = sum(windgentot.NameplateCapacity_MW_);   % MW total capacity

winddata = readtable("2024_isone_wind_ofsw_aggregated_power_data_2000_2023_attk_eleclosspwrts.csv");
windgen = winddata(:,9:end);
windgen = table2array(windgen);
windgen = reshape(windgen, [], 1);                 % capacity fraction (0–1)

% Per-unit optional
windgen_pu = windgen * windcap / bussys.baseMVA;

%% -------------------------------
% Load (optional, if you still want it)
%% -------------------------------
loaddata = readtable("ISONE_grossload_metdata_spliced_24yr_UTC.csv");

addload = loaddata.ISONE_grs_ld / bussys.baseMVA;
addload = addload(1:end-6);

%% -------------------------------
% Highest 30-day mean CF window
%% -------------------------------
window_hours = 30 * 24;

Ns = length(solargen);
bestSolarVal = -inf;
bestSolarIdx = 1;

for k = 1:(Ns - window_hours + 1)
    sVal = mean(solargen(k:k+window_hours-1));
    if sVal > bestSolarVal
        bestSolarVal = sVal;
        bestSolarIdx = k;
    end
end

Nw = length(windgen);
bestWindVal = -inf;
bestWindIdx = 1;

for k = 1:(Nw - window_hours + 1)
    wVal = mean(windgen(k:k+window_hours-1));
    if wVal > bestWindVal
        bestWindVal = wVal;
        bestWindIdx = k;
    end
end

fprintf('Best 30-day SOLAR window starts at index %d (mean CF %.4f)\n', ...
    bestSolarIdx, bestSolarVal);
fprintf('Best 30-day WIND window starts at index %d (mean CF %.4f)\n', ...
    bestWindIdx, bestWindVal);

%% -------------------------------
% 24-hour window + peak hour inside it
%% -------------------------------
if exist('bestSolarIdx','var') && exist('bestWindIdx','var')
    startIdx24 = min(bestSolarIdx, bestWindIdx);
elseif exist('bestSolarIdx','var')
    startIdx24 = bestSolarIdx;
elseif exist('bestWindIdx','var')
    startIdx24 = bestWindIdx;
else
    startIdx24 = 1;
end

T = 24;
Nseries = min(length(solargen), length(windgen));

if startIdx24 + T - 1 > Nseries
    startIdx24 = max(1, Nseries - T + 1);
end

hoursIdx24 = startIdx24:(startIdx24+T-1);

[~, iS] = max(solargen(hoursIdx24));
[~, iW] = max(windgen(hoursIdx24));

peakSolarIdx24 = hoursIdx24(iS);
peakWindIdx24  = hoursIdx24(iW);

peakSolarCF24 = solargen(peakSolarIdx24);
peakWindCF24  = windgen(peakWindIdx24);

fprintf('24-hour window start index: %d\n', startIdx24);
fprintf('Peak SOLAR hour in 24h window: idx %d (CF %.4f)\n', ...
    peakSolarIdx24, peakSolarCF24);
fprintf('Peak WIND hour in 24h window: idx %d (CF %.4f)\n', ...
    peakWindIdx24, peakWindCF24);


