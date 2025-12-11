% formatdata.m
define_constants;
bussys = loadcase('case39');   % MATPOWER named indices (BUS_I, PD, etc.)
bus    = bussys.bus;
branch = bussys.branch;
gen    = bussys.gen;
gencost= bussys.gencost;

nb = size(bus,    1);      % number of buses
nl = size(branch, 1);      % number of lines
ng = size(gen,    1);      % number of generators

%% -------------------------------
% Solar (capacity fraction series)
%% -------------------------------
solardata = readtable("ISONE_solarcapacity_modeled_24yr_EPT.csv");
solarcap  = solardata.Capacity(1);
solargen  = solardata.ISONE_pv_pwr * solarcap / bussys.baseMVA; % per-unit on system base
solargen  = solargen(6:end-1); % align to 1/1/2000–12/31/2023

% restrict to 672-hour window
solargen  = solargen(1:672);

% keep only non-zero hours, track indices
datatimes = find(solargen ~= 0);
solargen  = solargen(datatimes);

%% -------------------------------
% Wind (aggregated)
%% -------------------------------
windgentot = readtable("2024_isone_all_windplant_timeseries_locations.csv");
windcap    = sum(windgentot.NameplateCapacity_MW_); % total MW capacity

winddata = readtable("2024_isone_wind_ofsw_aggregated_power_data_2000_2023_attk_eleclosspwrts.csv");
windgen  = winddata(:,9:end);
windgen  = table2array(windgen);
windgen  = reshape(windgen, [], 1);
% windgen  = windgen * windcap / bussys.baseMVA;

windgen  = windgen(1:672);
windgen  = windgen(datatimes);  % aligned with solargen

%% -------------------------------
% Load (kept as in original)
%% -------------------------------
loaddata = readtable("ISONE_grossload_metdata_spliced_24yr_UTC.csv");
addload  = loaddata.ISONE_grs_ld / bussys.baseMVA;
addload  = addload(1:end-6);
addload  = addload(1:672);

%% -------------------------------
% Peak-hour indices (within current window)
%% -------------------------------
[~, idx_peak_solar] = max(solargen);  % index in aligned solargen/windgen/addload
[~, idx_peak_wind]  = max(windgen);
