define_constants;
bussys = loadcase('case39');   % MATPOWER named indices (BUS_I, PD, etc.)
bus    = bussys.bus;
branch = bussys.branch;
gen    = bussys.gen;
gencost= bussys.gencost;

nb = size(bus,    1);      % number of buses
nl = size(branch, 1);      % number of lines
ng = size(gen,    1);      % number of generators

solardata = readtable("2024_isone_variable_energy_resource_ver_data_series_2000-2023_rev1\solar\ISONE_solarcapacity_modeled_24yr_EPT.csv");
solarcap = solardata.Var4(1);
solargen = solardata.Var3 * solarcap / bussys.baseMVA; % Format generation to per-unit norm of system
solargen = solargen(6:end-1); % Capture same amount of data as wind data (starting from 1/1/2000, ending 12/31/2023)

solargen = solargen(1:100);
datatimes = find(solargen~=0);
solargen = solargen(datatimes);


windgentot = readtable("2024_isone_variable_energy_resource_ver_data_series_2000-2023_rev1\wind\2024_isone_all_windplant_timeseries_locations.csv");
windcap = sum(windgentot.NameplateCapacity_MW_); % Capacity as sum of nameplate capacities of wind farms in ISONE area

winddata = readtable("2024_isone_variable_energy_resource_ver_data_series_2000-2023_rev1\wind\2024_isone_wind_ofsw_aggregated_power_data_2000_2023_attk_eleclosspwrts.csv");
windgen = winddata(:,9:end);
windgen = table2array(windgen);
windgen = reshape(windgen, [], 1);
%windgen = windgen * windcap / bussys.baseMVA;

windgen = windgen(1:100);
windgen = windgen(datatimes);

loaddata = readtable("2024_isone_variable_energy_resource_ver_data_series_2000-2023_rev1\load\ISONE_grossload_metdata_spliced_24yr_UTC.csv");
addload = loaddata.ISONE_grs_ld / bussys.baseMVA;
addload = addload(1:end-6);
addload = addload(1:24);