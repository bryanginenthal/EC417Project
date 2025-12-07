function build_wind_timeseries()
% BUILD_WIND_TIMESERIES
% Reads ISO-NE onshore and offshore NET wind power data (2000–2023),
% each row = 1 day with HE1..HE24, and builds an hourly timetable:
%   windTT with variables:
%       OnshoreWind
%       OffshoreWind
%       TotalWind

    % Project root = folder containing this script
    projectRoot = fileparts(mfilename('fullpath'));
    dataFolder  = fullfile(projectRoot, 'wind_data');

    % File names (adjust if your names differ)
    onFile  = fullfile(dataFolder, ...
        '2024_isone_wind_onsw_aggregated_power_data_2000_2023_netpwrts.csv');
    offFile = fullfile(dataFolder, ...
        '2024_isone_wind_ofsw_aggregated_power_data_2000_2023_netpwrts.csv');

    if ~isfile(onFile)
        error('Onshore wind file not found: %s', onFile);
    end
    if ~isfile(offFile)
        error('Offshore wind file not found: %s', offFile);
    end

    fprintf('Reading onshore wind data...\n');
    tblOn  = readtable(onFile);    % default options are fine

    fprintf('Reading offshore wind data...\n');
    tblOff = readtable(offFile);

    % ----- Build hourly timestamps for onshore -----
    % Column 2 = Date (e.g. '1/1/2000')
    dateOn = datetime(tblOn{:, 2}, 'InputFormat', 'M/d/yyyy', 'TimeZone', 'UTC');

    % Columns 9..32 = HE1..HE24
    heCols = 9:32;
    onDaily = tblOn{:, heCols};    % size: [nDays x 24]

    nDaysOn = size(onDaily, 1);
    nHoursOn = nDaysOn * 24;

    dtOn   = NaT(nHoursOn, 1, 'TimeZone', 'UTC');
    onVals = zeros(nHoursOn, 1);

    idx = 1;
    for d = 1:nDaysOn
        for h = 0:23
            dtOn(idx)   = dateOn(d) + hours(h);
            onVals(idx) = onDaily(d, h+1);  % h+1 because HE1 is column 1 in onDaily
            idx = idx + 1;
        end
    end

    ttOn = timetable(dtOn, onVals, 'VariableNames', {'OnshoreWind'});

    % ----- Build hourly timestamps for offshore -----
    dateOff = datetime(tblOff{:, 2}, 'InputFormat', 'M/d/yyyy', 'TimeZone', 'UTC');
    offDaily = tblOff{:, heCols};

    nDaysOff = size(offDaily, 1);
    nHoursOff = nDaysOff * 24;

    dtOff    = NaT(nHoursOff, 1, 'TimeZone', 'UTC');
    offVals  = zeros(nHoursOff, 1);

    idx = 1;
    for d = 1:nDaysOff
        for h = 0:23
            dtOff(idx)   = dateOff(d) + hours(h);
            offVals(idx) = offDaily(d, h+1);
            idx = idx + 1;
        end
    end

    ttOff = timetable(dtOff, offVals, 'VariableNames', {'OffshoreWind'});

    % ----- Synchronize and add total -----
    windTT = synchronize(ttOn, ttOff);   % align by time
    windTT.TotalWind = windTT.OnshoreWind + windTT.OffshoreWind;

    % ----- Save result -----
    outFile = fullfile(projectRoot, 'wind_timeseries_UTC.mat');
    save(outFile, 'windTT');

    fprintf('Saved wind_timeseries_UTC.mat\n');
    disp(windTT.Properties.VariableNames');
end
