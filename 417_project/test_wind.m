build_wind_timeseries

load wind_timeseries_UTC.mat

% See first few rows
head(windTT)
% Rename row times to a common name 'Time'
windTT.Properties.DimensionNames{1} = 'Time';
%% Try to standardize all dtOn / dtOff variables to time later...

% Plot one year as a sanity check
y2020 = timerange('2020-01-01','2020-12-31');
plot(windTT.Time(y2020), windTT.TotalWind(y2020));
xlabel('Time (UTC)');
ylabel('Total wind (MW or p.u.)');
title('ISO-NE Total Wind – 2020');
