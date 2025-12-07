build_solar_timeseries

load solar_timeseries_UTC.mat

% Peek at first rows
solarTT(1:5, :)

% See all variable names
solarTT.Properties.VariableNames'

y2020 = timerange('2020-01-01','2020-12-31');
plot(solarTT.dt(y2020), solarTT.ISONE_pv_pwr(y2020));
xlabel('Time (UTC)');
ylabel('ISO-NE Solar [MW]');
title('ISO-NE Total Solar (ISONE\_pv\_pwr) - 2020');