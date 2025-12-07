load('load_timeseries_UTC.mat', 'loadTT')

% Plot system-level load for a sample year:
y2020 = timerange('2020-01-01','2020-12-31 23:59:59');

plot(loadTT.dt(y2020), loadTT.ISONE(y2020));
ylabel('ISONE gross load [MW]');
xlabel('Time (UTC)');
title('ISO-NE Gross Load — 2020');

load('load_timeseries_UTC.mat');
whos

loadTT.Properties.DimensionNames
loadTT.Properties.RowTimes

head(loadTT)

summary(loadTT)


