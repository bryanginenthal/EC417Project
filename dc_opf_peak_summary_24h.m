clear
clc
formatdata;

load('peak_placement_result.mat');

%% -------------------------------------------------
% (added step 1) 24-hour profiles
%% -------------------------------------------------
wind24_cf  = windgen(hoursIdx24);                      % (added step 1)
solar24_cf = solargen(hoursIdx24);                     % (added step 1)

wind24_MW  = wind24_cf  * windcap;                     % (added step 1)
solar24_MW = solar24_cf * solarcap;                    % (added step 1)

%% -------------------------------------------------
% (added step 1) Compact summary table
%% -------------------------------------------------
Summary = table( ...
    ["Wind"; "Solar"], ...
    [peakWindIdx24; peakSolarIdx24], ...
    [windbus; solarbus], ...
    [peakWindCF24; peakSolarCF24], ...
    [peakWindCF24*windcap; peakSolarCF24*solarcap], ...
    'VariableNames', ...
    {'Resource','PeakHourIndex','SelectedBus','PeakCF','PeakMW'});

disp(Summary);

%% -------------------------------------------------
% (added step 1) 24-hour CF plot
%% -------------------------------------------------
figure;
plot(1:24, wind24_cf, '-o'); hold on;
plot(1:24, solar24_cf, '-o');
xlabel('Hour in 24-hour window');
ylabel('Capacity fraction (0-1)');
legend('Wind CF','Solar CF');
title('(added step 1) 24-hour renewable capacity profiles');
grid on;

%% -------------------------------------------------
% (added step 1) 24-hour MW plot
%% -------------------------------------------------
figure;
plot(1:24, wind24_MW, '-o'); hold on;
plot(1:24, solar24_MW, '-o');
xlabel('Hour in 24-hour window');
ylabel('MW');
legend('Wind MW','Solar MW');
title('(added step 1) 24-hour renewable MW profiles');
grid on;

