% build_renewables_package.m
% Assemble 24-hour renewable injections for case39.

clc; clear; close all;

% Load data (defines bussys, windgen, solargen, hoursIdx24, etc.)
formatdata;

Nbus = nb;

% Bus selections (edit as needed)
windbus  = 30;
solarbus = 31;

% Fallback if not set in formatdata
if ~exist('hoursIdx24','var')
    hoursIdx24 = 1:24;
end

T = length(hoursIdx24);

% 24-hour MW profiles
Pwind_MW  = windgen(hoursIdx24)  * windcap;
Psolar_MW = solargen(hoursIdx24) * solarcap;

% Bus-level matrices
Pwind_bus_MW  = zeros(Nbus, T);
Psolar_bus_MW = zeros(Nbus, T);

% Wind placement
if numel(windbus) == 1
    Pwind_bus_MW(windbus, :) = Pwind_MW.';
else
    share = Pwind_MW.' / numel(windbus);
    for k = 1:numel(windbus)
        Pwind_bus_MW(windbus(k), :) = Pwind_bus_MW(windbus(k), :) + share;
    end
end

% Solar placement
if numel(solarbus) == 1
    Psolar_bus_MW(solarbus, :) = Psolar_MW.';
else
    share = Psolar_MW.' / numel(solarbus);
    for k = 1:numel(solarbus)
        Psolar_bus_MW(solarbus(k), :) = Psolar_bus_MW(solarbus(k), :) + share;
    end
end

% Save package
save('ren_profiles_case39.mat', ...
    'windbus','solarbus', ...
    'hoursIdx24', ...
    'Pwind_MW','Psolar_MW', ...
    'Pwind_bus_MW','Psolar_bus_MW', ...
    'windcap','solarcap');

disp('Saved: ren_profiles_case39.mat');

% Table for quick viewing
Hour    = (1:T).';
WindMW  = Pwind_MW(:);
SolarMW = Psolar_MW(:);
disp(table(Hour, WindMW, SolarMW));

% Plot profiles
figure;
plot(1:T, Pwind_MW,  '-o', 'MarkerFaceColor','auto'); hold on;
plot(1:T, Psolar_MW, '-o', 'MarkerFaceColor','auto');
xlabel('Hour');
ylabel('MW');
title('24-Hour Wind and Solar Profiles');
legend('Wind','Solar','Location','best');
grid on;
