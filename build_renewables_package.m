clc; clear; close all;

formatdata;
load('peak_placement_result.mat');   % windbus, solarbus, startIdx24, hoursIdx24

% Number of buses
Nbus = nb;

% Use a clean window (24-hour window already defined)
if ~exist('hoursIdx24','var')
    hoursIdx24 = 1:24;
end

T = length(hoursIdx24);

% Renewable MW profiles over the window
Pwind_MW  = windgen(hoursIdx24)  * windcap;
Psolar_MW = solargen(hoursIdx24) * solarcap;

% Bus injection matrices
Pwind_bus_MW  = zeros(Nbus, T);
Psolar_bus_MW = zeros(Nbus, T);

% Map to selected buses
Pwind_bus_MW(windbus, :)   = Pwind_MW.';
Psolar_bus_MW(solarbus, :) = Psolar_MW.';

% Save mapping package
save('ren_profiles_case39.mat', ...
    'windbus','solarbus', ...
    'hoursIdx24', ...
    'Pwind_MW','Psolar_MW', ...
    'Pwind_bus_MW','Psolar_bus_MW', ...
    'windcap','solarcap');

% Display results
disp('Task 4 Package Saved: ren_profiles_case39.mat');
disp('Includes bus mapping and 24-hour aligned injection matrices.');

%% Maybe delete, or use for report

% table for the 24h window
T = length(hoursIdx24);
Hour = (1:T).';
WindMW  = Pwind_MW(:);
SolarMW = Psolar_MW(:);

ProfileTable = table(Hour, WindMW, SolarMW);
disp(ProfileTable)

% bus injection matrices
disp('Wind MW (24h):');  disp(Pwind_MW.')
disp('Solar MW (24h):'); disp(Psolar_MW.')

load('ren_profiles_case39.mat');

figure;
plot(1:length(Pwind_MW), Pwind_MW, '-o', 'MarkerFaceColor','auto'); hold on;
plot(1:length(Psolar_MW), Psolar_MW, '-o', 'MarkerFaceColor','auto');
xlabel('Hour In 24-Hour Window');
ylabel('MW');
title('Task 4: 24-Hour Wind And Solar Injections');
legend('Wind','Solar');
grid on;


