% build_load_timeseries.m
%
% Convert ISO-NE UTC gross load CSV files into a single timetable
% with one column per region (including ISONE) for use in the 39-bus project.
%
% Assumes this file lives in the 417_project folder and the CSVs are in:
%   417_project/load_data/Region_grossload_metdata_spliced_24yr_UTC.csv

function build_load_timeseries()
    % Get project root as folder containing this script
    projectRoot = fileparts(mfilename('fullpath'));
    dataFolder  = fullfile(projectRoot, 'load_data');
    
    % List of regions corresponding to your files
    regions = {'CT','ME','NH','RI','VT','NEMA','WCMA','SEMA','ISONE'};
    
    % Initialize combined timetable
    loadTT = timetable();
    
    for rIdx = 1:numel(regions)
        region  = regions{rIdx};
        fileName = sprintf('%s_grossload_metdata_spliced_24yr_UTC.csv', region);
        fullPath = fullfile(dataFolder, fileName);
        
        if ~isfile(fullPath)
            error('File not found: %s', fullPath);
        end
        
        fprintf('Reading %s ...\n', fullPath);
        
        % Read the CSV
        opts = detectImportOptions(fullPath);
        tbl  = readtable(fullPath, opts);
        
        % --- Build UTC datetime vector ---
        % Assumes columns named 'Date' and 'Hour_Ending'
        
        % Let MATLAB parse the date string (ISO-NE uses standard date formats)
        dt = datetime(tbl.Date, 'TimeZone', 'UTC');
        
        % Hour_Ending might be numeric (1..24) or strings (e.g. 'HE01')
        he = tbl.Hour_Ending;
        if iscell(he) || isstring(he)
            % Strip non-digits and convert to number
            heNum = str2double(regexprep(string(he), '\D', ''));
        else
            heNum = he;
        end
        
        % Define timestamp at the *end* of the hour.
        % Example: Date = 2000-01-01, Hour_Ending = 1 -> 2000-01-01 01:00 UTC
        dt = dt + hours(heNum);
        
        % --- Extract load column for this region ---
        loadVarName = sprintf('%s_grs_ld', region);
        if ~ismember(loadVarName, tbl.Properties.VariableNames)
            error('Could not find variable %s in %s', loadVarName, fileName);
        end
        
        loadMW = tbl.(loadVarName);  % gross load in MW
        
        % Create a timetable with one column named after the region
        ttRegion = timetable(dt, loadMW, 'VariableNames', {region});
        
        % Combine with global timetable (intersection ensures common timestamps)
        if rIdx == 1
            loadTT = ttRegion;
        else
            loadTT = synchronize(loadTT, ttRegion, 'intersection');
        end
    end
    
    % Save the result for later tasks
    outFile = fullfile(projectRoot, 'load_timeseries_UTC.mat');
    save(outFile, 'loadTT', 'regions');
    
    fprintf('Saved combined load timetable to %s\n', outFile);
end


