function build_solar_timeseries()
    % BUILD_SOLAR_TIMESERIES
    % ----------------------
    % Reads all regional solar CSVs and builds a single timetable "solarTT"
    % with:
    %   - Row times: dt (UTC)
    %   - Columns: all *_pv_pwr series (regional + plant-level)
    %
    % For EACH CSV:
    %   Col 1: Date            (format 'M/d/yyyy', e.g. '1/1/2000')
    %   Col 2: Hour_Ending     (1..24 or string like "HE01")
    %   Col 3+: solar MW series
    %
    % We DO NOT rely on the file's header names (Var1, Var2, etc.).
    % Instead, we hard-code the intended variable names by region.
    %
    % File pattern:
    %   Region_solarcapacity_modeled_24yr_UTC.csv
    % where Region = CT, ME, NH, RI, VT, NEMA, WCMA, SEMA, ISONE.

    % Root of the project (where this script lives)
    projectRoot = fileparts(mfilename('fullpath'));
    dataFolder  = fullfile(projectRoot, 'solar_data');

    % Regions and their intended solar variable names (in order of columns 3:end)
    regionList = { ...
        'CT',   {'CT_pv_pwr',   'HartfordSolar_pv_pwr'}; ...
        'ME',   {'ME_pv_pwr',   'HancockSolar_pv_pwr', 'YorkSolar_pv_pwr'}; ...
        'NH',   {'NH_pv_pwr',   'CarrollSolar_pv_pwr'}; ...
        'RI',   {'RI_pv_pwr',   'CranstonSolar_pv_pwr'}; ...
        'VT',   {'VT_pv_pwr',   'AddisonSolar_pv_pwr', 'CoolidgeSolar_pv_pwr'}; ...
        'NEMA', {'NEMA_pv_pwr'}; ...
        'WCMA', {'WCMA_pv_pwr', 'SpencerSolar_pv_pwr'}; ...
        'SEMA', {'SEMA_pv_pwr'}; ...
        'ISONE',{'ISONE_pv_pwr'} ...
    };

    solarTT = timetable();  % combined result

    for rIdx = 1:size(regionList, 1)
        region   = regionList{rIdx, 1};
        pvVars   = regionList{rIdx, 2};   % cell array of names for col3+

        fileName = sprintf('%s_solarcapacity_modeled_24yr_UTC.csv', region);
        fullPath = fullfile(dataFolder, fileName);

        if ~isfile(fullPath)
            error('File not found: %s', fullPath);
        end

        fprintf('Reading %s ...\n', fullPath);

        % Read the CSV. Headers might be "Var1", "Var2", etc; we're not using them.
        tbl = readtable(fullPath, 'ReadVariableNames', true, ...
                                   'TextType', 'string', ...
                                   'VariableNamingRule', 'preserve');

        % --- Build datetime vector (UTC) ---
        % Col 1: Date (e.g. '1/1/2000')
        dateCol = tbl{:, 1};

        if isstring(dateCol) || iscellstr(dateCol)
            dt = datetime(dateCol, ...
                'InputFormat', 'M/d/yyyy', ...
                'TimeZone', 'UTC');
        elseif isdatetime(dateCol)
            dt = dateCol;
            dt.TimeZone = 'UTC';
        else
            error('Unexpected date column type in %s', fileName);
        end

        % Col 2: Hour_Ending (1..24 or strings like 'HE01')
        hourCol = tbl{:, 2};
        if isstring(hourCol) || iscell(hourCol)
            heNum = str2double(regexprep(string(hourCol), '\D', ''));
        else
            heNum = hourCol;
        end

        dt = dt + hours(heNum);

        % --- Extract solar data (columns 3:end) ---
        nCols    = width(tbl);
        nSolar   = nCols - 2;
        nNames   = numel(pvVars);

        if nSolar ~= nNames
            error('In %s: expected %d solar columns (3:end), but found %d.', ...
                   fileName, nNames, nSolar);
        end

        solarData = tbl{:, 3:end};   % numeric MW matrix (N x nSolar)

        % Build a small table and rename its variables to our intended names
        subTbl = array2table(solarData, 'VariableNames', pvVars);

        % Attach dt and convert to timetable
        subTbl.dt = dt;
        ttRegion = table2timetable(subTbl, 'RowTimes', 'dt');

        % Merge with global timetable using common timestamps
        if rIdx == 1
            solarTT = ttRegion;
        else
            solarTT = synchronize(solarTT, ttRegion, 'intersection');
        end
    end

    % Save combined timetable
    outFile = fullfile(projectRoot, 'solar_timeseries_UTC.mat');
    save(outFile, 'solarTT');

    fprintf('Saved combined solar timetable to %s\n', outFile);
    fprintf('Total variables in solarTT: %d\n', width(solarTT));
end
