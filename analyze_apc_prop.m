function outFiles = analyze_apc_prop(dataFolder, propSpec)
%ANALYZE_APC_PROP Fit APC prop C_T(J,Re) / C_P(J,Re) and write JSON.
%
%   analyze_apc_prop()
%       Prompts in the command window for folder + prop selection.
%
%   analyze_apc_prop(dataFolder)
%       If the folder has exactly one propeller (performance+geometry pair),
%       analyzes it. If multiple are present, lists them and asks you to
%       re-run with an explicit propSpec.
%
%   analyze_apc_prop(dataFolder, propSpec)
%       propSpec examples:
%           '5x11E'
%           '5x11E,7x11E'
%           {'5x11E','7x15E'}
%           'all'
%
%   outFiles = analyze_apc_prop(...)
%       Returns cell array of written JSON paths.
%
% Expected files in dataFolder (same naming scheme):
%   <name>_performance.dat
%   <name>_geometry.dat
%
% Output:
%   <dataFolder>/<name>_analysed.json

    if nargin < 1 || isempty(dataFolder)
        dataFolder = input('Path to folder with APC .dat files: ', 's');
    end
    if nargin < 2
        propSpec = [];
    end

    dataFolder = char(dataFolder);
    assert(isfolder(dataFolder), 'Data folder not found: %s', dataFolder);

    available = discoverProps(dataFolder);
    assert(~isempty(available), ...
        'No *_performance.dat files found in %s', dataFolder);

    propsToRun = resolvePropSelection(available, propSpec);
    if isempty(propsToRun)
        outFiles = {};
        return;
    end

    % Bivariate polynomial degrees (editable)
    degreeJ  = 4;      % powers of advance ratio J
    degreeRe = 2;      % powers of scaled Reynolds number
    ReScale  = 1e5;    % Re_hat = Re / ReScale  (keeps matrix well-conditioned)

    outFiles = cell(numel(propsToRun), 1);
    fprintf('Analyzing %d propeller(s) in %s\n', numel(propsToRun), dataFolder);
    fprintf('Fit: bivariate poly, deg_J=%d, deg_Re=%d, Re_hat = Re / %.0e\n', ...
        degreeJ, degreeRe, ReScale);

    for k = 1:numel(propsToRun)
        propName = propsToRun{k};
        fprintf('\n---- %s ----\n', propName);
        result = analyzeOneProp(dataFolder, propName, degreeJ, degreeRe, ReScale);
        outPath = fullfile(dataFolder, [propName '_analysed.json']);
        writeJson(outPath, result);
        outFiles{k} = outPath;
        fprintf('Wrote %s\n', outPath);
        fprintf('  CT RMSE = %.4g | CP RMSE = %.4g | N = %d pts\n', ...
            result.fit_model.rmse_CT, result.fit_model.rmse_CP, ...
            result.fit_model.n_points_used);
        fprintf('  J in [%.3f, %.3f] | Re in [%.0f, %.0f]\n', ...
            result.validity.J_min, result.validity.J_max, ...
            result.validity.Re_min, result.validity.Re_max);
    end
end

%% -------- selection / discovery --------

function names = discoverProps(dataFolder)
    files = dir(fullfile(dataFolder, '*_performance.dat'));
    names = cell(numel(files), 1);
    for i = 1:numel(files)
        names{i} = erase(files(i).name, '_performance.dat');
    end
    names = unique(names, 'stable');
end

function propsToRun = resolvePropSelection(available, propSpec)
    propsToRun = {};

    if isempty(propSpec)
        if numel(available) == 1
            propsToRun = available;
            fprintf('One propeller found (%s). Running analysis.\n', available{1});
            return;
        end
        fprintf('Multiple propellers found in folder:\n');
        for i = 1:numel(available)
            fprintf('  %d) %s\n', i, available{i});
        end
        fprintf('\nRe-run and choose which to analyze, for example:\n');
        fprintf('  analyze_apc_prop(folder, ''%s'')\n', available{1});
        if numel(available) >= 2
            fprintf('  analyze_apc_prop(folder, ''%s,%s'')\n', available{1}, available{2});
        end
        fprintf('  analyze_apc_prop(folder, ''all'')\n');
        return;
    end

    if ischar(propSpec) || isstring(propSpec)
        propSpec = char(propSpec);
        if strcmpi(strtrim(propSpec), 'all')
            propsToRun = available;
            return;
        end
        parts = strsplit(propSpec, {',', ' ', ';', '\t'});
        parts = parts(~cellfun(@isempty, strtrim(parts)));
        propsToRun = strtrim(parts);
    elseif iscell(propSpec)
        propsToRun = cellfun(@char, propSpec, 'UniformOutput', false);
    else
        error('propSpec must be a string, char, or cell array of names (or ''all'').');
    end

    missing = setdiff(propsToRun, available);
    if ~isempty(missing)
        error('Prop(s) not found in folder: %s\nAvailable: %s', ...
            strjoin(missing, ', '), strjoin(available, ', '));
    end
end

%% -------- per-prop analysis --------

function result = analyzeOneProp(dataFolder, propName, degreeJ, degreeRe, ReScale)
    perfPath = fullfile(dataFolder, [propName '_performance.dat']);
    geomPath = fullfile(dataFolder, [propName '_geometry.dat']);
    assert(isfile(perfPath), 'Missing performance file: %s', perfPath);

    raw = parsePropPerformance(perfPath);
    assert(strcmpi(raw.name, propName) || startsWith(raw.name, propName), ...
        'File prop name "%s" does not match requested "%s"', raw.name, propName);

    if isfile(geomPath)
        geom = parsePropGeometry(geomPath);
        chord = chordStats(geom);
    else
        warning('No geometry file for %s (%s). Chord fields will be null.', ...
            propName, geomPath);
        geom = struct('name', propName, 'Dia_in', raw.Dia_in, 'pitch_in', raw.pitch_in);
        chord = emptyChordStats();
    end

    J  = raw.J;
    Re = raw.Re;
    CT = raw.CT;
    CP = raw.CP;

    ReHat = Re / ReScale;
    [CT_coeffs, basis] = fitBivariatePoly(J, ReHat, CT, degreeJ, degreeRe);
    [CP_coeffs, ~]     = fitBivariatePoly(J, ReHat, CP, degreeJ, degreeRe);

    CTpred = evalBivariatePoly(J, ReHat, CT_coeffs, degreeJ, degreeRe);
    CPpred = evalBivariatePoly(J, ReHat, CP_coeffs, degreeJ, degreeRe);
    rmseCT = sqrt(mean((CT - CTpred).^2));
    rmseCP = sqrt(mean((CP - CPpred).^2));
    r2CT = 1 - sum((CT - CTpred).^2) / sum((CT - mean(CT)).^2);
    r2CP = 1 - sum((CP - CPpred).^2) / sum((CP - mean(CP)).^2);

    orderingText = sprintf([ ...
        'Coefficients are ordered to match basis[]. ' ...
        'Each entry multiplies J^i * Re_hat^j, with Re_hat = Re / %.0e. ' ...
        'CT = sum_k CT_coefficients[k] * J^(basis[k].i) * Re_hat^(basis[k].j). ' ...
        'Same ordering for CP_coefficients.'], ReScale);

    result = struct();
    result.schema_version = '1.0.0';
    result.prop_id = sprintf('APC-%s', propName);
    result.manufacturer = 'APC';
    result.model = propName;
    result.diameter_in = raw.Dia_in;
    result.diameter_m = raw.Dia_in * 0.0254;
    result.pitch_in = raw.pitch_in;
    result.blades = getfielddef(geom, 'blades', NaN);

    result.source = struct();
    result.source.data_folder = dataFolder;
    result.source.performance_file = [propName '_performance.dat'];
    result.source.performance_path = perfPath;
    if isfile(geomPath)
        result.source.geometry_file = [propName '_geometry.dat'];
        result.source.geometry_path = geomPath;
    else
        result.source.geometry_file = '';
        result.source.geometry_path = '';
    end

    result.fit_generated_utc = char(datetime('now', 'TimeZone', 'UTC', ...
        'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z'''));

    result.fit_model = struct();
    result.fit_model.type = 'bivariate_polynomial';
    result.fit_model.description = ...
        'Response surface: C(J, Re_hat) = sum_{i=0}^{nJ} sum_{j=0}^{nRe} a_ij * J^i * Re_hat^j';
    result.fit_model.degree_J = degreeJ;
    result.fit_model.degree_Re = degreeRe;
    result.fit_model.Re_scaling = struct( ...
        'method', 'divide', ...
        'scale', ReScale, ...
        'formula', sprintf('Re_hat = Re / %.0e', ReScale));
    result.fit_model.coefficient_ordering = orderingText;
    result.fit_model.basis = basis;
    result.fit_model.CT_coefficients = CT_coeffs(:).';
    result.fit_model.CP_coefficients = CP_coeffs(:).';
    result.fit_model.n_points_used = numel(J);
    result.fit_model.rmse_CT = rmseCT;
    result.fit_model.rmse_CP = rmseCP;
    result.fit_model.r2_CT = r2CT;
    result.fit_model.r2_CP = r2CP;

    result.validity = struct();
    result.validity.J_min = min(J);
    result.validity.J_max = max(J);
    result.validity.Re_min = min(Re);
    result.validity.Re_max = max(Re);
    result.validity.note = ...
        'Do not extrapolate outside these ranges. Re is APC-reported 75%-span Reynolds number.';

    result.units_and_conventions = struct();
    result.units_and_conventions.J = 'V/(n*D), advance ratio (dimensionless)';
    result.units_and_conventions.CT = 'T/(rho*n^2*D^4); positive thrust along prop axis';
    result.units_and_conventions.CP = 'P/(rho*n^3*D^5); positive when prop absorbs shaft power';
    result.units_and_conventions.Re = 'Reynolds number at 75% blade radius from APC performance tables';
    result.units_and_conventions.n = 'rotational speed in revolutions per second';
    result.units_and_conventions.D = 'propeller diameter; keep length units consistent with V';
    result.units_and_conventions.rho = 'air density consistent with force/power units';
    result.units_and_conventions.chord = 'meters in *_m fields; inches in *_in fields';
    result.units_and_conventions.sign_convention = ...
        'APC table CT/CP signs retained; static CT>0 and CP>0 for thrusting propeller';

    result.geometry = struct();
    result.geometry.radius_in = getfielddef(geom, 'radius_in', raw.Dia_in / 2);
    result.geometry.radius_m = result.geometry.radius_in * 0.0254;
    result.geometry.hub_radius_in = getfielddef(geom, 'hub_radius_in', NaN);
    result.geometry.hub_transition_in = getfielddef(geom, 'hub_transition_in', NaN);
    result.geometry.chord_min_in = chord.chord_min_in;
    result.geometry.chord_max_in = chord.chord_max_in;
    result.geometry.chord_avg_in = chord.chord_avg_in;
    result.geometry.chord_75R_in = chord.chord_75R_in;
    result.geometry.chord_min_m = chord.chord_min_m;
    result.geometry.chord_max_m = chord.chord_max_m;
    result.geometry.chord_avg_m = chord.chord_avg_m;
    result.geometry.chord_75R_m = chord.chord_75R_m;
    result.geometry.weight_kg = getfielddef(geom, 'weight_kg', NaN);
    result.geometry.moi_kg_m2 = getfielddef(geom, 'moi_kg_m2', NaN);
    result.geometry.activity_factor = getfielddef(geom, 'activity_factor', NaN);
    result.geometry.projected_area_in2 = getfielddef(geom, 'projected_area_in2', NaN);
end

%% -------- bivariate polynomial fit --------

function [coeffs, basis] = fitBivariatePoly(x, y, z, nX, nY)
    % Build basis list and design matrix for z ~= sum a_ij x^i y^j
    nTerms = (nX + 1) * (nY + 1);
    basis = repmat(struct('i', 0, 'j', 0), nTerms, 1);
    A = zeros(numel(z), nTerms);
    k = 0;
    for j = 0:nY
        for i = 0:nX
            k = k + 1;
            basis(k).i = i;
            basis(k).j = j;
            A(:, k) = (x.^i) .* (y.^j);
        end
    end
    coeffs = A \ z;   % least squares
end

function z = evalBivariatePoly(x, y, coeffs, nX, nY)
    z = zeros(size(x));
    k = 0;
    for j = 0:nY
        for i = 0:nX
            k = k + 1;
            z = z + coeffs(k) * (x.^i) .* (y.^j);
        end
    end
end

%% -------- geometry / chord helpers --------

function chord = chordStats(geom)
    chord = emptyChordStats();
    if ~isfield(geom, 'station_in') || isempty(geom.station_in)
        return;
    end

    blade = true(size(geom.chord_in));
    if isfield(geom, 'hub_transition_in') && ~isnan(geom.hub_transition_in)
        blade = blade & (geom.station_in >= geom.hub_transition_in);
    end
    if isfield(geom, 'radius_in') && ~isnan(geom.radius_in)
        blade = blade & (geom.station_in <= 0.95 * geom.radius_in);
        c75 = interp1(geom.station_in, geom.chord_in, 0.75 * geom.radius_in, ...
            'linear', 'extrap');
    else
        c75 = NaN;
    end

    c = geom.chord_in(blade);
    chord.chord_min_in = min(c);
    chord.chord_max_in = max(c);
    chord.chord_avg_in = mean(c);
    chord.chord_75R_in = c75;
    chord.chord_min_m = chord.chord_min_in * 0.0254;
    chord.chord_max_m = chord.chord_max_in * 0.0254;
    chord.chord_avg_m = chord.chord_avg_in * 0.0254;
    chord.chord_75R_m = chord.chord_75R_in * 0.0254;
end

function chord = emptyChordStats()
    chord = struct( ...
        'chord_min_in', NaN, 'chord_max_in', NaN, ...
        'chord_avg_in', NaN, 'chord_75R_in', NaN, ...
        'chord_min_m', NaN, 'chord_max_m', NaN, ...
        'chord_avg_m', NaN, 'chord_75R_m', NaN);
end

function v = getfielddef(s, name, default)
    if isstruct(s) && isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end

%% -------- parsers --------

function raw = parsePropPerformance(filePath)
    txt = fileread(filePath);
    lines = splitlines(txt);
    assert(contains(txt, 'PROP RPM'), ...
        '%s does not look like a performance file.', filePath);

    name = parsePropName(lines, filePath);
    [Dia_in, pitch_in] = parsePropSize(name);

    J = []; CT = []; CP = []; Re = []; RPM = [];
    currentRPM = NaN;

    for i = 1:numel(lines)
        line = lines{i};
        rpmTok = regexp(line, 'PROP\s+RPM\s*=\s*([0-9.]+)', 'tokens', 'once');
        if ~isempty(rpmTok)
            currentRPM = str2double(rpmTok{1});
            continue;
        end
        nums = sscanf(strtrim(line), '%f');
        if numel(nums) >= 14 && ~isnan(currentRPM)
            J(end+1,1) = nums(2); %#ok<AGROW>
            CT(end+1,1) = nums(4); %#ok<AGROW>
            CP(end+1,1) = nums(5); %#ok<AGROW>
            Re(end+1,1) = nums(14); %#ok<AGROW>
            RPM(end+1,1) = currentRPM; %#ok<AGROW>
        end
    end
    assert(~isempty(J), 'No performance rows in %s', filePath);

    raw = struct('name', name, 'Dia_in', Dia_in, 'pitch_in', pitch_in, ...
        'J', J, 'CT', CT, 'CP', CP, 'Re', Re, 'RPM', RPM);
end

function geom = parsePropGeometry(filePath)
    txt = fileread(filePath);
    lines = splitlines(txt);
    assert(contains(txt, 'AIRFOIL SUMMARY') || contains(txt, 'INERTIA AND AREA'), ...
        '%s does not look like a geometry file.', filePath);

    name = parsePropName(lines, filePath);
    [Dia_in, pitch_in] = parsePropSize(name);

    geom = struct();
    geom.name = name;
    geom.Dia_in = Dia_in;
    geom.pitch_in = pitch_in;
    geom.radius_in = NaN;
    geom.hub_radius_in = NaN;
    geom.hub_transition_in = NaN;
    geom.blades = NaN;
    geom.weight_kg = NaN;
    geom.moi_kg_m2 = NaN;
    geom.activity_factor = NaN;
    geom.projected_area_in2 = NaN;

    station = [];
    chord = [];

    for i = 1:numel(lines)
        line = lines{i};

        tok = regexp(line, 'RADIUS:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.radius_in = str2double(tok{1}); end
        tok = regexp(line, 'HUBRAD:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.hub_radius_in = str2double(tok{1}); end
        tok = regexp(line, 'HUBTRA:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.hub_transition_in = str2double(tok{1}); end
        tok = regexp(line, 'BLADES:\s*([0-9]+)', 'tokens', 'once');
        if ~isempty(tok), geom.blades = str2double(tok{1}); end
        tok = regexp(line, 'TOTAL WEIGHT \(Kg\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.weight_kg = str2double(tok{1}); end
        tok = regexp(line, 'MOMENT OF INERTIA \(Kg-M\*\*2\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.moi_kg_m2 = str2double(tok{1}); end
        tok = regexp(line, 'ACTIVITY FACTOR\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.activity_factor = str2double(tok{1}); end
        tok = regexp(line, 'TOTAL PROJECTED AREA \(IN\*\*2\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.projected_area_in2 = str2double(tok{1}); end

        nums = sscanf(strtrim(line), '%f');
        if numel(nums) >= 2 && contains(line, '.')
            if contains(upper(line), {'RADIUS', 'WEIGHT', 'INERTIA', 'FACTOR', 'LIMIT', 'HUB'})
                continue;
            end
            station(end+1,1) = nums(1); %#ok<AGROW>
            chord(end+1,1) = nums(2); %#ok<AGROW>
        end
    end

    geom.station_in = station;
    geom.chord_in = chord;
end

function name = parsePropName(lines, filePath)
    name = '';
    for i = 1:numel(lines)
        tok = regexp(strtrim(lines{i}), '^(\d+x\d+\w*)', 'tokens', 'once');
        if ~isempty(tok)
            name = tok{1};
            break;
        end
    end
    assert(~isempty(name), 'Could not find prop name in %s', filePath);
end

function [Dia_in, pitch_in] = parsePropSize(name)
    geom = regexp(name, '^(\d+)x(\d+)', 'tokens', 'once');
    assert(~isempty(geom), 'Could not parse diameter/pitch from "%s"', name);
    Dia_in = str2double(geom{1});
    pitch_in = str2double(geom{2});
end

%% -------- JSON I/O --------

function writeJson(outPath, result)
    % jsonencode does not like NaN well in all releases; convert NaN -> null via post-process
    txt = jsonencode(result, 'PrettyPrint', true);
    txt = regexprep(txt, '\bNaN\b', 'null');
    txt = regexprep(txt, '\bInf\b', 'null');
    fid = fopen(outPath, 'w');
    assert(fid > 0, 'Could not write %s', outPath);
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, txt, 'char');
end
