% Propeller CT / CP Polynomial Fits
% Reads APC/UIUC-style .dat files from prop_data, extracts advance-ratio
% performance data, fits polynomials for C_T(J) and C_P(J), plots
% raw data vs. fit, and saves results for use in other MATLAB scripts.
%
% Expected filenames in prop_data/:
%   <name>_performance.dat   - CT/CP vs J tables (required for fitting)
%   <name>_geometry.dat      - airfoil / inertia summary (optional)

clear; clc; close all;

% User settings
n = 5;              % polynomial order for C_T(J) and C_P(J) fits
Re_min = 0;         % keep points with Reyn >= Re_min  (set 0 for no lower bound)
Re_max = Inf;       % keep points with Reyn <= Re_max  (set Inf for no upper bound)

% Paths and file discovery
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
dataDir = fullfile(scriptDir, 'prop_data');
outFile = fullfile(scriptDir, 'prop_polyfit_data.mat');

perfFiles = dir(fullfile(dataDir, '*_performance.dat'));
assert(~isempty(perfFiles), ...
    'No *_performance.dat files found in %s', dataDir);

fprintf('Found %d performance file(s) in %s\n', numel(perfFiles), dataDir);
fprintf('Using polynomial order n = %d\n', n);
fprintf('Reynolds filter: %.0f <= Re <= %.0f\n', Re_min, Re_max);

% Process each propeller
props = struct([]);

for k = 1:numel(perfFiles)
    perfPath = fullfile(perfFiles(k).folder, perfFiles(k).name);
    propId = erase(perfFiles(k).name, '_performance.dat');
    geomPath = fullfile(dataDir, [propId '_geometry.dat']);

    fprintf('\nReading %s ...\n', perfFiles(k).name);
    raw = parsePropPerformance(perfPath);

    geom = struct();
    if isfile(geomPath)
        fprintf('  + geometry %s\n', [propId '_geometry.dat']);
        geom = parsePropGeometry(geomPath);
    else
        fprintf('  (no matching *_geometry.dat)\n');
    end

    % Filter by Reynolds number (75% span, as reported in the .dat files)
    keep = (raw.Re >= Re_min) & (raw.Re <= Re_max);
    J  = raw.J(keep);
    CT = raw.CT(keep);
    CP = raw.CP(keep);
    Re = raw.Re(keep);
    RPM = raw.RPM(keep);

    assert(numel(J) > n, ...
        '%s: only %d points after Re filter (need > n=%d). Widen Re_min/Re_max.', ...
        raw.name, numel(J), n);

    % Degree-n fits: C_T(J) and C_P(J)
    % polyfit returns [a_n ... a1 a0] for a_n*J^n + ... + a1*J + a0
    CT_poly = polyfit(J, CT, n);
    CP_poly = polyfit(J, CP, n);

    props(k).name     = raw.name;
    props(k).file     = perfFiles(k).name;
    props(k).geom_file = '';
    if isfile(geomPath)
        props(k).geom_file = [propId '_geometry.dat'];
    end
    props(k).Dia_in   = raw.Dia_in;          % inches (from DxP name)
    props(k).Dia_ft   = raw.Dia_in / 12;     % feet
    props(k).Dia_m    = raw.Dia_in * 0.0254; % meters
    props(k).pitch_in = raw.pitch_in;
    props(k).poly_n   = n;
    props(k).Re_min   = Re_min;
    props(k).Re_max   = Re_max;
    props(k).J        = J;
    props(k).CT       = CT;
    props(k).CP       = CP;
    props(k).Re       = Re;
    props(k).RPM      = RPM;
    props(k).CT_poly  = CT_poly;
    props(k).CP_poly  = CP_poly;
    props(k).J_range  = [min(J), max(J)];
    props(k).Re_range = [min(Re), max(Re)];
    props(k).geometry = geom;

    fprintf('  %s | Dia = %.2f in (%.4f m) | kept %d / %d pts | J in [%.3f, %.3f] | Re in [%.0f, %.0f]\n', ...
        props(k).name, props(k).Dia_in, props(k).Dia_m, numel(J), numel(raw.J), ...
        props(k).J_range(1), props(k).J_range(2), ...
        props(k).Re_range(1), props(k).Re_range(2));
    fprintf('  CT poly (n=%d): %s\n', n, mat2str(CT_poly, 6));
    fprintf('  CP poly (n=%d): %s\n', n, mat2str(CP_poly, 6));
end

% Plot C_T(J) - filtered data and fit
figure('Name', sprintf('C_T(J) order-%d fits', n), 'Color', 'w');
tiledlayout(1, numel(props), 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:numel(props)
    nexttile;
    Jfit = linspace(props(k).J_range(1), props(k).J_range(2), 200);
    CTfit = polyval(props(k).CT_poly, Jfit);

    plot(props(k).J, props(k).CT, '.', 'MarkerSize', 6, 'DisplayName', 'filtered');
    hold on;
    plot(Jfit, CTfit, 'LineWidth', 2, 'DisplayName', sprintf('n = %d fit', n));
    hold off;
    grid on;
    xlabel('J');
    ylabel('C_T');
    title(sprintf('%s  (D = %.0f in)', props(k).name, props(k).Dia_in));
    legend('Location', 'best');
end
sgtitle(sprintf('Thrust coefficient C_T(J)  |  n = %d  |  %.0f <= Re <= %.0f', n, Re_min, Re_max));

% Plot C_P(J) - filtered data and fit
figure('Name', sprintf('C_P(J) order-%d fits', n), 'Color', 'w');
tiledlayout(1, numel(props), 'Padding', 'compact', 'TileSpacing', 'compact');

for k = 1:numel(props)
    nexttile;
    Jfit = linspace(props(k).J_range(1), props(k).J_range(2), 200);
    CPfit = polyval(props(k).CP_poly, Jfit);

    plot(props(k).J, props(k).CP, '.', 'MarkerSize', 6, 'DisplayName', 'filtered');
    hold on;
    plot(Jfit, CPfit, 'LineWidth', 2, 'DisplayName', sprintf('n = %d fit', n));
    hold off;
    grid on;
    xlabel('J');
    ylabel('C_P');
    title(sprintf('%s  (D = %.0f in)', props(k).name, props(k).Dia_in));
    legend('Location', 'best');
end
sgtitle(sprintf('Power coefficient C_P(J)  |  n = %d  |  %.0f <= Re <= %.0f', n, Re_min, Re_max));

% Save results for later MATLAB use
% Load later with:
%   S = load('prop_polyfit_data.mat');
%   props = S.props;
%   CT = polyval(props(1).CT_poly, J);
%   CP = polyval(props(1).CP_poly, J);
%   D  = props(1).Dia_m;   % or Dia_in / Dia_ft
save(outFile, 'props', 'n', 'Re_min', 'Re_max');
fprintf('\nSaved %d propeller(s) (n = %d, Re filter) to:\n  %s\n', numel(props), n, outFile);

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
    assert(~isempty(geom), 'Could not parse diameter/pitch from name "%s"', name);
    Dia_in = str2double(geom{1});
    pitch_in = str2double(geom{2});
end

function raw = parsePropPerformance(filePath)
%PARSEPROPPERFORMANCE Parse APC performance tables (CT/CP vs J by RPM).

    txt = fileread(filePath);
    lines = splitlines(txt);

    assert(contains(txt, 'PROP RPM'), ...
        '%s does not look like a performance file (missing PROP RPM blocks).', filePath);

    name = parsePropName(lines, filePath);
    [Dia_in, pitch_in] = parsePropSize(name);

    J = [];
    CT = [];
    CP = [];
    Re = [];
    RPM = [];
    currentRPM = NaN;

    for i = 1:numel(lines)
        line = lines{i};

        rpmTok = regexp(line, 'PROP\s+RPM\s*=\s*([0-9.]+)', 'tokens', 'once');
        if ~isempty(rpmTok)
            currentRPM = str2double(rpmTok{1});
            continue;
        end

        % Numeric data rows: V, J, Pe, Ct, Cp, ..., Mach, Reyn, FOM
        nums = sscanf(strtrim(line), '%f');
        if numel(nums) >= 14 && ~isnan(currentRPM)
            J(end+1,1)   = nums(2); %#ok<AGROW>
            CT(end+1,1)  = nums(4); %#ok<AGROW>
            CP(end+1,1)  = nums(5); %#ok<AGROW>
            Re(end+1,1)  = nums(14); %#ok<AGROW>  % Reyn at 75% span
            RPM(end+1,1) = currentRPM; %#ok<AGROW>
        end
    end

    assert(~isempty(J), 'No performance data rows found in %s', filePath);

    raw.name     = name;
    raw.Dia_in   = Dia_in;
    raw.pitch_in = pitch_in;
    raw.J        = J;
    raw.CT       = CT;
    raw.CP       = CP;
    raw.Re       = Re;
    raw.RPM      = RPM;
end

function geom = parsePropGeometry(filePath)
%PARSEPROPGEOMETRY Parse APC geometry / airfoil / inertia summary file.

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
    geom.weight_lb = NaN;
    geom.weight_kg = NaN;
    geom.volume_in3 = NaN;
    geom.projected_area_in2 = NaN;
    geom.moi_kg_m2 = NaN;
    geom.activity_factor = NaN;

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

        tok = regexp(line, 'TOTAL WEIGHT \(LB\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.weight_lb = str2double(tok{1}); end

        tok = regexp(line, 'TOTAL WEIGHT \(Kg\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.weight_kg = str2double(tok{1}); end

        tok = regexp(line, 'TOTAL VOLUME \(IN\*\*3\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.volume_in3 = str2double(tok{1}); end

        tok = regexp(line, 'TOTAL PROJECTED AREA \(IN\*\*2\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.projected_area_in2 = str2double(tok{1}); end

        tok = regexp(line, 'MOMENT OF INERTIA \(Kg-M\*\*2\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.moi_kg_m2 = str2double(tok{1}); end

        tok = regexp(line, 'ACTIVITY FACTOR\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), geom.activity_factor = str2double(tok{1}); end
    end
end
