% Prop chord / geometry stats from *_geometry.dat files
% Useful for Reynolds number estimates: Re = V_tip * c / nu  (or omega*r*c/nu)

clear; clc;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir), scriptDir = pwd; end
dataDir = fullfile(scriptDir, 'prop_data');

files = dir(fullfile(dataDir, '*_geometry.dat'));
assert(~isempty(files), 'No *_geometry.dat files in %s', dataDir);

props = struct([]);

for k = 1:numel(files)
    p = parseGeometry(fullfile(files(k).folder, files(k).name));

    % Chord over the working blade: outside hub transition, before tip taper
    blade = true(size(p.chord_in));
    if ~isnan(p.hub_transition_in)
        blade = blade & (p.station_in >= p.hub_transition_in);
    end
    if ~isnan(p.radius_in)
        blade = blade & (p.station_in <= 0.95 * p.radius_in);
    end
    c_blade = p.chord_in(blade);

    % Chord at 75% radius (common Re reference station)
    r75 = 0.75 * p.radius_in;
    c75 = interp1(p.station_in, p.chord_in, r75, 'linear', 'extrap');

    props(k).name              = p.name;
    props(k).file              = files(k).name;
    props(k).Dia_in            = p.Dia_in;
    props(k).Dia_m             = p.Dia_in * 0.0254;
    props(k).radius_in         = p.radius_in;
    props(k).radius_m          = p.radius_in * 0.0254;
    props(k).hub_radius_in     = p.hub_radius_in;
    props(k).hub_transition_in = p.hub_transition_in;
    props(k).blades            = p.blades;
    props(k).chord_min_in      = min(c_blade);
    props(k).chord_max_in      = max(c_blade);
    props(k).chord_avg_in      = mean(c_blade);
    props(k).chord_75R_in      = c75;
    props(k).chord_min_m       = props(k).chord_min_in * 0.0254;
    props(k).chord_max_m       = props(k).chord_max_in * 0.0254;
    props(k).chord_avg_m       = props(k).chord_avg_in * 0.0254;
    props(k).chord_75R_m       = props(k).chord_75R_in * 0.0254;
    props(k).weight_kg         = p.weight_kg;
    props(k).moi_kg_m2         = p.moi_kg_m2;
    props(k).activity_factor   = p.activity_factor;
    props(k).station_in        = p.station_in;
    props(k).chord_in          = p.chord_in;

    fprintf('\n=== %s ===\n', props(k).name);
    fprintf('  Dia = %.2f in (%.4f m) | R = %.2f in | blades = %g\n', ...
        props(k).Dia_in, props(k).Dia_m, props(k).radius_in, props(k).blades);
    fprintf('  Chord (blade, in):  min = %.4f | max = %.4f | avg = %.4f | @0.75R = %.4f\n', ...
        props(k).chord_min_in, props(k).chord_max_in, props(k).chord_avg_in, props(k).chord_75R_in);
    fprintf('  Chord (blade, m):   min = %.5f | max = %.5f | avg = %.5f | @0.75R = %.5f\n', ...
        props(k).chord_min_m, props(k).chord_max_m, props(k).chord_avg_m, props(k).chord_75R_m);
    fprintf('  Weight = %.4f kg | MoI = %.3e kg-m^2 | AF = %.1f\n', ...
        props(k).weight_kg, props(k).moi_kg_m2, props(k).activity_factor);
end

outFile = fullfile(scriptDir, 'prop_chord_stats.mat');
save(outFile, 'props');
fprintf('\nSaved to %s\n', outFile);

function p = parseGeometry(filePath)
    lines = splitlines(fileread(filePath));

    p.name = '';
    for i = 1:numel(lines)
        tok = regexp(strtrim(lines{i}), '^(\d+x\d+\w*)', 'tokens', 'once');
        if ~isempty(tok), p.name = tok{1}; break; end
    end
    assert(~isempty(p.name), 'No prop name in %s', filePath);

    sz = regexp(p.name, '^(\d+)x(\d+)', 'tokens', 'once');
    p.Dia_in = str2double(sz{1});
    p.pitch_in = str2double(sz{2});

    p.radius_in = NaN;
    p.hub_radius_in = NaN;
    p.hub_transition_in = NaN;
    p.blades = NaN;
    p.weight_kg = NaN;
    p.moi_kg_m2 = NaN;
    p.activity_factor = NaN;

    station = [];
    chord = [];

    for i = 1:numel(lines)
        line = lines{i};

        tok = regexp(line, 'RADIUS:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), p.radius_in = str2double(tok{1}); end

        tok = regexp(line, 'HUBRAD:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), p.hub_radius_in = str2double(tok{1}); end

        tok = regexp(line, 'HUBTRA:\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), p.hub_transition_in = str2double(tok{1}); end

        tok = regexp(line, 'BLADES:\s*([0-9]+)', 'tokens', 'once');
        if ~isempty(tok), p.blades = str2double(tok{1}); end

        tok = regexp(line, 'TOTAL WEIGHT \(Kg\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), p.weight_kg = str2double(tok{1}); end

        tok = regexp(line, 'MOMENT OF INERTIA \(Kg-M\*\*2\)\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), p.moi_kg_m2 = str2double(tok{1}); end

        tok = regexp(line, 'ACTIVITY FACTOR\s*=\s*([0-9.eE+-]+)', 'tokens', 'once');
        if ~isempty(tok), p.activity_factor = str2double(tok{1}); end

        % Station table rows: station(in), chord(in), ...
        nums = sscanf(strtrim(line), '%f');
        if numel(nums) >= 2 && contains(line, '.')
            % Skip summary lines that also parse as numbers
            if contains(upper(line), {'RADIUS', 'WEIGHT', 'INERTIA', 'FACTOR', 'LIMIT', 'HUB'})
                continue;
            end
            station(end+1,1) = nums(1); %#ok<AGROW>
            chord(end+1,1)   = nums(2); %#ok<AGROW>
        end
    end

    assert(~isempty(chord), 'No chord stations found in %s', filePath);
    p.station_in = station;
    p.chord_in = chord;
end
