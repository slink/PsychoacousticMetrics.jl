% run_sqat_pa.m — Octave driver for the psychoacoustic-annoyance oracle rig.
%
% Runs the pinned SQAT reference implementation
% (github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb,
% CC BY-NC) in two modes, and vendors the numeric outputs as test
% fixtures. This script only CALLS the reference; no SQAT source is
% reproduced or shipped by this file. The checkout is never modified —
% every path used below is added via `addpath` (no shims needed, same as
% run_sqat_fs.m).
%
% Generated/run under: Octave (OCTAVE_VERSION), octave-signal package
% (see `pkg load signal; ver('signal')`). Both are printed to stdout by
% every invocation for provenance.
%
% Usage:
%   octave --no-gui run_sqat_pa.m formula <manifest.json> <output.json> [sqat_dir]
%   octave --no-gui run_sqat_pa.m signal  <manifest.json> <output_dir>  [sqat_dir]
%
% formula mode:
%   Calls PsychoacousticAnnoyance_Widmann1992_from_percentile(N,S,R,FS)
%   once per row of the manifest (the reference function is written for
%   scalar arguments — see pins file). Manifest: JSON array of objects
%   {"N": ..., "S": ..., "R": ..., "FS": ...}. Output: single JSON object
%   {"N": [...], "S": [...], "R": [...], "FS": [...], "pa": [...]}
%   (arrays in manifest order).
%
% signal mode:
%   Calls PsychoacousticAnnoyance_Widmann1992(insig, fs, LoudnessField=0,
%   time_skip=0, showPA=0, show=0) per case (LoudnessField=0 = free
%   field, matching this package's ZwickerLoudness default field_type =
%   :free; time_skip=0 because all cross-check signals are steady-state
%   with no onset transient to discard). Manifest: JSON array of objects
%   {"name": ..., "fs": ..., "signal_path": "<float64 LE .bin>"}. Per-case
%   output: <output_dir>/<name>.json with fields
%   {name, fs, pa, N5, S5, R5, FS5} (pa = OUT.ScalarPA; N5/S5/R5/FS5 =
%   OUT.L.N5 / OUT.S.S5 / OUT.R.R5 / OUT.FS.FS5, the percentile values SQAT
%   feeds its own formula call — vendored so any deviation is attributable
%   per-component before any wrapper tolerance is set).

args = argv();
if numel(args) < 3
    error('usage: octave run_sqat_pa.m <formula|signal> <manifest.json> <output> [sqat_dir]');
end
mode = args{1};
manifest_path = args{2};
output_path = args{3};
if numel(args) >= 4
    sqat_dir = args{4};
else
    sqat_dir = '/tmp/sqat-pinned';
end

pkg load signal;
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'PsychoacousticAnnoyance_Widmann1992'));
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'Loudness_ISO532_1'));
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'Sharpness_DIN45692'));
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'Roughness_Daniel1997'));
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'FluctuationStrength_Osses2016'));
addpath(fullfile(sqat_dir, 'utilities'));
addpath(fullfile(sqat_dir, 'sound_level_meter'));

printf('# SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb, Octave %s, octave-signal %s, %s\n', ...
       OCTAVE_VERSION, ver('signal').Version, datestr(now));

manifest = jsondecode(fileread(manifest_path));
if ~iscell(manifest)
    manifest = num2cell(manifest); % jsondecode gives a struct array for uniform objects
end

if strcmp(mode, 'formula')

    n = numel(manifest);
    N_out = zeros(n, 1); S_out = zeros(n, 1); R_out = zeros(n, 1);
    FS_out = zeros(n, 1); pa_out = zeros(n, 1);
    for i = 1:n
        c = manifest{i};
        N_out(i) = c.N; S_out(i) = c.S; R_out(i) = c.R; FS_out(i) = c.FS;
        pa_out(i) = PsychoacousticAnnoyance_Widmann1992_from_percentile(c.N, c.S, c.R, c.FS);
    end
    printf('formula mode: %d cases\n', n);

    result = struct();
    result.N = N_out; result.S = S_out; result.R = R_out; result.FS = FS_out;
    result.pa = pa_out;
    fid = fopen(output_path, 'w');
    fprintf(fid, '%s', jsonencode(result));
    fclose(fid);

elseif strcmp(mode, 'signal')

    output_dir = output_path;
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    for i = 1:numel(manifest)
        c = manifest{i};
        name = c.name;
        fs = c.fs;

        fid = fopen(c.signal_path, 'rb');
        insig = fread(fid, Inf, 'double');
        fclose(fid);

        printf('running %s: N=%d fs=%g\n', name, numel(insig), fs);

        OUT = PsychoacousticAnnoyance_Widmann1992(insig, fs, 0, 0, 0, 0);

        result = struct();
        result.name = name;
        result.fs = fs;
        result.pa = OUT.ScalarPA;
        result.N5 = OUT.L.N5;
        result.S5 = OUT.S.S5;
        result.R5 = OUT.R.R5;
        result.FS5 = OUT.FS.FS5;

        out_path = fullfile(output_dir, [name '.json']);
        fid = fopen(out_path, 'w');
        fprintf(fid, '%s', jsonencode(result));
        fclose(fid);
    end

else
    error('unknown mode: %s (expected formula or signal)', mode);
end
