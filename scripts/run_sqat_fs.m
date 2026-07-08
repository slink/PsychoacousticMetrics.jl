% run_sqat_fs.m — Octave driver for the fluctuation-strength oracle rig.
%
% Runs the pinned SQAT reference implementation
% (github.com/ggrecow/SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb,
% CC BY-NC) on signals synthesized by this package's own test generator
% (test/support/fs_stimuli.jl), and vendors the numeric outputs as test
% fixtures. This script only CALLS the reference; no SQAT source is
% reproduced or shipped by this file. The checkout is never modified —
% every path used below is added via `addpath` (no shims were needed).
%
% Generated/run under: Octave (OCTAVE_VERSION), octave-signal package
% (see `pkg load signal; ver('signal')`). Both are printed to stdout by
% every invocation for provenance.
%
% Usage:
%   octave --no-gui run_sqat_fs.m <manifest.json> <output_dir> [sqat_dir]
%
% Manifest: JSON array of objects
%   {"name": ..., "fs": ..., "method": ..., "signal_path": "<float64 LE .bin>"}
%
% Per-case output: <output_dir>/<name>.json with fields
%   {name, fs, method, fluct, t, fs_mean, specific_tavg}

args = argv();
if numel(args) < 2
    error('usage: octave run_sqat_fs.m <manifest.json> <output_dir> [sqat_dir]');
end
manifest_path = args{1};
output_dir = args{2};
if numel(args) >= 3
    sqat_dir = args{3};
else
    sqat_dir = '/tmp/sqat-pinned';
end

pkg load signal;
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'FluctuationStrength_Osses2016'));
addpath(fullfile(sqat_dir, 'utilities'));

printf('# SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb, Octave %s, octave-signal %s, %s\n', ...
       OCTAVE_VERSION, ver('signal').Version, datestr(now));

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

manifest = jsondecode(fileread(manifest_path));
if ~iscell(manifest)
    manifest = num2cell(manifest); % jsondecode gives a struct array for uniform objects
end

for i = 1:numel(manifest)
    c = manifest{i};
    name = c.name;
    fs = c.fs;
    method = c.method;

    fid = fopen(c.signal_path, 'rb');
    insig = fread(fid, Inf, 'double');
    fclose(fid);

    printf('running %s: N=%d fs=%g method=%d\n', name, numel(insig), fs, method);

    OUT = FluctuationStrength_Osses2016(insig, fs, method, 0);

    result = struct();
    result.name = name;
    result.fs = fs;
    result.method = method;
    result.fluct = OUT.InstantaneousFluctuationStrength(:);
    result.t = OUT.time(:);
    result.fs_mean = mean(OUT.InstantaneousFluctuationStrength);
    result.specific_tavg = OUT.TimeAveragedSpecificFluctuationStrength(:);

    out_path = fullfile(output_dir, [name '.json']);
    fid = fopen(out_path, 'w');
    fprintf(fid, '%s', jsonencode(result));
    fclose(fid);
end
