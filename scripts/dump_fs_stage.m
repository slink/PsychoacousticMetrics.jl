% dump_fs_stage.m — instrumented oracle reproduction. Calls the pinned SQAT
% checkout's own functions via addpath; no SQAT code is copied here.
%
% Produces stage_dumps.json: frame-1 intermediate quantities of the
% fluctuation-strength pipeline for the anchor_44k stimulus, so that the
% from-scratch Julia implementation (Tasks 2-4) can check its own
% intermediate stages against the running oracle.
%
% Provenance/approach: every reusable reference function is CALLED from
% the pinned SQAT checkout (cos_ramp, calculate_a0,
% TerhardtExcitationPatterns, Get_Hweight_fluctuation,
% Get_gzi_fluctuation — reachable via addpath, incl. the private/ dir,
% since Octave does not enforce MATLAB's private-directory scoping).
% The glue between those calls below encodes PINNED FACTS about the
% pipeline (recorded in docs/oracle-pins.md §2) in this script's own wording:
%  - frame windowing: 50 ms cosine on/off ramps (cos_ramp);
%  - a0 stage: linear-phase FIR (order K=2^12), group delay K/2 samples
%    compensated by zero-padding and dropping the first K/2 outputs;
%  - envelope stats per critical band: h0 = time-mean of |excitation|,
%    fluctuation = Hweight-band-passed (h0-removed) envelope, modulation
%    depth = rms(fluctuation)/h0 (0 where h0 == 0);
%  - modulation-depth compression: depths above 0.7 are mapped to
%    0.7 + 0.3*(excess), then clamped to <= 1  [Osses et al. 2016,
%    doi:10.1121/2.0000410, "Version 3" FM correction];
%  - inter-band similarity: Pearson correlation between band k and k+2
%    fluctuation envelopes (0 where either variance is 0); boundary
%    channels are nearest-neighbor COPIES — pinned fact from
%    docs/oracle-pins.md §2.4: under Octave the reference's boundary
%    extrapolation deterministically resolves to
%      ki1[46] = ki1[47] = ki1[45]  and  ki2[1] = ki2[2] = ki1[3],
%    with ki2[k] = ki1[k-2] for k = 3..47 — encoded here directly as
%    those assignments;
%  - specific FS: fi = 0.4980 * gzi .* md.^1.7 .* |k1.*k2|.^1.7
%    .* sign(k1.*k2)  [published model equation, cal doubled 15/06/2016
%    per the reference's parameter table].
%
% Usage:
%   octave --no-gui dump_fs_stage.m <signal.bin> <fs> <out.json> [sqat_dir]

args = argv();
if numel(args) < 3
    error('usage: octave dump_fs_stage.m <signal.bin> <fs> <out.json> [sqat_dir]');
end
signal_path = args{1};
fs = str2double(args{2});
out_path = args{3};
if numel(args) >= 4
    sqat_dir = args{4};
else
    sqat_dir = '/tmp/sqat-pinned';
end

pkg load signal;
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'FluctuationStrength_Osses2016'));
addpath(fullfile(sqat_dir, 'psychoacoustic_metrics', 'FluctuationStrength_Osses2016', 'private'));
addpath(fullfile(sqat_dir, 'utilities'));

fid = fopen(signal_path, 'rb');
insig = fread(fid, Inf, 'double');
fclose(fid);

% Frame 1: first 2 s (or the whole signal if shorter).
N = min(round(2 * fs), numel(insig));
frame = insig(1:N);

% Windowing: 50 ms cosine ramps at both ends.
ramp = cos_ramp(N, fs, 50, 50);
windowed = ramp(:) .* frame;

% a0 outer/middle-ear stage: linear-phase FIR, delay K/2 compensated.
K = 2^12;
B = calculate_a0(fs, K, 'fluctuationstrength_osses2016');
padded = [windowed(:).' zeros(1, K/2)];
post_a0 = filter(B, 1, padded);
post_a0 = post_a0(K/2+1:end);
post_a0_rms = sqrt(mean(post_a0.^2));

% Excitation patterns (47 half-Bark channels), dBFS = 94 (1 Pa <-> 94 dB).
ei = TerhardtExcitationPatterns(post_a0, fs, 94);

% Envelope statistics, time x channel orientation.
env = abs(ei).';               % time x 47
h0 = mean(env);                % 1 x 47, per-channel time-mean
Hweight = Get_Hweight_fluctuation(fs);
fluctenv = sosfilt(Hweight, env - h0);   % band-passed, h0 removed (broadcast)
fluctrms = rms(fluctenv, 1);
mdepth_pre = zeros(1, 47);
nz = h0 > 0;
mdepth_pre(nz) = fluctrms(nz) ./ h0(nz);

% Compression of modulation depth (pinned model equation, threshold 0.7).
mdepth_post = mdepth_pre;
over = mdepth_post > 0.7;
mdepth_post(over) = 0.7 + (1 - 0.7) * (mdepth_post(over) - 0.7);
mdepth_post = min(mdepth_post, 1);

% Pearson correlation between channels k and k+2.
ki1 = zeros(1, 47);
for k = 1:45
    c = cov(fluctenv(:, k).', fluctenv(:, k+2).');
    denom = sqrt(c(1,1) * c(2,2));
    if denom > 0
        ki1(k) = c(2,1) / denom;
    end
end
% Boundary channels: pinned nearest-neighbor copies (docs/oracle-pins.md §2.4).
ki1(46) = ki1(45);
ki1(47) = ki1(45);
ki2 = zeros(1, 47);
ki2(1) = ki1(3);
ki2(2) = ki1(3);
ki2(3:47) = ki1(1:45);

% Specific fluctuation strength per channel.
gzi = Get_gzi_fluctuation(47);
kp = ki1 .* ki2;
fi = 0.4980 * (gzi .* mdepth_post.^1.7 .* abs(kp).^1.7 .* sign(kp));

dump = struct();
dump.x_provenance = sprintf(['SQAT @ 00b449e40599f1c1ef4abe0596094552213d57eb, ' ...
    'Octave %s, octave-signal %s, %s, scripts/dump_fs_stage.m'], ...
    OCTAVE_VERSION, ver('signal').Version, datestr(now));
dump.h0 = h0(:);
dump.mdepth_pre = mdepth_pre(:);
dump.mdepth_post = mdepth_post(:);
dump.ki1 = ki1(:);
dump.ki2 = ki2(:);
dump.fi = fi(:);
dump.frame_first10 = windowed(1:10)(:);
dump.frame_last10 = windowed(end-9:end)(:);
dump.post_a0_rms = post_a0_rms;

txt = jsonencode(dump);
% jsonencode cannot emit a leading-underscore field name; patch it in.
txt = strrep(txt, '"x_provenance"', '"_provenance"');
fid = fopen(out_path, 'w');
fprintf(fid, '%s', txt);
fclose(fid);
printf('wrote %s\n', out_path);
