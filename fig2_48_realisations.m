function summary = fig2_48_realisations(outputDir, baseSeed, plotOnly, plotRealisations)
%FIG2_48_REALISATIONS Run all four Figure 2 conditions 48 times with PARFOR.
%
% Usage (put this file and your model functions on the MATLAB path):
%   summary = fig2_48_realisations;
%   summary = fig2_48_realisations('fig2_results', 20260928);
% Replot saved results from the 84-second version without running simulations:
%   summary = fig2_48_realisations('existing_results_folder', [], true);
% Choose the saved realisations for the three movement-trace rows:
%   summary = fig2_48_realisations('existing_results_folder', [], true, 'random');
%   summary = fig2_48_realisations('existing_results_folder', [], true, 7);
%   summary = fig2_48_realisations('existing_results_folder', [], true, [3 12 42]);
% The fourth argument is optional (default 1 for all three rows). A scalar
% selects the same run for every row; three integers select rows in order:
% Healthy, PD, PD severe. Values must be within 1:completed_realisations.
% 'random' draws independently for each row, so run numbers may repeat.
% It uses the client RNG; rng(123) before the call makes the choice repeatable.
% This option also works for new simulations (third argument false).
% The selected numbers are printed, shown in row titles, and saved in
% summary.plot_realisations. To repeat a random selection, pass this saved
% vector as the fourth argument. All aggregate panels still use all 48 runs.
% On the first update, plot-only mode reads all 48 realisation files to
% calculate single-programme latency (6-30 s). It saves the updated summary
% as summary_replotted.mat plus new summary/latency CSVs and FIG/PNG files.
% Original summary.mat and realisation files are not rewritten. Subsequent
% replot calls reuse summary_replotted.mat if it is newer than summary.mat.
% Frequency-weighted power is derived from the saved original power, so
% repeated replot calls never multiply an already weighted spectrum again.
%
% Each of the 48 independent realisations contains all four conditions:
% 1 Healthy, 2 PD, 3 PD mild adaptations, 4 PD severe (192 simulations).
% Severe retains the previous PD-adaptation parameters. Mild retains 40%
% of each PD-to-severe parameter change: mild = PD + 0.4*(severe - PD).
% Plots show indices [1 2 4] only: Healthy, PD and PD severe. All four
% conditions remain in the simulations, saved arrays, statistics and CSVs.
% Completed realisations are saved immediately, in separate MAT files:
%   realisation_001.mat ... realisation_048.mat  (variable: runData)
%   manifest.mat  (configuration, inputs, and the planned file list)
%   summary.mat   (run-level arrays, mean, STD, SE, and 95% CIs)
%   figure2_48_realisations.fig / .png
%   metrics_by_state.fig / .png (all six metrics in all four states)
%   metrics_per_realisation.csv / metrics_summary.csv
%   metrics_all_movement_per_realisation.csv / metrics_all_movement_summary.csv
%   latency_single_per_realisation.csv / latency_single_summary.csv
%
% Load one realisation without loading the other 47:
%   S = load(fullfile(summary.output_dir, 'realisation_001.mat'), 'runData');
%   pd = S.runData.conditions{2};
%   plot(S.runData.inputs.ttime, pd.muscleposition);
%   singleError = squeeze(summary.metrics.error(:, 2, :)); % 4-by-48
%   twoProgrammeError = squeeze(summary.metrics.error(:, 3, :)); % 4-by-48
%   allMovementError = summary.movement_all.metrics.error; % 4-by-48
%   rigidity = summary.metrics.rigidity_covariance(2, 1, 7); % PD/rest/run 7
%   rigidityMean = summary.stats.rigidity_covariance.mean;  % 4-by-4
%   rigiditySE = summary.stats.rigidity_covariance.se;
%   rigiditySTD = summary.stats.rigidity_covariance.std;
%   lowerCI = summary.stats.rigidity_covariance.ci95_lower;
%   upperCI = summary.stats.rigidity_covariance.ci95_upper;
%   f = summary.joint_spectra.frequency_hz;
%   meanPDrest = squeeze(summary.joint_spectra.stats.power.mean(2, 1, :));
%   plot(f, meanPDrest);
%   weightedPDrest = squeeze(summary.joint_spectra.stats.frequency_weighted_power.mean(2, 1, :));
%   plot(f, weightedPDrest); % f .* original power, as shown in Figure 2
%
% Array order:
%   summary.metrics.<name>(condition, state, realisation) = 4-by-4-by-48.
%   summary.stats.<name>.<statistic>(condition, state) = 4-by-4.
%   summary.joint_spectra.<quantity>(condition,state,realisation,frequency).
%   summary.joint_spectra.stats.<quantity>.<statistic>(condition,state,frequency).
% Conditions: Healthy, PD, PD mild adaptations, PD severe.
% States: 1 rest (0-6 s), 2 movement_single (6-30 s),
% 3 movement_two (30-78 s), 4 hyperdirect (78-84 s).
% Whole-movement (6-78 s) metrics and spectra are also saved under
% summary.movement_all, separately from the four disjoint states.
% Metrics: error, velocity, tracking_covariance,
% tracking_correlation, rigidity_covariance, rigidity_correlation.
% Correlations are Pearson r. Constant signals give NaN, not zero. The
% statistics report n_valid, n_missing, and n_expected, and omit undefined
% values. All error bars/bands are mean +/- 1.96*SE (pointwise approximate
% 95% normal intervals for the mean across runs); correlation
% values are not Fisher transformed. STD is the sample SD across runs.
% SE = STD/sqrt(n_valid), normally using 48 realisations. Correlation
% intervals are not clipped to [-1,1]. Spectra show frequency .* power on
% LINEAR axes, with one condition per panel and an independent y-axis range
% for every panel. SE and interval bounds are scaled by the same frequency.
% Original power is retained; frequency_weighted_power and its statistics
% are saved separately for every condition/state/run, also for whole movement.
% Rest and movement flexor-extensor correlations use separate bar panels
% with independent y-axis ranges; movement groups single and two programmes.
% The upper figure uses equal total widths for trajectories/velocity on the
% left and three spectral columns on the right. Both figures use 12-point
% text and major/minor x/y gridlines with a line width of 1 point.
% Figure 2 is 160 mm wide by 200 mm high (16 x 20 cm), including its saved
% FIG geometry and PNG print area. The separate metrics figure keeps its size.
% Error is mean(abs(position(idx) - active_target(idx))) in joint-position
% units, calculated separately within each state. Both signals use the same
% timestamp t=(i-1)/fs: position(i+1) belongs to the next time step and must
% not be compared with target(i). Programme 1 is selected before 54 s and
% programme 2 from 54 s. At rest/hyperdirect, the scheduled moving target is
% still the reference, so error is not a measure of drift from a fixed pose.
%
% Requirements: Parallel Computing Toolbox, Wavelet Toolbox, Signal
% Processing Toolbox, and the five model functions checked below.
% Use a LOCAL PROCESS pool. The model functions must use the supplied state
% arguments and MATLAB's current random stream, without resetting RNG or
% retaining unreset hidden simulation state between calls.
%
% Healthy, PD and severe retain their condition-specific parameters.
% Mild is interpolated as above. cerebellum_CM_weight = [0.025 0.025 0.035 0.05].
% 21000 steps / 250 Hz = 84 s. Movement is extended from 48 to 72 s;
% rest and hyperdirect activation remain 6 s each. SMA drive is 0.4 for
% the selected programme and 0.15 for the unselected enabled programme;
% resting drive remains 0.05. Programme 2 is completely disabled from
% 6 to 30 s: no model calls, muscle output or CM history. Shared Golgi/L5
% state updates continue only through programme 1. At 30 s programme 2
% is enabled with empty CM history; the selected target switches at 54 s.
% Sinusoidal target phases stay continuous on the 84-second time axis.
%
% Changes needed for independent repeats and data-based plotting:
%   * Reset Golgi, CM, L5 and joint states for every condition and realisation.
%   * Use a separate reproducible Threefry substream per run/condition.
%     Healthy, PD and severe retain their assigned substream IDs; mild
%     has 48 additional IDs. The changed schedule changes the trajectories.
%   * Save raw outputs before smoothing; no plotting inside PARFOR.
%   * Use non-overlapping rest/single/two-programme/hyperdirect windows.
%   * Calculate six scalar metrics in EACH complete state of EACH run.
%     Error is mean absolute target-position error. Velocity retains the
%     phase-wise 10-sample moving mean of motor output, then mean(abs)*fs.
%     Tracking uses the scheduled active target even when movement is
%     inhibited. Covariance uses sample normalisation (N-1). Rigidity uses
%     the UNSMOOTHED flexor/extensor outputs, each summed over programmes.
%     These are the original gain-scaled muscle contributions, not forces
%     in physical units. Each correlation/covariance is one calculation
%     over the whole state, not a mean of sliding-window coefficients.
%   * The plotted latency is the maximum mean-centred cross-correlation lag
%     from the single-programme window (6-30 s), saved under latency_single.
%     Positive latency means achieved position lags the target; constant
%     signals give NaN. Periodic targets can give ambiguous correlation lags.
%     The older two-programme latency remains available under legacy only.
%   * Remove the artificial +0.001 velocity offset and fixed latency bars.
%   * Add Morse CWT of JOINT POSITION over the full record, then average
%     abs(CWT) and abs(CWT).^2 separately within each complete state.
%     Average the per-run distributions across runs, never the raw signals.
%     Gamma=3, time-bandwidth=60, voices/octave=10 preserve MATLAB defaults.
%     In beta/gamma notation this is beta=20, not beta=60.
%     The frequency grid is common to all conditions/states/realisations.
%     Full-state distributions include boundary/transition coefficients;
%     an additional power_coi distribution excludes recording-edge effects
%     using the cone of influence. This does not remove state-transition
%     mixing. Low-frequency estimates from six-second states need caution.
%     Joint CWT power is not a PSD in units per Hz or a probability density.
%     Frequency weighting is pointwise f .* power, not a sum over frequencies
%     or a normalisation. Weighted power has units position^2 Hz.
%   * Retain prior scalar metrics and thalamic spectra under summary.legacy.
%
% Reproducibility and saving pattern follow MathWorks guidance:
% https://www.mathworks.com/help/parallel-computing/repeat-random-numbers-in-parfor-loops.html
% https://www.mathworks.com/help/parallel-computing/save-variables-in-parfor-loop.html
% https://www.mathworks.com/help/wavelet/ref/cwt.html

if nargin < 3 || isempty(plotOnly)
    plotOnly = false;
end
validateattributes(plotOnly, {'logical'}, {'scalar'});
if nargin < 4 || isempty(plotRealisations)
    plotRealisations = 1;
end
plotConditionIndices = [1 2 4];
if plotOnly
    if nargin < 1 || isempty(outputDir)
        error('Fig2:MissingResultsFolder', 'Provide the folder containing the saved 84-second results.');
    end
    outputDir = char(outputDir);
    originalSummary = fullfile(outputDir, 'summary.mat');
    updatedSummary = fullfile(outputDir, 'summary_replotted.mat');
    sourceSummary = originalSummary;
    if exist(updatedSummary, 'file')
        originalInfo = dir(originalSummary);
        updatedInfo = dir(updatedSummary);
        if isempty(originalInfo) || updatedInfo.datenum >= originalInfo.datenum
            sourceSummary = updatedSummary;
        end
    end
    saved = load(sourceSummary, 'summary');
    summary = saved.summary;
    expectedStates = {'rest', 'movement_single', 'movement_two', 'hyperdirect'};
    if ~isequal(summary.state_names, expectedStates) || numel(summary.condition_labels) ~= 4
        error('Fig2:IncompatibleSavedResults', ...
            'Plot-only mode requires the four-condition, four-state (84-second) results.');
    end
    summary.config.plot_condition_indices = plotConditionIndices;
    [summary.config.plot_realisations, summary.config.plot_realisation_selection] = ...
        resolve_plot_realisations(plotRealisations, summary.completed_realisations, numel(plotConditionIndices));
    summary.plot_condition_indices = plotConditionIndices;
    summary.plot_condition_labels = summary.condition_labels(plotConditionIndices);
    summary.output_dir = outputDir;
    summary = ensure_single_programme_latency(summary, outputDir);
    summary = refresh_summary_intervals(summary);
    summary = add_frequency_weighted_power(summary);
    summary = record_plot_realisations(summary);
    save(updatedSummary, 'summary', '-v7.3');
    writetable(summary.metric_statistics_table, fullfile(outputDir, 'metrics_summary_replotted.csv'));
    writetable(summary.movement_all.metric_statistics_table, ...
        fullfile(outputDir, 'metrics_all_movement_summary_replotted.csv'));
    save_latency_tables(summary, outputDir);
    save_result_figures(summary, outputDir);
    fprintf('Figures updated using saved results; no simulations were run.\n');
    return;
end

if nargin < 1 || isempty(outputDir)
    outputDir = fullfile(pwd, ['fig2_results_' datestr(now, 'yyyymmdd_HHMMSS')]);
end
if nargin < 2 || isempty(baseSeed)
    baseSeed = 20260928;
end
validateattributes(baseSeed, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', '<=', 2^32-1});
outputDir = char(outputDir);

%% Settings
cfg.n_realisations = 48;
cfg.circuits = 500;
cfg.fs = 250;
cfg.duration_s.rest = 6;
cfg.duration_s.movement_single = 24;
cfg.duration_s.movement_two = 48;
cfg.duration_s.hyperdirect = 6;
cfg.total_duration_s = cfg.duration_s.rest + cfg.duration_s.movement_single ...
    + cfg.duration_s.movement_two + cfg.duration_s.hyperdirect;
cfg.repetitions = round(cfg.total_duration_s * cfg.fs);
cfg.sma_rest = 0.05;
cfg.sma_active = 0.4;
cfg.sma_inactive = 0.15;
cfg.base_seed = double(baseSeed);
cfg.random_generator = 'Threefry';
cfg.condition_labels = {'Healthy', 'PD', 'PD mild adaptations', 'PD severe'};
cfg.plot_condition_indices = plotConditionIndices;
[cfg.plot_realisations, cfg.plot_realisation_selection] = ...
    resolve_plot_realisations(plotRealisations, cfg.n_realisations, numel(plotConditionIndices));
cfg.condition_short_labels = {'Healthy', 'PD', 'PD mild', 'PD severe'};
cfg.condition_colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
    0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250];
cfg.mild_adaptation_fraction = 0.4;
cfg.state_names = {'rest', 'movement_single', 'movement_two', 'hyperdirect'};
cfg.state_labels = {'Rest', 'Single programme', 'Two programmes', 'Hyperdirect'};
cfg.metric_names = {'error', 'velocity', 'tracking_covariance', ...
    'tracking_correlation', 'rigidity_covariance', 'rigidity_correlation'};
cfg.schema_version = 8;
cfg.ci_multiplier = 1.96;
cfg.neuron_fire = 0.80;
cfg.resistance = 0;
cfg.mod_inhibition = -0.5;
cfg.mod_excitation = 0.5;
cfg.muscleposition_start = 0.5;
cfg.gain = 1.5 / cfg.circuits;
cfg.cerebellar_circuits = 100;
cfg.cm_buffer_length = fix(0.200 * cfg.fs);
cfg.smoothing_samples = 10;
cfg.target_amplitude = 0.15;
cfg.target_period = 10;
cfg.target_centre = 0.5;
cfg.make_figure = true;
cfg.wavelet.name = 'morse';
cfg.wavelet.gamma = 3;
cfg.wavelet.time_bandwidth = 60;
cfg.wavelet.beta = cfg.wavelet.time_bandwidth / cfg.wavelet.gamma;
cfg.wavelet.voices_per_octave = 10;
% Raw time series and frequency summaries are always saved. These optional
% arrays substantially increase disk and worker-memory requirements:
cfg.save_full_cwt = false;
cfg.save_neuron_thalamus = false;

% Original profiles are [Healthy, PD, severe]. The loop below inserts mild
% as the third entry, keeping all three original profiles unchanged.
cfg.cerebellum_VL_weight = [0.0125 0.0125 0.04];
cfg.cerebellum_CM_weight = [0.025 0.025 0.05];
cfg.stn_bg_connectivity = [0.25 0.25 2];
cfg.bg_cer_connectivity = [1.8 1.8 1.2];
cfg.CM_striatum_w = [0.004 0.004 0.24];
cfg.CM_acth_modifier = [0.5 0.5 1];
cfg.ach_levels = [2 1.5 1.5];
cfg.dopamine_released_fig = [2 0.8 0.8];
cfg.L5_CM_same_weight = [0.5 0.5 0.5];
cfg.L5_CM_opp_weight = [0.2 0.2 0.2];
cfg.stn_fire_rate = [0.22 0.22 0.3];
cfg.cortex_acth_modifier = [0.3 0.3 0.75];
cfg.BG_CM_weight = [1.2 1.2 1];
cfg.BG_VL_weight = [3 3 2.5];
conditionParameters = {'cerebellum_VL_weight', 'cerebellum_CM_weight', ...
    'stn_bg_connectivity', 'bg_cer_connectivity', 'CM_striatum_w', ...
    'CM_acth_modifier', 'ach_levels', 'dopamine_released_fig', ...
    'L5_CM_same_weight', 'L5_CM_opp_weight', 'stn_fire_rate', ...
    'cortex_acth_modifier', 'BG_CM_weight', 'BG_VL_weight'};
for parameterIndex = 1:numel(conditionParameters)
    name = conditionParameters{parameterIndex};
    profile = cfg.(name);
    mild = profile(2) + cfg.mild_adaptation_fraction * (profile(3) - profile(2));
    cfg.(name) = [profile(1:2), mild, profile(3)];
end

% Keep each condition/realisation's previously assigned stream ID. The new
% duration, drive and disabling schedule intentionally change trajectories.
% All 192 condition/realisations still have distinct substreams.
originalSubstreams = reshape(1:(3*cfg.n_realisations), 3, []).';
mildSubstreams = (3*cfg.n_realisations + (1:cfg.n_realisations)).';
cfg.random_substreams = [originalSubstreams(:, 1:2), mildSubstreams, originalSubstreams(:, 3)];
cfg.matlab_version = version;
cfg.matlab_release = version('-release');
cfg.created_at = datestr(now, 30);

%% Check dependencies and destination before starting the simulations
modelFunctions = {'basal_ganglia_fcn_pts_rtl', 'hbm_vnr', ...
    'M1_fcn_pts_rtl', 'spinal_circuits_fcn_pts', 'put_in_range'};
requiredFunctions = [modelFunctions, {'parpool', 'cwt', 'xcorr'}];
for k = 1:numel(requiredFunctions)
    if isempty(which(requiredFunctions{k}))
        error('Fig2:MissingDependency', ...
            'Required function %s is not on the MATLAB path.', requiredFunctions{k});
    end
end
cfg.model_function_paths = cellfun(@which, modelFunctions, 'UniformOutput', false);
if exist(outputDir, 'file') && ~exist(outputDir, 'dir')
    error('Fig2:OutputPath', 'The output path is an existing file.');
end
if exist(outputDir, 'dir')
    listing = dir(outputDir);
    listing = listing(~ismember({listing.name}, {'.', '..'}));
    if ~isempty(listing)
        error('Fig2:OutputExists', ...
            'Choose a new or empty output folder. Existing results will not be overwritten.');
    end
else
    [ok, msg] = mkdir(outputDir);
    if ~ok
        error('Fig2:OutputFolder', '%s', msg);
    end
end
[ok, attributes] = fileattrib(outputDir);
if ~ok
    error('Fig2:OutputFolder', 'Could not resolve the output folder.');
end
outputDir = attributes.Name; % Absolute path is shared with local workers.

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local');
end
if isa(pool, 'parallel.ThreadPool')
    error('Fig2:ProcessPoolRequired', ...
        'Use a local process pool: delete(gcp(''nocreate'')); parpool(''local'');');
end
if ~isa(pool.Cluster, 'parallel.cluster.Local')
    error('Fig2:LocalPoolRequired', ...
        'This script saves to the local disk. Start a local process pool with parpool(''local'').');
end

inputs = make_inputs(cfg);
nRuns = cfg.n_realisations;
runFiles = cell(nRuns, 1);
for r = 1:nRuns
    runFiles{r} = fullfile(outputDir, sprintf('realisation_%03d.mat', r));
end
save(fullfile(outputDir, 'manifest.mat'), 'cfg', 'inputs', 'runFiles', '-v7.3');

%% Independent runs; each worker writes only its own file
smallResults = cell(nRuns, 1);
progress = parallel.pool.DataQueue;
afterEach(progress, @(r) fprintf('Saved realisation %02d of %02d.\n', r, nRuns));
fprintf('Running %d realisations, each containing %d conditions.\n', ...
    nRuns, numel(cfg.condition_labels));
fprintf('Output folder: %s\n', outputDir);

parfor realisation = 1:nRuns
    runData = simulate_realisation(realisation, cfg, inputs);
    save_realisation(runFiles{realisation}, runData);
    % Return only compact summaries, not all raw data, to the client.
    smallResults{realisation} = compact_results(runData);
    send(progress, realisation);
end

%% Aggregate after all workers finish; raw files already exist independently
summary = aggregate_results(smallResults, cfg, runFiles, outputDir);
summary = add_frequency_weighted_power(summary);
summary = record_plot_realisations(summary);
save(fullfile(outputDir, 'summary.mat'), 'summary', '-v7.3');
try
    writetable(summary.metric_values_table, fullfile(outputDir, 'metrics_per_realisation.csv'));
    writetable(summary.metric_statistics_table, fullfile(outputDir, 'metrics_summary.csv'));
    writetable(summary.movement_all.metric_values_table, ...
        fullfile(outputDir, 'metrics_all_movement_per_realisation.csv'));
    writetable(summary.movement_all.metric_statistics_table, ...
        fullfile(outputDir, 'metrics_all_movement_summary.csv'));
    save_latency_tables(summary, outputDir);
catch exportError
    warning('Fig2:CSVExportFailed', ...
        'All data and tables are saved in MAT files; CSV export failed: %s', exportError.message);
end
fprintf('All %d realisations saved. Summary saved to summary.mat.\n', nRuns);

if cfg.make_figure
    % Plotting failure must not prevent delivery of completed simulations.
    try
        save_result_figures(summary, outputDir);
    catch plotError
        warning('Fig2:PlotFailed', ...
            'Simulation data are saved, but plotting failed: %s', plotError.message);
    end
end
end

function inputs = make_inputs(cfg)
N = cfg.repetitions;
n = N + 1;
durations = cellfun(@(name) cfg.duration_s.(name), cfg.state_names);
samples = durations * cfg.fs;
assert(all(durations > 0) && all(abs(samples-round(samples)) < 1e-9), ...
    'Each state must contain a positive integer number of samples.');
samples = round(samples);
edges = [0 cumsum(samples)];
assert(edges(end) == N, 'State durations must sum to the simulation duration.');
assert(mod(samples(3), 2) == 0, 'Two-programme movement must split into equal sample counts.');
movementStart = edges(2) + 1;
switchSample = edges(3) + samples(3)/2;
inputs.ttime = (0:N) / cfg.fs;
inputs.sample_time = inputs.ttime(1:N);
f = 1 / cfg.target_period;
inputs.desiredposition = [ ...
    cfg.target_centre + cfg.target_amplitude * sin(2*pi*f*inputs.ttime); ...
    cfg.target_centre + 1.5*cfg.target_amplitude * sin(2*pi*f*inputs.ttime + pi/2)];
% Explicit times keep rest and hyperdirect at 6 s when movement is extended.
% Intervals are [start,end); the extra endpoint is only for plotting.
for s = 1:numel(cfg.state_names)
    inputs.phase_indices.(cfg.state_names{s}) = edges(s)+1:edges(s+1);
end
inputs.phase_indices.movement = edges(2)+1:edges(4); % Full 6-78 s summary.
inputs.phase_times_s = [edges(1:end-1).', edges(2:end).'] / cfg.fs;
inputs.movement_times_s = [edges(2) edges(4)] / cfg.fs;
inputs.two_programme_enable_s = edges(3) / cfg.fs;
inputs.target_switch_s = switchSample / cfg.fs;
inputs.is_rest = false(1, N);
inputs.is_rest(inputs.phase_indices.rest) = true;

inputs.program_enabled = true(2, n);
inputs.program_enabled(2, inputs.phase_indices.movement_single) = false;
inputs.programme_disable_note = ['Programme 2 is not simulated during single-programme movement; ' ...
    'all of its outputs are zero and its programme-specific CM buffers are cleared.'];
inputs.motor_drive = cfg.sma_rest * ones(2, n);
inputs.motor_drive(:, movementStart:n) = cfg.sma_inactive;
inputs.motor_drive(1, movementStart:switchSample) = cfg.sma_active;
inputs.motor_drive(2, switchSample+1:n) = cfg.sma_active;
inputs.motor_drive(~inputs.program_enabled) = 0;
inputs.hyperdirect_firing = zeros(1, n);
inputs.hyperdirect_firing(edges(4)+1:n) = 1;
inputs.dopamine_initial = 0.75 * ones(1, n);
inputs.dopamine_change_idx = movementStart:n;
inputs.motorgroupselected = ones(1, N);
inputs.motorgroupselected(switchSample+1:N) = 2;
activeIndices = sub2ind([2 n], inputs.motorgroupselected, 1:N);
inputs.active_target = inputs.desiredposition(activeIndices);
% Shift the original legacy windows by the inserted 24 s. Latency uses
% 30-54 s (programme 1 selected with both programmes present). The narrow
% legacy thalamic spectrum uses 36-42 s. New state metrics use full states.
inputs.tracking_indices = edges(3)+1:switchSample;
inputs.tracking_times_s = [edges(3) switchSample] / cfg.fs;
inputs.spectrum_indices = inputs.phase_indices;
inputs.spectrum_indices.movement = edges(3)+round(6*cfg.fs)+1:edges(3)+round(12*cfg.fs);
end

function runData = simulate_realisation(realisation, cfg, inputs)
runData.realisation = realisation;
runData.config = cfg;
runData.inputs = inputs;
runData.conditions = cell(1, numel(cfg.condition_labels));
runData.started_at = datestr(now, 30);
timer = tic;
for condition = 1:numel(cfg.condition_labels)
    runData.conditions{condition} = simulate_condition(realisation, condition, cfg, inputs);
end
% One scalar for every metric, condition and state in this realisation.
for m = 1:numel(cfg.metric_names)
    name = cfg.metric_names{m};
    values = nan(numel(cfg.condition_labels), numel(cfg.state_names));
    allMovement = nan(numel(cfg.condition_labels), 1);
    for condition = 1:numel(cfg.condition_labels)
        values(condition, :) = runData.conditions{condition}.state_metrics.(name);
        allMovement(condition) = runData.conditions{condition}.movement_all_metrics.(name);
    end
    runData.metrics.(name) = values;
    runData.movement_all_metrics.(name) = allMovement;
end
runData.metrics_layout = 'condition x state; one value per complete state';
runData.elapsed_seconds = toc(timer);
runData.completed_at = datestr(now, 30);
end

function data = simulate_condition(realisation, condition, cfg, inputs)
% All evolving states are local to one condition. The supplied script kept
% several of these states outside its condition loop, causing carry-over.
N = cfg.repetitions;
P = size(inputs.desiredposition, 1);
C = cfg.circuits;
substream = cfg.random_substreams(realisation, condition);
stream = RandStream(cfg.random_generator, 'Seed', cfg.base_seed);
stream.Substream = substream;
oldStream = RandStream.getGlobalStream;
cleanup = onCleanup(@() RandStream.setGlobalStream(oldStream)); %#ok<NASGU>
RandStream.setGlobalStream(stream);

data.condition_index = condition;
data.condition_label = cfg.condition_labels{condition};
data.rng.generator = cfg.random_generator;
data.rng.seed = cfg.base_seed;
data.rng.substream = substream;
data.rng.initial_state = stream.State;
data.dopamine_levels = inputs.dopamine_initial;
data.dopamine_levels(inputs.dopamine_change_idx) = cfg.dopamine_released_fig(condition);
data.ach_level = cfg.ach_levels(condition);

weight_params.keepdefaults = true;
weightNames = {'cerebellum_VL_weight', 'cerebellum_CM_weight', ...
    'CM_acth_modifier', 'CM_striatum_w', 'L5_CM_same_weight', ...
    'L5_CM_opp_weight', 'stn_fire_rate', 'cortex_acth_modifier', ...
    'BG_CM_weight', 'BG_VL_weight'};
for k = 1:numel(weightNames)
    values = cfg.(weightNames{k});
    weight_params.(weightNames{k}) = values(condition);
end
pathologies.keepdefaults = true;
data.weight_params = weight_params;
data.pathologies = pathologies;
data.stn_bg_connectivity = cfg.stn_bg_connectivity(condition);
data.bg_cer_connectivity = cfg.bg_cer_connectivity(condition);

golgi_cellfire_flexor = zeros(cfg.cerebellar_circuits, 1);
golgi_cellfire_extensor = zeros(cfg.cerebellar_circuits, 1);
cerebellar_agonist = zeros(cfg.cerebellar_circuits, 1);
cerebellar_antagonist = zeros(cfg.cerebellar_circuits, 1);
CM_nucleus.agonist = zeros(P, cfg.cm_buffer_length);
CM_nucleus.antagonist = zeros(P, cfg.cm_buffer_length);
CM_index = 1;
% L5 and Golgi state remain shared and are updated sequentially by ENABLED
% programmes. Disabling programme 2 must not reset programme 1's state.
L5_flexor_prev = zeros(1, C);
L5_extensor_prev = zeros(1, C);
data.muscleposition = cfg.muscleposition_start * ones(1, N+1);
data.output_raw = zeros(1, N);
data.correct_motor = false(1, N);
data.error_observed = zeros(1, N);
data.CM_agonist_keep = zeros(P, N);
data.CM_antagonist_keep = zeros(P, N);
data.gpi_agonist_keep = zeros(P, N);
data.gpi_antagonist_keep = zeros(P, N);
data.acth_rr_agonist = zeros(P, N);
data.acth_rr_antagonist = zeros(P, N);
data.L5_flexor_mean = zeros(P, N);
data.L5_extensor_mean = zeros(P, N);
movement.muscle_flexor = zeros(P, N);
movement.muscle_extensor = zeros(P, N);
movement.thalamus_flexor = zeros(P, N);
movement.thalamus_extensor = zeros(P, N);
if cfg.save_neuron_thalamus
    data.thalamus_flexor_by_program = zeros(C, N, P);
    data.thalamus_extensor_by_program = zeros(C, N, P);
end
output_mp = zeros(P, 1);

for i = 1:N
    % Clear every sample so disabling a programme cannot reuse its previous
    % muscle contribution. Its preallocated activity/output traces stay zero.
    output_mp(:) = 0;
    for motor_programs = 1:P
        if ~inputs.program_enabled(motor_programs, i)
            CM_nucleus.agonist(motor_programs, :) = 0;
            CM_nucleus.antagonist(motor_programs, :) = 0;
            continue; % Skip BG, cerebellar, cortical and spinal model calls.
        end
        bg_agonist = basal_ganglia_fcn_pts_rtl( ...
            inputs.motor_drive(motor_programs, i), cfg.ach_levels(condition), ...
            data.dopamine_levels(i), weight_params, pathologies, 10, ...
            inputs.hyperdirect_firing(i), CM_nucleus.agonist(motor_programs, :), CM_index);
        bg_antagonist = basal_ganglia_fcn_pts_rtl( ...
            inputs.motor_drive(motor_programs, i), cfg.ach_levels(condition), ...
            data.dopamine_levels(i), weight_params, pathologies, 10, ...
            inputs.hyperdirect_firing(i), CM_nucleus.antagonist(motor_programs, :), CM_index);
        data.acth_rr_agonist(motor_programs, i) = bg_agonist.acth_release_rate;
        data.acth_rr_antagonist(motor_programs, i) = bg_antagonist.acth_release_rate;

        desired_input = 0.04 + 0.32*inputs.desiredposition(motor_programs, i);
        muscleposition_input = 0.04 + 0.32*data.muscleposition(i);
        desired_input_rev = (0.04 + 0.32) - desired_input + 0.04;
        muscleposition_input_rev = (0.04 + 0.32) - muscleposition_input + 0.04;
        desired_input = put_in_range(desired_input ...
            - cfg.bg_cer_connectivity(condition)*rand*bg_agonist.gpi_dep_probability, 0, 0.36);
        muscleposition_input = put_in_range(muscleposition_input ...
            + cfg.stn_bg_connectivity(condition)*rand*bg_antagonist.stn_dep_probability ...
            - 0.075*cfg.bg_cer_connectivity(condition)*rand*bg_antagonist.gpi_dep_probability, 0, 0.36);
        desired_input_rev = put_in_range(desired_input_rev ...
            - cfg.bg_cer_connectivity(condition)*rand*bg_antagonist.gpi_dep_probability, 0, 0.36);
        muscleposition_input_rev = put_in_range(muscleposition_input_rev ...
            + cfg.stn_bg_connectivity(condition)*rand*bg_agonist.stn_dep_probability ...
            - 0.075*cfg.bg_cer_connectivity(condition)*rand*bg_agonist.gpi_dep_probability, 0, 0.36);
        if inputs.is_rest(i)
            desired_input = 0;
            desired_input_rev = 0;
        end
        for y = 1:cfg.cerebellar_circuits
            [cerebellar_agonist(y), golgi_cellfire_flexor(y), ~] = hbm_vnr( ...
                desired_input, muscleposition_input_rev, golgi_cellfire_flexor(y), ...
                cfg.mod_inhibition, cfg.mod_excitation, cfg.neuron_fire);
            [cerebellar_antagonist(y), golgi_cellfire_extensor(y), ~] = hbm_vnr( ...
                desired_input_rev, muscleposition_input, golgi_cellfire_extensor(y), ...
                cfg.mod_inhibition, cfg.mod_excitation, cfg.neuron_fire);
        end

        % Preserve the original CM buffer-index update order.
        CM_index = mod(i-1, cfg.cm_buffer_length) + 1;
        [L5_flexor, L5_extensor, thalamus_flexor, thalamus_extensor, CM_tmp] = ...
            M1_fcn_pts_rtl(cerebellar_agonist, cerebellar_antagonist, ...
            bg_agonist.gpi_firings, bg_antagonist.gpi_firings, ...
            inputs.motor_drive(motor_programs, i), C, ...
            L5_flexor_prev, L5_extensor_prev, weight_params);
        L5_flexor_prev = L5_flexor(:).';
        L5_extensor_prev = L5_extensor(:).';
        data.L5_flexor_mean(motor_programs, i) = mean(L5_flexor(:));
        data.L5_extensor_mean(motor_programs, i) = mean(L5_extensor(:));
        data.gpi_agonist_keep(motor_programs, i) = mean(bg_agonist.gpi_firings(:));
        data.gpi_antagonist_keep(motor_programs, i) = mean(bg_antagonist.gpi_firings(:));
        CM_nucleus.agonist(motor_programs, CM_index) = CM_tmp.CM_agonist;
        CM_nucleus.antagonist(motor_programs, CM_index) = CM_tmp.CM_antagonist;
        data.CM_agonist_keep(motor_programs, i) = CM_tmp.CM_agonist;
        data.CM_antagonist_keep(motor_programs, i) = CM_tmp.CM_antagonist;

        spinal = spinal_circuits_fcn_pts(L5_flexor, L5_extensor, ...
            data.muscleposition(i), weight_params);
        flexor_output = sum(spinal.flexorpower(:) * cfg.gain/10);
        extensor_output = sum(spinal.extensorpower(:) * cfg.gain/10);
        output_mp(motor_programs) = flexor_output - extensor_output;
        movement.muscle_flexor(motor_programs, i) = flexor_output;
        movement.muscle_extensor(motor_programs, i) = extensor_output;
        movement.thalamus_flexor(motor_programs, i) = sum(thalamus_flexor(:));
        movement.thalamus_extensor(motor_programs, i) = sum(thalamus_extensor(:));
        if cfg.save_neuron_thalamus
            data.thalamus_flexor_by_program(:, i, motor_programs) = thalamus_flexor(:);
            data.thalamus_extensor_by_program(:, i, motor_programs) = thalamus_extensor(:);
        end
    end
    data.output_raw(i) = sum(output_mp);
    % Both values are at t=(i-1)/fs. nextPosition below is at t+1/fs.
    % Keep this signed trace; take abs BEFORE averaging for reported MAE.
    data.error_observed(i) = data.muscleposition(i) - inputs.active_target(i);
    data.correct_motor(i) = (-data.error_observed(i) * data.output_raw(i)) > 0;
    nextPosition = data.muscleposition(i) + data.output_raw(i) - cfg.resistance;
    if ~isfinite(nextPosition)
        error('Fig2:NonfinitePosition', ...
            'Nonfinite joint position in realisation %d, condition %d, step %d.', ...
            realisation, condition, i);
    end
    if nextPosition < 0
        nextPosition = 0;
    elseif nextPosition > 1
        nextPosition = 1;
    end
    data.muscleposition(i+1) = nextPosition;
end

data.movement_covariance_check = movement;
data.target_pos_correlation = inputs.active_target;
data.position_correlation = data.muscleposition(1:N);
data.output = movmean(data.output_raw, cfg.smoothing_samples);
data.moving_index = movmean(double(data.output ~= 0), cfg.smoothing_samples);
% The supplied script overwrote thalamus_flexor/extensor on each programme
% call, so its post-loop activity_thalamus represented programme 2 only.
data.activity_thalamus = movement.thalamus_flexor(end, :) + movement.thalamus_extensor(end, :);
data.activity_thalamus_all_programs = sum(movement.thalamus_flexor + movement.thalamus_extensor, 1);
data.thalamus_spectrum_signal = 'activity_thalamus (last motor programme, as in source)';
phase = inputs.phase_indices;
data.power_rest = movmean(data.output_raw(phase.rest), cfg.smoothing_samples);
data.power_movement = movmean(data.output_raw(phase.movement), cfg.smoothing_samples);
data.hyperdirect_movement = movmean(data.output_raw(phase.hyperdirect), cfg.smoothing_samples);
data.muscle_flexor_total = sum(movement.muscle_flexor, 1);
data.muscle_extensor_total = sum(movement.muscle_extensor, 1);
data.state_metrics = calculate_state_metrics(data, inputs, cfg);
data.movement_all_metrics = calculate_state_metrics(data, inputs, cfg, {'movement'});
data.state_metric_samples = cellfun(@(name) numel(phase.(name)), cfg.state_names);
data.latency_single_s = calculate_single_programme_latency(data, inputs, cfg.fs);

metrics.velocity_rest = mean(abs(data.power_rest)) * cfg.fs;
metrics.velocity_movement = mean(abs(data.power_movement)) * cfg.fs;
metrics.velocity_hyperdirect = mean(abs(data.hyperdirect_movement)) * cfg.fs;
metrics.movement_error = mean(abs(data.error_observed(phase.movement)));
metrics.error_all_phases = mean(abs(data.error_observed));
metrics.correct_movement_fraction = mean(data.correct_motor(phase.movement));
metrics.zero_movement_fraction = mean(data.output_raw(phase.movement) == 0);
actual = data.muscleposition(inputs.tracking_indices);
target = inputs.desiredposition(1, inputs.tracking_indices);
covariance = cov(actual, target);
metrics.tracking_covariance = covariance(1, 2);
[crossCovariance, lags] = xcorr(actual - mean(actual), target - mean(target));
[~, peakIndex] = max(crossCovariance);
metrics.latency_s = lags(peakIndex) / cfg.fs;
data.lag_keep = lags(peakIndex);
[originalCorrelation, originalLags] = xcorr(actual, target);
[~, originalPeak] = max(originalCorrelation);
data.legacy_uncentred_lag_samples = originalLags(originalPeak);
data.latency_method = 'Maximum mean-centred xcorr; positive means achieved position lags target';
data.latency_note = 'Periodic targets can yield ambiguous lags; this is not a circular-phase analysis.';
data.latency_times_s = inputs.tracking_times_s;
data.metrics = metrics;

% Retain the original default CWT call and amplitude-distribution definition.
% No boundary exclusion is silently introduced; save COI for later analysis.
[coefficients, frequency, coi] = cwt(data.activity_thalamus, cfg.fs);
magnitude = abs(coefficients);
data.power_condition.freq = frequency;
data.power_condition.coi = coi;
data.power_condition.method = 'cwt(activity_thalamus, fs), MATLAB default Morse wavelet';
data.power_condition.quantity = 'Mean absolute CWT coefficient (magnitude, not squared power)';
data.power_condition.boundary_exclusion_applied = false;
names = {'rest', 'movement', 'hyperdirect'};
for k = 1:numel(names)
    name = names{k};
    idx = inputs.spectrum_indices.(name);
    data.power_condition.([name 'distribution']) = mean(magnitude(:, idx), 2);
    data.power_condition.mean_squared_magnitude.(name) = mean(magnitude(:, idx).^2, 2);
    data.power_condition.analysis_times_s.(name) = inputs.sample_time(idx([1 end]));
end
if cfg.save_full_cwt
    data.power_condition.coefficients = coefficients;
end
% Joint distributions use all four complete states; the narrow 36-42 s
% legacy thalamic window is not used for joint-position analysis.
data.joint_spectra = calculate_joint_spectra(data.muscleposition(1:N), inputs, cfg);
data.final_state.CM_nucleus = CM_nucleus;
data.final_state.CM_index = CM_index;
data.final_state.golgi_cellfire_flexor = golgi_cellfire_flexor;
data.final_state.golgi_cellfire_extensor = golgi_cellfire_extensor;
data.final_state.L5_flexor_prev = L5_flexor_prev;
data.final_state.L5_extensor_prev = L5_extensor_prev;
data.rng.final_state = stream.State;
end

function metrics = calculate_state_metrics(data, inputs, cfg, stateNames)
if nargin < 4
    stateNames = cfg.state_names;
end
S = numel(stateNames);
for m = 1:numel(cfg.metric_names)
    metrics.(cfg.metric_names{m}) = nan(1, S);
end
for s = 1:S
    idx = inputs.phase_indices.(stateNames{s});
    achieved = data.muscleposition(idx);
    target = inputs.active_target(idx);
    % Same-time position and selected target; no smoothing, lag alignment,
    % amplitude normalisation or averaging over inactive programmes.
    metrics.error(s) = mean(abs(achieved - target));
    % Keep the existing velocity definition and its within-state smoothing.
    smoothedOutput = movmean(data.output_raw(idx), cfg.smoothing_samples);
    metrics.velocity(s) = mean(abs(smoothedOutput)) * cfg.fs;
    [metrics.tracking_covariance(s), metrics.tracking_correlation(s)] = ...
        covariance_and_correlation(target, achieved);
    [metrics.rigidity_covariance(s), metrics.rigidity_correlation(s)] = ...
        covariance_and_correlation(data.muscle_flexor_total(idx), ...
        data.muscle_extensor_total(idx));
end
end

function latencySeconds = calculate_single_programme_latency(data, inputs, fs)
% Retain the existing lag method, restricted to single-programme movement.
% Same-time samples are used; positive lag means position follows target.
idx = inputs.phase_indices.movement_single;
actual = data.muscleposition(idx);
target = inputs.active_target(idx);
actual = actual(:); target = target(:);
latencySeconds = NaN;
if numel(actual) < 2 || numel(actual) ~= numel(target) ...
        || any(~isfinite(actual)) || any(~isfinite(target)) ...
        || all(actual == actual(1)) || all(target == target(1))
    return;
end
[crossCovariance, lags] = xcorr(actual-mean(actual), target-mean(target));
[~, peakIndex] = max(crossCovariance);
latencySeconds = lags(peakIndex) / fs;
end

function [sampleCovariance, pearsonR] = covariance_and_correlation(x, y)
% Sample covariance and Pearson r over all paired samples in one state.
x = x(:); y = y(:);
sampleCovariance = NaN; pearsonR = NaN;
if numel(x) ~= numel(y) || numel(x) < 2 || any(~isfinite(x)) || any(~isfinite(y))
    return;
end
if all(x == x(1)) || all(y == y(1))
    sampleCovariance = 0;
    return; % Correlation of a constant signal is undefined.
end
xc = x - mean(x); yc = y - mean(y);
crossProduct = sum(xc .* yc);
sampleCovariance = crossProduct / (numel(x) - 1);
denominator = sqrt(sum(xc.^2)) * sqrt(sum(yc.^2));
if denominator > 0
    pearsonR = max(-1, min(1, crossProduct / denominator));
end
end

function spectra = calculate_joint_spectra(position, inputs, cfg)
% Transform the full record once to retain a common grid across the states.
[coefficients, frequency, coi] = cwt(position, cfg.wavelet.name, cfg.fs, ...
    'WaveletParameters', [cfg.wavelet.gamma cfg.wavelet.time_bandwidth], ...
    'VoicesPerOctave', cfg.wavelet.voices_per_octave);
[frequency, order] = sort(frequency(:), 'ascend');
coefficients = coefficients(order, :);
magnitude = abs(coefficients);
S = numel(cfg.state_names);
F = numel(frequency);
analysisNames = [cfg.state_names, {'movement'}];
A = numel(analysisNames);
spectra.frequency_hz = frequency.';
spectra.state_names = cfg.state_names;
spectra.state_times_s = inputs.phase_times_s;
spectra.layout = 'state x frequency';
spectra.magnitude = nan(A, F);
spectra.power = nan(A, F);
spectra.power_coi = nan(A, F);
spectra.coi_valid_time_samples = zeros(A, F);
spectra.state_sample_count = zeros(1, A);
spectra.sample_times_s = zeros(A, 2);
spectra.coi_hz = coi(:).';
spectra.signal = 'muscleposition(1:repetitions), no smoothing or detrending';
spectra.wavelet = cfg.wavelet;
spectra.magnitude_definition = 'Mean over state time samples of abs(CWT)';
spectra.power_definition = 'Mean over state time samples of abs(CWT).^2; not PSD per Hz';
spectra.boundary_note = ['magnitude and power use all state samples. power_coi uses only ' ...
    'frequencies above the recording cone-of-influence boundary. This extra mask ' ...
    'does not remove wavelet mixing across internal state transitions.'];
for s = 1:A
    idx = inputs.phase_indices.(analysisNames{s});
    stateMagnitude = magnitude(:, idx);
    statePower = stateMagnitude.^2;
    spectra.magnitude(s, :) = mean(stateMagnitude, 2).';
    spectra.power(s, :) = mean(statePower, 2).';
    % COI returned with fs is in Hz. At a given time, lower frequencies
    % than this boundary are affected by the edges of the full recording.
    valid = bsxfun(@gt, frequency, spectra.coi_hz(idx));
    counts = sum(valid, 2);
    statePower(~valid) = 0;
    edgeExcluded = sum(statePower, 2) ./ max(counts, 1);
    edgeExcluded(counts == 0) = NaN;
    spectra.power_coi(s, :) = edgeExcluded.';
    spectra.coi_valid_time_samples(s, :) = counts.';
    spectra.state_sample_count(s) = numel(idx);
    spectra.sample_times_s(s, :) = inputs.sample_time(idx([1 end]));
end
% Keep the overlapping whole-movement summary outside the four-state axis.
rowFields = {'magnitude', 'power', 'power_coi', 'coi_valid_time_samples', 'sample_times_s'};
for k = 1:numel(rowFields)
    name = rowFields{k};
    values = spectra.(name);
    spectra.movement_all.(name) = values(A, :);
    spectra.(name) = values(1:S, :);
end
spectra.movement_all.state_sample_count = spectra.state_sample_count(A);
spectra.movement_all.interval_s = inputs.movement_times_s;
spectra.state_sample_count = spectra.state_sample_count(1:S);
if cfg.save_full_cwt
    spectra.coefficients = coefficients;
end
end

function save_realisation(filename, runData)
% SAVE is deliberately inside a helper, not directly in the PARFOR body.
% A final filename appears only after its complete MAT file has been written.
temporaryFile = [filename '.partial.mat'];
save(temporaryFile, 'runData', '-v7.3');
[ok, msg] = movefile(temporaryFile, filename);
if ~ok
    error('Fig2:SaveFailed', 'Could not finalise %s: %s', filename, msg);
end
end

function compact = compact_results(runData)
C = numel(runData.conditions);
compact.metrics = runData.metrics;
compact.movement_all_metrics = runData.movement_all_metrics;
compact.state_times_s = runData.inputs.phase_times_s;
compact.movement_times_s = runData.inputs.movement_times_s;
compact.tracking_times_s = runData.inputs.tracking_times_s;
singleState = find(strcmp(runData.config.state_names, 'movement_single'));
compact.latency_single_window_s = runData.inputs.phase_times_s(singleState, :);
compact.latency_single_s = nan(C, 1);
compact.legacy_metrics = cell(1, C);
compact.legacy_spectra = cell(1, C);
compact.joint_spectra = cell(1, C);
for c = 1:C
    compact.latency_single_s(c) = runData.conditions{c}.latency_single_s;
    compact.legacy_metrics{c} = runData.conditions{c}.metrics;
    pc = runData.conditions{c}.power_condition;
    compact.legacy_spectra{c}.freq = pc.freq;
    compact.legacy_spectra{c}.restdistribution = pc.restdistribution;
    compact.legacy_spectra{c}.movementdistribution = pc.movementdistribution;
    compact.legacy_spectra{c}.hyperdirectdistribution = pc.hyperdirectdistribution;
    js = runData.conditions{c}.joint_spectra;
    % Do not send full wavelet coefficients or time series back to client.
    compact.joint_spectra{c}.frequency_hz = js.frequency_hz;
    compact.joint_spectra{c}.magnitude = js.magnitude;
    compact.joint_spectra{c}.power = js.power;
    compact.joint_spectra{c}.power_coi = js.power_coi;
    compact.joint_spectra{c}.coi_valid_time_samples = js.coi_valid_time_samples;
    compact.joint_spectra{c}.state_sample_count = js.state_sample_count;
    compact.joint_spectra{c}.movement_all = js.movement_all;
end
end

function summary = aggregate_results(results, cfg, runFiles, outputDir)
R = numel(results);
C = numel(cfg.condition_labels);
S = numel(cfg.state_names);
summary.config = cfg;
summary.output_dir = outputDir;
summary.files = runFiles;
summary.condition_labels = cfg.condition_labels;
summary.plot_condition_indices = cfg.plot_condition_indices;
summary.plot_condition_labels = cfg.condition_labels(cfg.plot_condition_indices);
summary.state_names = cfg.state_names;
summary.state_labels = cfg.state_labels;
summary.state_times_s = results{1}.state_times_s;
summary.metric_names = cfg.metric_names;
summary.completed_realisations = R;
summary.metrics_layout = sprintf('condition x state x realisation (%d x %d x %d)', C, S, R);
summary.statistics_layout = 'condition x state';
summary.interval_method = ['Arithmetic mean, sample STD, SE=STD/sqrt(n_valid), ' ...
    'and pointwise mean +/- 1.96*SE (approximate 95% normal interval)'];
summary.undefined_correlation_policy = ['Constant signals give NaN. Aggregation omits ' ...
    'undefined values and reports n_valid and n_missing; SE uses n_valid.'];
summary.correlation_averaging = 'Arithmetic mean of per-realisation Pearson r values; no Fisher transform';
summary.comparison_note = 'Conditions use independent substreams; no paired-condition design or significance tests are implied.';
summary.figure_note = ['Trajectories show the realisations identified in their row titles; joint spectra and six state metrics ' ...
    'summarise each of the four states separately. Plotted latency uses single-programme ' ...
    'movement (6-30 s). Plots show condition indices 1, 2 and 4, with mean +/- 1.96*SE. ' ...
    'All four conditions are saved, with whole movement also under movement_all.'];
summary.metric_definitions.error = 'Mean absolute same-time achieved-minus-selected-target error in joint-position units';
summary.error_alignment = 'position(i) and active_target(i) are both at t=(i-1)/fs; the extra final position sample is excluded';
summary.error_reference = ['Scheduled programme 1 target before 54 s; programme 2 from 54 s. ' ...
    'The scheduled moving target is retained during rest and hyperdirect activation.'];
summary.metric_definitions.velocity = 'Mean absolute 10-sample moving-mean motor output within state, multiplied by fs';
summary.metric_definitions.tracking_covariance = 'Sample covariance of scheduled active target and achieved position over the whole state';
summary.metric_definitions.tracking_correlation = 'Pearson correlation of scheduled active target and achieved position over the whole state';
summary.metric_definitions.rigidity_covariance = 'Sample covariance of summed unsmoothed gain-scaled flexor and extensor outputs over the whole state';
summary.metric_definitions.rigidity_correlation = 'Pearson correlation of summed unsmoothed flexor and extensor outputs over the whole state';
summary.movement_all.condition_labels = cfg.condition_labels;
summary.movement_all.state_names = {'movement'};
summary.movement_all.state_labels = {'All movement'};
summary.movement_all.state_times_s = results{1}.movement_times_s;
summary.movement_all.metric_names = cfg.metric_names;
summary.movement_all.completed_realisations = R;
summary.movement_all.metrics_layout = 'condition x realisation';
summary.movement_all.statistics_layout = 'condition x 1';
summary.movement_all.note = ['Metrics are recalculated over the full 6-78 s interval, ' ...
    'not averaged from the two movement states.'];
summary.movement_all.metric_definitions = summary.metric_definitions;
latencyValues = nan(C, R);
for r = 1:R
    latencyValues(:, r) = results{r}.latency_single_s;
end
summary.latency_single = summarise_single_programme_latency( ...
    latencyValues, results{1}.latency_single_window_s, cfg.condition_labels);
for k = 1:numel(cfg.metric_names)
    field = cfg.metric_names{k};
    values = nan(C, S, R);
    for r = 1:R
        assert(isequal(size(results{r}.metrics.(field)), [C S]), ...
            'Unexpected shape for per-realisation state metrics.');
        values(:, :, r) = results{r}.metrics.(field);
    end
    summary.metrics.(field) = values;
    stats = statistics_over_runs(values, 3);
    summary.stats.(field) = stats;
    % Convenience views with the same condition x state indexing.
    summary.mean.(field) = stats.mean;
    summary.std.(field) = stats.std;
    summary.se.(field) = stats.se;
    summary.ci95_lower.(field) = stats.ci95_lower;
    summary.ci95_upper.(field) = stats.ci95_upper;
    summary.ci95_halfwidth.(field) = stats.ci95_halfwidth;
    allMovement = nan(C, 1, R);
    for r = 1:R
        allMovement(:, 1, r) = results{r}.movement_all_metrics.(field);
    end
    summary.movement_all.metrics.(field) = reshape(allMovement, C, R);
    summary.movement_all.stats.(field) = statistics_over_runs(allMovement, 3);
end

% Joint frequency distributions: keep condition and state separate and
% average each frequency bin over independent per-realisation spectra.
frequency = results{1}.joint_spectra{1}.frequency_hz;
F = numel(frequency);
summary.joint_spectra.frequency_hz = frequency;
summary.joint_spectra.layout = 'condition x state x realisation x frequency';
summary.joint_spectra.statistics_layout = 'condition x state x frequency';
summary.joint_spectra.wavelet = cfg.wavelet;
summary.joint_spectra.state_sample_count = results{1}.joint_spectra{1}.state_sample_count;
summary.joint_spectra.state_names = cfg.state_names;
summary.joint_spectra.state_times_s = summary.state_times_s;
summary.joint_spectra.magnitude_definition = 'Time mean of abs(CWT(joint position)) within each state';
summary.joint_spectra.power_definition = 'Time mean of abs(CWT(joint position)).^2; not PSD per Hz';
summary.joint_spectra.averaging = 'Each run contributes one distribution per state; runs are equally weighted';
summary.joint_spectra.boundary_note = ['magnitude/power include all state samples; power_coi excludes ' ...
    'full-record edge effects but not mixing across state transitions. Inspect valid sample counts.'];
summary.movement_all.joint_spectra.frequency_hz = frequency;
summary.movement_all.joint_spectra.wavelet = cfg.wavelet;
summary.movement_all.joint_spectra.layout = 'condition x realisation x frequency';
summary.movement_all.joint_spectra.statistics_layout = 'condition x frequency';
summary.movement_all.joint_spectra.interval_s = results{1}.movement_times_s;
summary.movement_all.joint_spectra.state_sample_count = ...
    results{1}.joint_spectra{1}.movement_all.state_sample_count;
summary.movement_all.joint_spectra.power_definition = summary.joint_spectra.power_definition;
summary.movement_all.joint_spectra.boundary_note = summary.joint_spectra.boundary_note;
quantities = {'magnitude', 'power', 'power_coi'};
for q = 1:numel(quantities)
    name = quantities{q};
    values = nan(C, S, R, F);
    allMovement = nan(C, R, F);
    for r = 1:R
        for c = 1:C
            js = results{r}.joint_spectra{c};
            assert(isequal(js.frequency_hz, frequency), ...
                'Joint frequency grids differ between runs or conditions.');
            assert(isequal(size(js.(name)), [S F]), 'Unexpected joint spectrum shape.');
            values(c, :, r, :) = reshape(js.(name), [1 S 1 F]);
            allMovement(c, r, :) = reshape(js.movement_all.(name), [1 1 F]);
        end
    end
    summary.joint_spectra.(name) = values;
    summary.joint_spectra.stats.(name) = statistics_over_runs(values, 3);
    summary.movement_all.joint_spectra.(name) = allMovement;
    summary.movement_all.joint_spectra.stats.(name) = statistics_over_runs(allMovement, 2);
end
counts = zeros(C, S, R, F);
allMovementCounts = zeros(C, R, F);
for r = 1:R
    for c = 1:C
        counts(c, :, r, :) = reshape(results{r}.joint_spectra{c}.coi_valid_time_samples, [1 S 1 F]);
        allMovementCounts(c, r, :) = reshape( ...
            results{r}.joint_spectra{c}.movement_all.coi_valid_time_samples, [1 1 F]);
    end
end
summary.joint_spectra.coi_valid_time_samples = counts;
summary.movement_all.joint_spectra.coi_valid_time_samples = allMovementCounts;

% Retain the previous scalar definitions and last-programme thalamic signal.
% The first-programme window is shifted by the inserted 24-second interval.
summary.legacy.tracking_times_s = results{1}.tracking_times_s;
fields = fieldnames(results{1}.legacy_metrics{1});
for k = 1:numel(fields)
    field = fields{k};
    values = zeros(R, C);
    for r = 1:R
        for c = 1:C
            values(r, c) = results{r}.legacy_metrics{c}.(field);
        end
    end
    summary.legacy.metrics.(field) = values;
    summary.legacy.stats.(field) = statistics_over_runs(values, 1);
end
names = {'restdistribution', 'movementdistribution', 'hyperdirectdistribution'};
for c = 1:C
    frequency = results{1}.legacy_spectra{c}.freq;
    summary.legacy.thalamus_spectra(c).frequency = frequency;
    for k = 1:numel(names)
        values = zeros(R, numel(frequency));
        for r = 1:R
            assert(isequal(results{r}.legacy_spectra{c}.freq, frequency), ...
                'Frequency grids differ between realisations.');
            values(r, :) = results{r}.legacy_spectra{c}.(names{k})(:).';
        end
        summary.legacy.thalamus_spectra(c).runs.(names{k}) = values;
        summary.legacy.thalamus_spectra(c).stats.(names{k}) = statistics_over_runs(values, 1);
    end
end
[summary.metric_values_table, summary.metric_statistics_table] = metric_tables(summary);
% Reuse the table builder with a singleton state axis for whole movement.
allMovementTableInput = summary.movement_all;
for k = 1:numel(cfg.metric_names)
    name = cfg.metric_names{k};
    allMovementTableInput.metrics.(name) = reshape(summary.movement_all.metrics.(name), C, 1, R);
end
[summary.movement_all.metric_values_table, summary.movement_all.metric_statistics_table] = ...
    metric_tables(allMovementTableInput);
end

function stats = statistics_over_runs(values, runDimension)
% Reduce only the realisation axis. Keep all conditions/states/frequencies.
shape = size(values);
nd = max(numel(shape), runDimension);
shape(end+1:nd) = 1;
otherDimensions = setdiff(1:nd, runDimension, 'stable');
outputShape = shape(otherDimensions);
if numel(outputShape) == 1
    outputShape = [1 outputShape];
end
ordered = permute(values, [runDimension otherDimensions]);
flat = reshape(ordered, shape(runDimension), []);
valid = isfinite(flat);
counts = sum(valid, 1);
clean = flat;
clean(~valid) = 0;
average = sum(clean, 1) ./ max(counts, 1);
average(counts == 0) = NaN;
deviations = bsxfun(@minus, flat, average);
deviations(~valid) = 0;
sampleSTD = sqrt(sum(deviations.^2, 1) ./ max(counts-1, 1));
sampleSTD(counts < 2) = NaN;
standardError = sampleSTD ./ sqrt(max(counts, 1));
halfwidth = 1.96 * standardError;
stats.mean = reshape(average, outputShape);
stats.std = reshape(sampleSTD, outputShape);
stats.se = reshape(standardError, outputShape);
stats.ci95_lower = reshape(average-halfwidth, outputShape);
stats.ci95_upper = reshape(average+halfwidth, outputShape);
stats.ci95_halfwidth = reshape(halfwidth, outputShape);
stats.n_valid = reshape(counts, outputShape);
stats.n_missing = reshape(shape(runDimension)-counts, outputShape);
stats.n_expected = shape(runDimension);
stats.ci_multiplier = 1.96;
stats.ci_method = 'Mean +/- 1.96*SE (pointwise approximate 95% normal interval)';
end

function latency = summarise_single_programme_latency(values, window, conditionLabels)
C = size(values, 1);
R = size(values, 2);
latency.values_s = values;
latency.layout = 'condition x realisation';
latency.condition_labels = conditionLabels;
latency.state_name = 'movement_single';
latency.window_s = window;
latency.method_version = 1;
latency.method = 'Lag at maximum mean-centred xcorr(position,target), divided by fs; positive means position lags target';
latency.note = 'Constant/nonfinite signals give NaN. Periodic targets can give ambiguous lags; no phase wrapping or lag restriction is applied.';
latency.stats = statistics_over_runs(reshape(values, C, 1, R), 3);

[ConditionIndex, Realisation] = ndgrid(1:C, 1:R);
ConditionIndex = ConditionIndex(:); Realisation = Realisation(:);
Condition = reshape(conditionLabels(ConditionIndex), [], 1);
Latency_s = values(:);
latency.values_table = table(ConditionIndex, Condition, Realisation, Latency_s);
ConditionIndex = (1:C).';
Condition = reshape(conditionLabels, [], 1);
stats = latency.stats;
numbers = [stats.mean(:), stats.std(:), stats.se(:), ...
    stats.ci95_lower(:), stats.ci95_upper(:), stats.n_valid(:), ...
    stats.n_missing(:), repmat(R, C, 1), repmat(window, C, 1)];
latency.statistics_table = [table(ConditionIndex, Condition), ...
    array2table(numbers, 'VariableNames', {'Mean_s', 'STD_s', 'SE_s', ...
    'CI95_lower_s', 'CI95_upper_s', 'N_valid', 'N_missing', 'N_expected', 'Start_s', 'End_s'})];
end

function summary = ensure_single_programme_latency(summary, outputDir)
% Recover the new 6-30 s metric from saved trajectories, without model calls.
C = numel(summary.condition_labels);
R = summary.completed_realisations;
singleState = find(strcmp(summary.state_names, 'movement_single'));
window = summary.state_times_s(singleState, :);
hasValues = isfield(summary, 'latency_single') ...
    && isfield(summary.latency_single, 'values_s') ...
    && isfield(summary.latency_single, 'method_version') ...
    && summary.latency_single.method_version == 1 ...
    && isfield(summary.latency_single, 'window_s') ...
    && isequal(summary.latency_single.window_s, window) ...
    && isequal(size(summary.latency_single.values_s), [C R]);
if hasValues
    values = summary.latency_single.values_s;
else
    files = cell(R, 1);
    for r = 1:R
        files{r} = fullfile(outputDir, sprintf('realisation_%03d.mat', r));
        if ~exist(files{r}, 'file')
            error('Fig2:MissingLatencyData', ...
                ['Single-programme latency needs the saved trajectory from every run. ' ...
                'Restore the missing file: %s. No simulations have been started.'], files{r});
        end
    end
    if isempty(which('xcorr'))
        error('Fig2:MissingLatencyDependency', 'xcorr is required to calculate latency from saved traces.');
    end
    values = nan(C, R);
    fprintf('Calculating single-programme latency from %d saved realisations.\n', R);
    for r = 1:R
        savedRun = load(files{r}, 'runData');
        runData = savedRun.runData;
        assert(runData.realisation == r && numel(runData.conditions) == C, ...
            'Saved realisation index or condition count does not match summary.mat.');
        assert(runData.config.fs == summary.config.fs ...
            && runData.config.base_seed == summary.config.base_seed ...
            && isequal(runData.config.condition_labels, summary.condition_labels), ...
            'Saved realisation configuration does not match the summary.');
        assert(isequal(runData.inputs.phase_times_s(singleState, :), window), ...
            'Single-programme windows differ between the saved run and summary.');
        for c = 1:C
            values(c, r) = calculate_single_programme_latency( ...
                runData.conditions{c}, runData.inputs, runData.config.fs);
        end
        if mod(r, 8) == 0 || r == R
            fprintf('Analysed saved realisation %d of %d.\n', r, R);
        end
    end
end
summary.latency_single = summarise_single_programme_latency(values, window, summary.condition_labels);
end

function summary = refresh_summary_intervals(summary)
% Update old t intervals everywhere; retain the saved means, STD and SE.
summary = normal_interval_tree(summary);
summary.analysis_schema_version = 8;
summary.config.ci_multiplier = 1.96;
summary.interval_method = ['Arithmetic mean, sample STD, SE=STD/sqrt(n_valid), ' ...
    'and pointwise mean +/- 1.96*SE (approximate 95% normal interval)'];
summary.undefined_correlation_policy = ['Constant signals give NaN. Aggregation omits ' ...
    'undefined values and reports n_valid and n_missing; SE uses n_valid.'];
summary.figure_note = ['Trajectories show the realisations identified in their row titles. Plots show condition indices 1, 2 and 4. ' ...
    'All error bars/bands are mean +/- 1.96*SE across realisations. ' ...
    'The latency panel uses single-programme movement only.'];
for k = 1:numel(summary.metric_names)
    name = summary.metric_names{k};
    stats = summary.stats.(name);
    summary.ci95_lower.(name) = stats.ci95_lower;
    summary.ci95_upper.(name) = stats.ci95_upper;
    summary.ci95_halfwidth.(name) = stats.ci95_halfwidth;
end
[summary.metric_values_table, summary.metric_statistics_table] = metric_tables(summary);
C = numel(summary.condition_labels);
R = summary.completed_realisations;
allMovement = summary.movement_all;
for k = 1:numel(summary.metric_names)
    name = summary.metric_names{k};
    allMovement.metrics.(name) = reshape(summary.movement_all.metrics.(name), C, 1, R);
end
[summary.movement_all.metric_values_table, summary.movement_all.metric_statistics_table] = ...
    metric_tables(allMovement);
end

function summary = add_frequency_weighted_power(summary)
% Derive from ORIGINAL per-run power on every call; this is safe to repeat.
% Multiplication by a fixed positive frequency scales mean, STD, SE and
% interval bounds equally, without changing the number of valid runs.
spectra = summary.joint_spectra;
frequency = spectra.frequency_hz(:).';
assert(all(isfinite(frequency) & frequency > 0), ...
    'Frequency weighting requires finite, positive frequencies in Hz.');
assert(size(spectra.power, 4) == numel(frequency), ...
    'The saved state power frequency axis does not match frequency_hz.');
spectra.frequency_weighted_power = bsxfun(@times, spectra.power, ...
    reshape(frequency, 1, 1, 1, []));
spectra.stats.frequency_weighted_power = ...
    statistics_over_runs(spectra.frequency_weighted_power, 3);
spectra.frequency_weighted_power_definition = ...
    'frequency_hz .* power at each frequency; original power retained; units position^2 Hz';
spectra.plot_quantity = 'frequency_weighted_power';
spectra.plot_note = ['Linear axes; one condition per panel; each panel has its own y-axis range. ' ...
    'Bands are mean +/- 1.96*SE across runs, without clipping negative lower bounds. ' ...
    'Frequency weighting is not a normalisation or a sum across frequency bins.'];
if isfield(spectra, 'log_plot_note')
    spectra = rmfield(spectra, 'log_plot_note');
end
summary.joint_spectra = spectra;

allMovement = summary.movement_all.joint_spectra;
assert(isequal(allMovement.frequency_hz(:).', frequency), ...
    'Whole-movement and state power must use the same frequency grid.');
assert(size(allMovement.power, 3) == numel(frequency), ...
    'The saved whole-movement power frequency axis does not match frequency_hz.');
allMovement.frequency_weighted_power = bsxfun(@times, allMovement.power, ...
    reshape(frequency, 1, 1, []));
allMovement.stats.frequency_weighted_power = ...
    statistics_over_runs(allMovement.frequency_weighted_power, 2);
allMovement.frequency_weighted_power_definition = spectra.frequency_weighted_power_definition;
summary.movement_all.joint_spectra = allMovement;
summary.analysis_schema_version = 8;
summary.figure_note = ['Trajectories show the realisations identified in their row titles. Plots show condition indices 1, 2 and 4. ' ...
    'Spectra show frequency-weighted power on independent linear y-axes for each condition/state. ' ...
    'Rest and movement flexor-extensor correlations have separate, independent y-axis ranges. ' ...
    'All error bars/bands are mean +/- 1.96*SE across realisations. ' ...
    'The latency panel uses single-programme movement only. All four conditions are saved, ' ...
    'with whole movement also under movement_all.'];
end

function node = normal_interval_tree(node)
% Numeric statistic structs appear under metrics, spectra, legacy and latency.
for j = 1:numel(node)
    if isfield(node, 'mean') && isfield(node, 'se') ...
            && isnumeric(node(j).mean) && isnumeric(node(j).se)
        node(j).ci95_halfwidth = 1.96 * node(j).se;
        node(j).ci95_lower = node(j).mean - node(j).ci95_halfwidth;
        node(j).ci95_upper = node(j).mean + node(j).ci95_halfwidth;
        node(j).ci_multiplier = 1.96;
        node(j).ci_method = 'Mean +/- 1.96*SE (pointwise approximate 95% normal interval)';
    end
    names = fieldnames(node(j));
    for k = 1:numel(names)
        name = names{k};
        if isstruct(node(j).(name))
            node(j).(name) = normal_interval_tree(node(j).(name));
        end
    end
end
end

function save_latency_tables(summary, outputDir)
writetable(summary.latency_single.values_table, fullfile(outputDir, 'latency_single_per_realisation.csv'));
writetable(summary.latency_single.statistics_table, fullfile(outputDir, 'latency_single_summary.csv'));
end

function [valueTable, statisticsTable] = metric_tables(summary)
C = numel(summary.condition_labels);
S = numel(summary.state_names);
R = summary.completed_realisations;
M = numel(summary.metric_names);
N = C*S*R;
ConditionIndex = zeros(N, 1); StateIndex = zeros(N, 1); Realisation = zeros(N, 1);
metricValues = nan(N, M);
row = 0;
for c = 1:C
    for s = 1:S
        for r = 1:R
            row = row + 1;
            ConditionIndex(row) = c; StateIndex(row) = s; Realisation(row) = r;
            for m = 1:M
                metricValues(row, m) = summary.metrics.(summary.metric_names{m})(c, s, r);
            end
        end
    end
end
Condition = reshape(summary.condition_labels(ConditionIndex), [], 1);
State = reshape(summary.state_names(StateIndex), [], 1);
valueTable = [table(ConditionIndex, Condition, StateIndex, State, Realisation), ...
    array2table(metricValues, 'VariableNames', summary.metric_names)];

N = C*S*M;
ConditionIndex = zeros(N, 1); StateIndex = zeros(N, 1); Metric = cell(N, 1);
numbers = nan(N, 8);
row = 0;
for c = 1:C
    for s = 1:S
        for m = 1:M
            row = row + 1;
            name = summary.metric_names{m};
            stat = summary.stats.(name);
            ConditionIndex(row) = c; StateIndex(row) = s; Metric{row} = name;
            numbers(row, :) = [stat.mean(c,s), stat.std(c,s), stat.se(c,s), ...
                stat.ci95_lower(c,s), stat.ci95_upper(c,s), ...
                stat.n_valid(c,s), stat.n_missing(c,s), stat.n_expected];
        end
    end
end
Condition = reshape(summary.condition_labels(ConditionIndex), [], 1);
State = reshape(summary.state_names(StateIndex), [], 1);
statisticsTable = [table(ConditionIndex, Condition, StateIndex, State, Metric), ...
    array2table(numbers, 'VariableNames', ...
    {'Mean', 'STD', 'SE', 'CI95_lower', 'CI95_upper', 'N_valid', 'N_missing', 'N_expected'})];
end

function [indices, mode] = resolve_plot_realisations(selection, nRuns, nRows)
% Choose illustrative traces only. Simulation RNG substreams are unchanged.
if ischar(selection) || isstring(selection)
    isSingleText = (ischar(selection) && isrow(selection)) ...
        || (isstring(selection) && isscalar(selection));
    if ~isSingleText || ~strcmpi(strtrim(char(selection)), 'random')
        error('Fig2:InvalidPlotRealisations', ...
            'Use ''random'', one realisation number, or %d numbers in displayed-condition order.', nRows);
    end
    indices = randi(nRuns, 1, nRows);
    mode = 'random';
    return;
end
validateattributes(selection, {'numeric'}, ...
    {'vector', 'nonempty', 'real', 'finite', 'integer', '>=', 1, '<=', nRuns}, ...
    mfilename, 'plotRealisations');
if isscalar(selection)
    indices = repmat(double(selection), 1, nRows);
elseif numel(selection) == nRows
    indices = double(selection(:).');
else
    error('Fig2:InvalidPlotRealisations', ...
        'Specify one realisation number for all rows, or exactly %d numbers, one per displayed condition.', nRows);
end
mode = 'specified';
end

function summary = record_plot_realisations(summary)
% Save the exact choices so a random figure can be reproduced with a vector.
indices = summary.config.plot_realisations;
conditions = summary.config.plot_condition_indices;
labels = summary.condition_labels(conditions);
summary.plot_realisations = indices;
summary.plot_realisation_selection = summary.config.plot_realisation_selection;
summary.figure_examples.row_indices = 1:numel(conditions);
summary.figure_examples.condition_indices = conditions;
summary.figure_examples.condition_labels = labels;
summary.figure_examples.realisation_indices = indices;
summary.figure_examples.selection_mode = summary.plot_realisation_selection;
summary.figure_examples.note = ['These choices affect the three illustrative trajectory rows only. ' ...
    'Spectra and bar plots still summarise all realisations. Random draws may repeat.'];
summary.figure_examples.files = cell(1, numel(conditions));
for c = 1:numel(conditions)
    summary.figure_examples.files{c} = fullfile(summary.output_dir, ...
        sprintf('realisation_%03d.mat', indices(c)));
end
fprintf('Trajectory realisations (%s): %s.\n', strjoin(labels, ' / '), mat2str(indices));
end

function examples = load_figure_examples(summary, outputDir)
% Load each selected run once; retain only the traces and common inputs.
% Do not retain all conditions/raw arrays from three complete run files.
indices = summary.plot_realisations;
conditions = summary.config.plot_condition_indices;
selectedRuns = unique(indices, 'stable');
for r = selectedRuns
    filename = fullfile(outputDir, sprintf('realisation_%03d.mat', r));
    if ~exist(filename, 'file')
        error('Fig2:MissingPlotRealisation', ...
            'Selected trajectory realisation %d is missing. Restore this saved file: %s', r, filename);
    end
end
examples.muscleposition = cell(1, numel(conditions));
examples.realisation_indices = indices;
examples.condition_indices = conditions;
examples.inputs = [];
for r = selectedRuns
    saved = load(fullfile(outputDir, sprintf('realisation_%03d.mat', r)), 'runData');
    runData = saved.runData;
    assert(runData.realisation == r && numel(runData.conditions) == numel(summary.condition_labels), ...
        'Selected run index or condition count does not match the summary.');
    assert(runData.config.fs == summary.config.fs ...
        && runData.config.base_seed == summary.config.base_seed ...
        && isequal(runData.config.condition_labels, summary.condition_labels) ...
        && isequal(runData.inputs.phase_times_s, summary.state_times_s), ...
        'Selected run configuration or state windows do not match the summary.');
    if isempty(examples.inputs)
        examples.inputs = runData.inputs;
    else
        assert(isequaln(examples.inputs, runData.inputs), ...
            'Selected trajectory realisations have different simulation inputs.');
    end
    for c = find(indices == r)
        position = runData.conditions{conditions(c)}.muscleposition;
        assert(numel(position) == numel(examples.inputs.ttime), ...
            'The selected joint-position trajectory does not match the saved time axis.');
        examples.muscleposition{c} = position;
    end
    clear saved runData;
end
end

function save_result_figures(summary, outputDir)
% Shared by a new simulation and plot-only mode. Numerical data are read-only.
examples = load_figure_examples(summary, outputDir);
fig = plot_figure2(examples, summary);
set_figure_style(fig, 12, 1);
savefig(fig, fullfile(outputDir, 'figure2_48_realisations.fig'));
print(fig, fullfile(outputDir, 'figure2_48_realisations.png'), '-dpng', '-r200');
metricFig = plot_state_metrics(summary);
set_figure_style(metricFig, 12, 1);
savefig(metricFig, fullfile(outputDir, 'metrics_by_state.fig'));
print(metricFig, fullfile(outputDir, 'metrics_by_state.png'), '-dpng', '-r200');
end

function fig = plot_figure2(examples, summary)
cfg = summary.config;
inputs = examples.inputs;
t = inputs.ttime;
plotConditions = cfg.plot_condition_indices;
colors = cfg.condition_colors(plotConditions, :);
labels = cfg.condition_labels(plotConditions);
shortLabels = cfg.condition_short_labels(plotConditions);
C = numel(plotConditions);
S = numel(cfg.state_names);
% Three displayed trajectories and velocity on the left; one spectral row
% per state on the right, with Healthy / PD / PD severe in separate panels.
% Four bottom panels separate rest rigidity from movement rigidity.
topRows = max(C+1, S);
columns = 12;
leftSpan = 6;
spectrumSpan = 2;
metricSpan = 3;
% Set both screen geometry and print area so exports use 160 x 200 mm.
fig = figure('Color', 'w', 'WindowStyle', 'normal', ...
    'Units', 'centimeters', 'Position', [2 2 16 20], ...
    'PaperUnits', 'centimeters', 'PaperSize', [16 20], ...
    'PaperPosition', [0 0 16 20], 'PaperPositionMode', 'manual');
layout = tiledlayout(fig, topRows+1, columns, 'TileSpacing', 'compact', 'Padding', 'compact');
phase = inputs.phase_times_s;
restState = find(strcmp(cfg.state_names, 'rest'));
hyperdirectState = find(strcmp(cfg.state_names, 'hyperdirect'));
for c = 1:C
    ax = nexttile(layout, 1 + (c-1)*columns, [1 leftSpan]);
    hold(ax, 'on');
    patch(ax, phase(restState, [1 2 2 1]), ...
        [0 0 1 1], [0.5 0 0], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    patch(ax, phase(hyperdirectState, [1 2 2 1]), ...
        [0 0 1 1], [0.5 0 0], 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    idx = 1:5:numel(t);
    plot(ax, t(idx), examples.muscleposition{c}(idx), ...
        'Color', colors(c,:), 'LineWidth', 1.5);
    styles = {'--', ':'};
    step = 500;
    for p = 1:2
        for k = 1:step:numel(t)-1
            last = min(k+step, numel(t));
            priority = mean(inputs.motor_drive(p, k:last)) / cfg.sma_active;
            priority = max(0, min(1, priority));
            targetSegment = inputs.desiredposition(p, k:last);
            targetSegment(~inputs.program_enabled(p, k:last)) = NaN;
            plot(ax, t(k:last), targetSegment, ...
                'Color', [1-priority priority 0], 'LineStyle', styles{p}, 'LineWidth', 1.2);
        end
    end
    xline(ax, inputs.two_programme_enable_s, ':k', 'Programme 2 enabled');
    xline(ax, inputs.target_switch_s, '--k', 'Target switch');
    ylim(ax, [0 1]); xlim(ax, [t(1) t(end)]);
    title(ax, sprintf('%s - illustrative realisation %d', labels{c}, examples.realisation_indices(c)));
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Joint position');
end

titles = cfg.state_labels;
frequency = summary.joint_spectra.frequency_hz;
selected = find(frequency >= 1 & frequency <= 15);
frequency = frequency(selected);
for k = 1:S
    for c = 1:C
        conditionIndex = plotConditions(c);
        tileIndex = leftSpan + 1 + (k-1)*columns + (c-1)*spectrumSpan;
        ax = nexttile(layout, tileIndex, [1 spectrumSpan]); hold(ax, 'on');
        set(ax, 'YScale', 'linear');
        stats = summary.joint_spectra.stats.frequency_weighted_power;
        m = reshape(stats.mean(conditionIndex, k, selected), 1, []);
        ci = 1.96 * reshape(stats.se(conditionIndex, k, selected), 1, []);
        lower = m-ci; upper = m+ci;
        % Use only this condition/state to set the linear y-axis limits.
        % Keep finite negative CI limits; do not clip bands at zero.
        yLimits = independent_linear_limits([m lower upper]);
        patch(ax, [3 7 7 3], yLimits([1 1 2 2]), [0.7 0.7 0.7], ...
            'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot_linear_band(ax, frequency, lower, upper, colors(c,:));
        plot(ax, frequency, m, 'Color', colors(c,:), 'LineWidth', 1.5);
        xlim(ax, [1 15]); ylim(ax, yLimits);
        xticks(ax, [1 5 10 15]);
        title(ax, {titles{k}, labels{c}});
        xlabel(ax, 'Frequency (Hz)');
        ylabel(ax, 'f \times power (position^2 Hz)');
    end
end

ax = nexttile(layout, 1 + C*columns, [1 leftSpan]);
means = summary.stats.velocity.mean(plotConditions, :).';
ci = 1.96 * summary.stats.velocity.se(plotConditions, :).';
bars = bar(ax, means, 'grouped'); hold(ax, 'on');
for c = 1:C
    bars(c).FaceColor = colors(c,:);
    errorbar(ax, bars(c).XEndPoints, means(:,c).', ci(:,c).', ...
        'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
end
xticks(ax, 1:S); xticklabels(ax, cfg.state_labels);
title(ax, 'Velocity'); ylabel(ax, 'Units/s');
legend(ax, bars, shortLabels, 'Location', 'best');
metricNames = {'rigidity_correlation', 'rigidity_correlation', ...
    'tracking_correlation', 'latency_single'};
metricTitles = {{'Flexor-extensor correlation', 'Rest'}, ...
    {'Flexor-extensor correlation', 'Movement'}, 'Target-position correlation', ...
    sprintf('Single-programme latency (%g-%g s)', summary.latency_single.window_s)};
metricUnits = {'Pearson r', 'Pearson r', 'Pearson r', 'Seconds'};
panelStates = {restState, ...
    find(ismember(cfg.state_names, {'movement_single', 'movement_two'})), ...
    find(ismember(cfg.state_names, {'movement_single', 'movement_two'})), []};
for k = 1:numel(metricNames)
    tileIndex = 1 + topRows*columns + (k-1)*metricSpan;
    ax = nexttile(layout, tileIndex, [1 metricSpan]); hold(ax, 'on');
    if strcmp(metricNames{k}, 'latency_single')
        values = reshape(summary.latency_single.stats.mean(plotConditions), 1, []);
        interval = 1.96 * reshape(summary.latency_single.stats.se(plotConditions), 1, []);
    else
        selectedStates = panelStates{k};
        values = summary.stats.(metricNames{k}).mean(plotConditions, selectedStates).';
        interval = 1.96 * summary.stats.(metricNames{k}).se(plotConditions, selectedStates).';
    end
    if size(values, 1) == 1
        % A single state uses one flat-colour series with three condition bars.
        b = bar(ax, 1:C, values, 'FaceColor', 'flat'); b.CData = colors;
        errorbar(ax, 1:C, values, interval, 'k', 'LineStyle', 'none');
        xticks(ax, 1:C); xticklabels(ax, shortLabels);
    else
        b = bar(ax, values, 'grouped');
        for c = 1:C
            b(c).FaceColor = colors(c,:);
            errorbar(ax, b(c).XEndPoints, values(:,c).', interval(:,c).', ...
                'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
        end
        xticks(ax, 1:numel(selectedStates));
        xticklabels(ax, cfg.state_labels(selectedStates));
    end
    % Rest and movement get independent limits, including their own CIs.
    % Keeping zero visible avoids a truncated baseline for the bars.
    ylim(ax, independent_linear_limits([values(:); values(:)-interval(:); values(:)+interval(:)]));
    title(ax, metricTitles{k}); ylabel(ax, metricUnits{k});
end
title(layout, sprintf('Figure 2: mean +/- 1.96 SE across %d realisations; trace runs %s (Healthy / PD / PD severe)', ...
    summary.completed_realisations, mat2str(examples.realisation_indices)));
end

function fig = plot_state_metrics(summary)
plotConditions = summary.config.plot_condition_indices;
colors = summary.config.condition_colors(plotConditions, :);
labels = summary.condition_labels(plotConditions);
names = summary.metric_names;
titles = {'Target error', 'Velocity', 'Target-position covariance', ...
    'Target-position correlation', 'Flexor-extensor covariance', 'Flexor-extensor correlation'};
units = {'Mean absolute error', 'Units/s', 'Covariance', 'Pearson r', 'Covariance', 'Pearson r'};
fig = figure('Color', 'w', 'Position', [70 70 1800 950]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(names)
    ax = nexttile(layout); hold(ax, 'on');
    means = summary.stats.(names{k}).mean(plotConditions, :).';
    intervals = 1.96 * summary.stats.(names{k}).se(plotConditions, :).';
    bars = bar(ax, means, 'grouped');
    for c = 1:numel(plotConditions)
        bars(c).FaceColor = colors(c, :);
        errorbar(ax, bars(c).XEndPoints, means(:, c).', intervals(:, c).', ...
            'k', 'LineStyle', 'none', 'HandleVisibility', 'off');
    end
    xticks(ax, 1:numel(summary.state_labels));
    xticklabels(ax, summary.state_labels);
    title(ax, titles{k}); ylabel(ax, units{k});
    if k == 1
        legend(ax, bars, labels, 'Location', 'best');
    end
end
title(layout, sprintf('State metrics: mean +/- 1.96 SE across %d runs; undefined correlations omitted', ...
    summary.completed_realisations));
end

function set_figure_style(fig, fontSize, gridLineWidth)
% Apply exact point sizes to axes, labels, titles, legends and line labels.
% GridLineWidth/MinorGridLineWidth are available from MATLAB R2023a.
% Earlier releases use axes LineWidth for the grid, ticks and outline.
% https://www.mathworks.com/help/matlab/ref/matlab.graphics.axis.axes-properties.html
set(findall(fig, '-property', 'FontUnits'), 'FontUnits', 'points');
axesHandles = findall(fig, 'Type', 'axes');
set(axesHandles, 'LabelFontSizeMultiplier', 1, 'TitleFontSizeMultiplier', 1);
set(findall(fig, '-property', 'FontSize'), 'FontSize', fontSize);
for k = 1:numel(axesHandles)
    ax = axesHandles(k);
    set(ax, 'XGrid', 'on', 'YGrid', 'on', ...
        'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
        'GridLineStyle', '-', 'MinorGridLineStyle', ':');
    if isprop(ax, 'GridLineWidth') && isprop(ax, 'MinorGridLineWidth')
        ax.GridLineWidth = gridLineWidth;
        ax.MinorGridLineWidth = gridLineWidth;
    else
        ax.LineWidth = gridLineWidth;
    end
end
end

function limits = independent_linear_limits(values)
% Each call sees one panel only. Include zero and all finite interval bounds.
values = values(isfinite(values));
lower = min([0; values(:)]);
upper = max([0; values(:)]);
span = upper-lower;
if span == 0
    limits = [0 1];
    return;
end
if lower < 0
    lower = lower - 0.08*span;
end
limits = [lower upper+0.08*span];
end

function plot_linear_band(ax, frequency, lower, upper, color)
% Include finite bounds of either sign on the linear axes. Separate valid
% segments so a patch cannot bridge a missing/undefined frequency bin.
valid = isfinite(frequency) & isfinite(lower) & isfinite(upper);
transitions = diff([false valid false]);
starts = find(transitions == 1);
stops = find(transitions == -1) - 1;
for k = 1:numel(starts)
    idx = starts(k):stops(k);
    if numel(idx) >= 2
        fill(ax, [frequency(idx) fliplr(frequency(idx))], ...
            [lower(idx) fliplr(upper(idx))], color, 'FaceAlpha', 0.15, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    else
        plot(ax, [frequency(idx) frequency(idx)], [lower(idx) upper(idx)], ...
            'Color', color, 'HandleVisibility', 'off');
    end
end
end
