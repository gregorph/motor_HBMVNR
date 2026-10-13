function results = for_network_training_corrected(options)
%FOR_NETWORK_TRAINING_CORRECTED Objective-guided search with current DBS APIs.
%   results = for_network_training_corrected();
%   cfg = for_network_training_corrected('defaults');
%   cfg.use_parallel = true;   % requires Parallel Computing Toolbox
%   results = for_network_training_corrected(cfg);
%
% Baseline: 8 independent trials of identical initial settings.
% Search: 24 candidates, each changing one parameter from the retained state.
% Validation: 8 independent trials of the identical selected candidate.
% The retained state changes only after validation succeeds. Pending validation
% is completed even if the nominal iteration budget has just been reached.
%
% DBS normal proposals are increments. Redraw until the resulting amplitude
% is strictly positive; negative increments may reduce an existing amplitude.
% Zero amplitudes in the initial state continue to represent no stimulation.
%
% Required project helpers: spinal_circuits_fcn_pts, put_in_range and
% direction_dbs_logistic_curve. See README.md for the bundled DBS modules.

cfg = default_configuration();
if nargin == 1 && (ischar(options) || (isstring(options) && isscalar(options)))
    if strcmp(options, 'defaults')
        results = cfg;
        return
    end
    error('OGL:InvalidOptions', 'The only text option is ''defaults''.');
end
if nargin >= 1 && ~isempty(options)
    validateattributes(options, {'struct'}, {'scalar'});
    fields = fieldnames(options);
    for k = 1:numel(fields)
        name = fields{k};
        if ~isfield(cfg, name)
            error('OGL:UnknownOption', 'Unknown configuration field: %s', name);
        end
        cfg.(name) = options.(name);
    end
end
cfg = validate_configuration(cfg);

parameter_names = { ...
    'mod_excitation_bg', 'mod_inhibition_bg', ...
    'mod_excitation_cortex', 'mod_inhibition_cortex', ...
    'mod_excitation_spinal', 'mod_inhibition_spinal', ...
    'mod_excitation_cerebellum', 'mod_inhibition_cerebellum', ...
    'mod_excitation_brainstem', 'mod_inhibition_brainstem', ...
    'striatum_dbs', 'gpe_dbs', 'stn_dbs', 'gpi_dbs', ...
    'thalamus_dbs', 'L4_dbs', 'L23_dbs', 'L5_dbs', 'DCN_dbs'};
metric_names = {'movement_error', 'movement_covariance', 'rest_movement', ...
    'movement_variance', 'rest_variance'};
number_parameters = numel(parameter_names);
dbs_indices = 11:19;

caller_rng = rng;
rng_cleanup = onCleanup(@() rng(caller_rng)); %#ok<NASGU>
rng(cfg.random_seed, 'twister');
distribution_details = repmat([cfg.initial_mu, cfg.initial_sigma], number_parameters, 1);
number_tickets = ceil(cfg.number_tickets/number_parameters)*number_parameters;
randomisation_array = randi(number_parameters, 1, number_tickets);
current_values = cfg.initial_values;
pending_values = [];
trial_counter = 0;

baseline_values = repmat(current_values, cfg.baseline_trials, 1);
[baseline_metrics, baseline_seeds, trial_counter] = run_batch( ...
    baseline_values, cfg, trial_counter);
initial_metrics = mean(baseline_metrics, 1);
reference_metrics = initial_metrics;
current_metrics = initial_metrics;
history = make_record(1, 'baseline', baseline_values, baseline_metrics, ...
    zeros(cfg.baseline_trials, 1), zeros(cfg.baseline_trials, 1), ...
    baseline_seeds, reference_metrics);
history.reference_after = reference_metrics;
if cfg.verbose
    fprintf('Baseline: %d trials with identical initial parameters.\n', cfg.baseline_trials);
end

iteration = 2;
while iteration <= cfg.iterations || ~isempty(pending_values)
    if isempty(pending_values)
        % Create the entire batch before updating any learning distributions.
        candidates = repmat(current_values, cfg.candidates_per_iteration, 1);
        changed_parameters = zeros(cfg.candidates_per_iteration, 1);
        changes = zeros(cfg.candidates_per_iteration, 1);
        draws = zeros(cfg.candidates_per_iteration, 1);
        for j = 1:cfg.candidates_per_iteration
            p = randomisation_array(randi(numel(randomisation_array)));
            changed_parameters(j) = p;
            [candidates(j, :), changes(j), draws(j)] = propose_parameter_change( ...
                current_values, p, distribution_details, dbs_indices, cfg.max_proposal_draws);
        end
        [metrics, seeds, trial_counter] = run_batch(candidates, cfg, trial_counter);
        relative_scores = objective_scores(metrics, reference_metrics, cfg);
        global_scores = objective_scores(metrics, initial_metrics, cfg);
        record = make_record(iteration, 'exploration', candidates, metrics, ...
            relative_scores, global_scores, seeds, reference_metrics);
        record.changed_parameters = changed_parameters;
        record.changes = changes;
        record.proposal_draws = draws;
        record.ticket_updates = zeros(cfg.candidates_per_iteration, 1);
        % Every exploratory candidate contributes, including candidate 1.
        for j = 1:cfg.candidates_per_iteration
            [distribution_details, randomisation_array, record.ticket_updates(j)] = ...
                update_distributions_learning_corrected(relative_scores(j), ...
                changed_parameters(j), changes(j), cfg.learning_rate, ...
                cfg.learning_rate_mu, cfg.learning_rate_sigma, ...
                randomisation_array, distribution_details, cfg.sigma_bounds);
        end
        [best_score, best_index] = max(relative_scores);
        if best_score > cfg.acceptance_threshold
            pending_values = candidates(best_index, :);
            record.selected_candidate = best_index;
        end
        if cfg.verbose
            fprintf('Iteration %d: best exploratory score %.4f; selected %d.\n', ...
                iteration, best_score, record.selected_candidate);
        end
    else
        % No mutation here: every validation row is exactly the same candidate.
        candidates = repmat(pending_values, cfg.validation_trials, 1);
        [metrics, seeds, trial_counter] = run_batch(candidates, cfg, trial_counter);
        relative_scores = objective_scores(metrics, reference_metrics, cfg);
        global_scores = objective_scores(metrics, initial_metrics, cfg);
        record = make_record(iteration, 'validation', candidates, metrics, ...
            relative_scores, global_scores, seeds, reference_metrics);
        record.beyond_iteration_budget = iteration > cfg.iterations;
        record.accepted = mean(relative_scores) > cfg.acceptance_threshold;
        if record.accepted
            current_values = pending_values;
            current_metrics = mean(metrics, 1);
            % Lower is better for all metrics except covariance (column 2).
            reference_metrics([1, 3, 4, 5]) = min( ...
                reference_metrics([1, 3, 4, 5]), current_metrics([1, 3, 4, 5]));
            reference_metrics(2) = max(reference_metrics(2), current_metrics(2));
        end
        % Rejection leaves the last validated state and its references intact.
        pending_values = [];
        if cfg.verbose
            fprintf('Iteration %d: mean validation score %.4f; accepted %d.\n', ...
                iteration, mean(relative_scores), record.accepted);
        end
    end
    record.reference_after = reference_metrics;
    history(end+1) = record; %#ok<AGROW>
    iteration = iteration + 1;
end

results.configuration = cfg;
results.parameter_names = parameter_names;
results.metric_names = metric_names;
results.initial_values = cfg.initial_values;
results.final_values = current_values;
results.final_parameters = cell2struct(num2cell(current_values), parameter_names, 2);
results.parameter_table = table(parameter_names(:), cfg.initial_values(:), ...
    current_values(:), 'VariableNames', {'Parameter', 'Initial', 'Retained'});
results.initial_metrics = initial_metrics;
results.final_metrics = current_metrics;
results.reference_metrics = reference_metrics;
results.distribution_details = distribution_details;
results.randomisation_array = randomisation_array;
results.selection_probabilities = accumarray(randomisation_array(:), 1, ...
    [number_parameters, 1])/numel(randomisation_array);
results.history = history;
results.total_trials = trial_counter;
results.completed_iterations = numel(history);
end

function cfg = default_configuration()
cfg.iterations = 50;                  % includes baseline and validation batches
cfg.baseline_trials = 8;
cfg.validation_trials = 8;
cfg.candidates_per_iteration = 24;
cfg.acceptance_threshold = 0.05;
cfg.initial_values = [repmat([0.5, -0.5], 1, 5), zeros(1, 9)];
cfg.number_tickets = 10000;
cfg.initial_mu = 0;
cfg.initial_sigma = 0.05;
cfg.sigma_bounds = [0.01, 0.2];
cfg.learning_rate = 200;
cfg.learning_rate_mu = 0.01;
cfg.learning_rate_sigma = 0.005;
cfg.max_proposal_draws = 1000000;
cfg.score_weights = [5, 5, 5, 2];      % error, covariance, rest movement, variance
cfg.score_divisor = 15;               % intentional; not the sum of the weights
cfg.score_epsilon = 1e-6;
cfg.tremor_metric = 'rest';           % 'rest' or 'movement' in the weighted score
cfg.require_rest_tremor_improvement = true;
cfg.minimum_rest_tremor_improvement = 0.025;
cfg.failed_gate_score = -0.05;
cfg.dbs_frequency_hz = 130;           % scalar, or 9 values in parameter order
cfg.dcn_activation_scale = 0.08;      % explicit analogue of BG/VL activation
cfg.random_seed = 1;
cfg.use_parallel = false;
cfg.verbose = true;
cfg.samples = 5000;
cfg.sampling_frequency = 250;
cfg.cortical_circuits = 1000;
cfg.cerebellar_circuits = 100;
cfg.basal_ganglia_circuits = 10;
cfg.neuron_firing = 0.8;
cfg.initial_position = 0.45;
cfg.muscle_gain = 1.5;                % divided by cortical_circuits, then 10
cfg.resistance = 0;
cfg.desired_amplitude = 0.3;
cfg.desired_period = 10;
cfg.desired_centre = 0.5;
cfg.motor_drive = 0.4;
cfg.movement_start_fraction = 0.4;
cfg.analysis_start_fraction = 0.5;
cfg.variance_window = 100;
cfg.acth_release_rate = 2;
cfg.movement_dopamine = 2;
cfg.parkinsonian_level = 0.4;
cfg.bg_cerebellum_connectivity = 0.4;
cfg.stn_brainstem_connectivity = 1.8;
cfg.cerebellum_VL_weight = 0.036;
cfg.cm_history_samples = 50;
cfg.hyperdirect_firing = 0;           % retains the original driver's absent input
cfg.bg_weights = struct();
cfg.cortex_weights = struct();
cfg.spinal_weights = struct();
cfg.trial_function = [];             % test hook: metrics = f(values, cfg)
end

function cfg = validate_configuration(cfg)
positive_integers = {'iterations', 'baseline_trials', 'validation_trials', ...
    'candidates_per_iteration', 'number_tickets', 'max_proposal_draws', ...
    'samples', 'cortical_circuits', 'cerebellar_circuits', ...
    'basal_ganglia_circuits', 'variance_window', 'cm_history_samples'};
for k = 1:numel(positive_integers)
    name = positive_integers{k};
    validateattributes(cfg.(name), {'numeric'}, {'scalar', 'real', 'finite', ...
        'integer', 'positive'}, mfilename, name);
end
positive_scalars = {'initial_sigma', 'score_divisor', 'score_epsilon', ...
    'sampling_frequency', 'desired_period'};
for k = 1:numel(positive_scalars)
    name = positive_scalars{k};
    validateattributes(cfg.(name), {'numeric'}, {'scalar', 'real', 'finite', 'positive'}, ...
        mfilename, name);
end
nonnegative_scalars = {'learning_rate', 'learning_rate_mu', 'learning_rate_sigma', ...
    'dcn_activation_scale', 'acceptance_threshold'};
for k = 1:numel(nonnegative_scalars)
    name = nonnegative_scalars{k};
    validateattributes(cfg.(name), {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}, ...
        mfilename, name);
end
finite_scalars = {'initial_mu', 'failed_gate_score', 'muscle_gain', 'resistance', ...
    'desired_amplitude', 'desired_centre', 'acth_release_rate', 'movement_dopamine', ...
    'parkinsonian_level', 'bg_cerebellum_connectivity', 'stn_brainstem_connectivity', ...
    'cerebellum_VL_weight'};
for k = 1:numel(finite_scalars)
    name = finite_scalars{k};
    validateattributes(cfg.(name), {'numeric'}, {'scalar', 'real', 'finite'}, mfilename, name);
end
probabilities = {'initial_position', 'neuron_firing', 'motor_drive', ...
    'hyperdirect_firing', 'minimum_rest_tremor_improvement'};
for k = 1:numel(probabilities)
    name = probabilities{k};
    validateattributes(cfg.(name), {'numeric'}, ...
        {'scalar', 'real', 'finite', '>=', 0, '<=', 1}, mfilename, name);
end
validateattributes(cfg.random_seed, {'numeric'}, ...
    {'scalar', 'integer', '>=', 0, '<', 2^32}, mfilename, 'random_seed');
validateattributes(cfg.initial_values, {'numeric'}, ...
    {'vector', 'numel', 19, 'real', 'finite'}, mfilename, 'initial_values');
cfg.initial_values = reshape(cfg.initial_values, 1, []);
if any(cfg.initial_values(11:19) < 0)
    error('OGL:NegativeInitialDBS', 'Initial DBS amplitudes must be nonnegative.');
end
validateattributes(cfg.sigma_bounds, {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'finite', 'positive'});
cfg.sigma_bounds = reshape(cfg.sigma_bounds, 1, []);
if cfg.sigma_bounds(1) > cfg.sigma_bounds(2) || ...
        cfg.initial_sigma < cfg.sigma_bounds(1) || cfg.initial_sigma > cfg.sigma_bounds(2)
    error('OGL:SigmaBounds', 'Initial sigma must lie within ordered sigma_bounds.');
end
validateattributes(cfg.score_weights, {'numeric'}, ...
    {'vector', 'numel', 4, 'real', 'finite', 'nonnegative'});
cfg.score_weights = reshape(cfg.score_weights, 1, []);
cfg.tremor_metric = validatestring(cfg.tremor_metric, {'rest', 'movement'});
validateattributes(cfg.dbs_frequency_hz, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonnegative'});
if isscalar(cfg.dbs_frequency_hz)
    cfg.dbs_frequency_hz = repmat(cfg.dbs_frequency_hz, 1, 9);
elseif numel(cfg.dbs_frequency_hz) == 9
    cfg.dbs_frequency_hz = reshape(cfg.dbs_frequency_hz, 1, []);
else
    error('OGL:FrequencyCount', 'dbs_frequency_hz must contain one or nine values.');
end
validateattributes(cfg.movement_start_fraction, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', 0, '<', 1});
validateattributes(cfg.analysis_start_fraction, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>', cfg.movement_start_fraction, '<', 1});
if floor(cfg.samples*cfg.movement_start_fraction) < 2 || ...
        cfg.samples-floor(cfg.samples*cfg.analysis_start_fraction)+1 < 2
    error('OGL:AnalysisWindow', 'Rest and movement windows need at least two samples.');
end
flags = {'use_parallel', 'verbose', 'require_rest_tremor_improvement'};
for k = 1:numel(flags)
    validateattributes(cfg.(flags{k}), {'logical', 'numeric'}, ...
        {'scalar', 'binary'}, mfilename, flags{k});
end
validateattributes(cfg.bg_weights, {'struct'}, {'scalar'});
validateattributes(cfg.cortex_weights, {'struct'}, {'scalar'});
validateattributes(cfg.spinal_weights, {'struct'}, {'scalar'});
if ~isempty(cfg.trial_function) && ~isa(cfg.trial_function, 'function_handle')
    error('OGL:TrialFunction', 'trial_function must be empty or a function handle.');
end
if cfg.use_parallel && ~license('test', 'Distrib_Computing_Toolbox')
    error('OGL:ParallelToolbox', 'Set use_parallel=false or enable Parallel Computing Toolbox.');
end
if isempty(cfg.trial_function)
    dependencies = {'basal_ganglia_fcn_dbs', 'M1_fcn_pts_dbs_corrected', ...
        'hbm_vnr_dbs', 'spinal_circuits_fcn_pts', 'put_in_range', ...
        'direction_dbs_logistic_curve'};
    for k = 1:numel(dependencies)
        if isempty(which(dependencies{k}))
            error('OGL:MissingDependency', 'Add %s.m to the MATLAB path.', dependencies{k});
        end
    end
end
end

function [metrics, seeds, trial_counter] = run_batch(values, cfg, trial_counter)
number_trials = size(values, 1);
% Distinct deterministic seeds; no worker inherits another trial's state.
seeds = mod(cfg.random_seed + trial_counter + (1:number_trials), 2^32);
trial_counter = trial_counter + number_trials;
metrics = zeros(number_trials, 5);
if cfg.use_parallel
    parfor j = 1:number_trials
        metrics(j, :) = run_trial(values(j, :), cfg, seeds(j));
    end
else
    for j = 1:number_trials
        metrics(j, :) = run_trial(values(j, :), cfg, seeds(j));
    end
end
end

function metrics = run_trial(values, cfg, seed)
previous_rng = rng;
cleanup_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(seed, 'twister');
if isempty(cfg.trial_function)
    metrics = simulate_trial(values, cfg);
else
    metrics = cfg.trial_function(values, cfg);
end
validateattributes(metrics, {'numeric'}, {'vector', 'numel', 5, 'real', 'finite'});
metrics = reshape(metrics, 1, []);
if any(metrics([1, 3, 4, 5]) < 0)
    error('OGL:InvalidMetrics', 'Error, movement magnitude and variances must be nonnegative.');
end
end

function scores = objective_scores(metrics, reference, cfg)
if strcmp(cfg.tremor_metric, 'rest')
    variance_column = 5;
else
    variance_column = 4;
end
columns = [1, 2, 3, variance_column];
% abs(reference) preserves the correct direction even for negative covariance.
denominators = max(abs(reference(columns)), cfg.score_epsilon);
improvements = (metrics(:, columns)-reference(columns))./denominators;
improvements = improvements.*[-1, 1, -1, -1];
scores = (improvements*cfg.score_weights(:))/cfg.score_divisor;
if cfg.require_rest_tremor_improvement
    failed = metrics(:, 5) > ...
        (1-cfg.minimum_rest_tremor_improvement)*reference(5);
    scores(failed) = cfg.failed_gate_score;
end
if any(~isfinite(scores))
    error('OGL:NonfiniteScore', 'Objective score overflow; check weights and metric scaling.');
end
end

function record = make_record(iteration, phase, values, metrics, relative, global_scores, seeds, reference)
record.iteration = iteration;
record.phase = phase;
record.parameters = values;
record.metrics = metrics;
record.relative_scores = relative;
record.global_scores = global_scores;
record.trial_seeds = seeds;
record.reference_before = reference;
record.reference_after = reference;
record.changed_parameters = [];
record.changes = [];
record.proposal_draws = [];
record.ticket_updates = [];
record.selected_candidate = 0;
record.accepted = false;
record.beyond_iteration_budget = false;
end

function metrics = simulate_trial(p, cfg)
% One original motor program, with the new function signatures wired explicitly.
N = cfg.samples;
circuits = cfg.cortical_circuits;
time = (0:N)/cfg.sampling_frequency;
desired = cfg.desired_centre + cfg.desired_amplitude*sin(2*pi*time/cfg.desired_period);
movement_start = floor(cfg.movement_start_fraction*N);
movement_analysis_start = floor(cfg.analysis_start_fraction*N);
dopamine = zeros(1, N);
dopamine(movement_start:end) = cfg.movement_dopamine;
muscleposition = ones(1, N+1)*cfg.initial_position;
velocity = zeros(1, N);
golgi_flexor = zeros(1, cfg.cerebellar_circuits);
golgi_extensor = zeros(1, cfg.cerebellar_circuits);
cerebellar_agonist = zeros(1, cfg.cerebellar_circuits);
cerebellar_antagonist = zeros(1, cfg.cerebellar_circuits);
L5_flexor_prev = zeros(1, circuits);
L5_extensor_prev = zeros(1, circuits);
CM_history = zeros(1, cfg.cm_history_samples);

bg = cfg.bg_weights;
bg.keepdefaults = true;
bg.neuron_firing = cfg.neuron_firing;
bg.glutamine_mod = p(1);
bg.gaba_mod = p(2);
bg.dbs_striatum_a = p(11); bg.dbs_striatum_f = cfg.dbs_frequency_hz(1);
bg.dbs_gpe_a = p(12);      bg.dbs_gpe_f = cfg.dbs_frequency_hz(2);
bg.dbs_stn_a = p(13);      bg.dbs_stn_f = cfg.dbs_frequency_hz(3);
bg.dbs_gpi_a = p(14);      bg.dbs_gpi_f = cfg.dbs_frequency_hz(4);
pathologies = struct('keepdefaults', true);

cortex = cfg.cortex_weights;
cortex.keepdefaults = true;
cortex.neuron_fire = cfg.neuron_firing;
cortex.glutamine_mod = p(3);
cortex.gaba_mod = p(4);
cortex.cerebellum_VL_weight = cfg.cerebellum_VL_weight;
cortex.dbs_thalamus_a = p(15); cortex.dbs_thalamus_f = cfg.dbs_frequency_hz(5);
cortex.dbs_L4_a = p(16);       cortex.dbs_L4_f = cfg.dbs_frequency_hz(6);
cortex.dbs_L23_a = p(17);      cortex.dbs_L23_f = cfg.dbs_frequency_hz(7);
cortex.dbs_L5_a = p(18);       cortex.dbs_L5_f = cfg.dbs_frequency_hz(8);

spinal = cfg.spinal_weights;
spinal.keepdefaults = true;
spinal.mod_excitation = p(5);
spinal.mod_inhibition = p(6);
dcn = struct('dbs_DCN_a', p(19), 'dbs_DCN_f', cfg.dbs_frequency_hz(9), ...
    'dbs_DCN_activation_scale', cfg.dcn_activation_scale);
gain = cfg.muscle_gain/circuits;

for i = 1:N
    % Read the oldest CM history slot before replacing it with this sample.
    % The full history is also supplied for the BG cholinergic modifier.
    current_index = mod(i-1, cfg.cm_history_samples)+1;
    bg_output = basal_ganglia_fcn_dbs(cfg.motor_drive, cfg.acth_release_rate, ...
        dopamine(i)*cfg.parkinsonian_level, bg, pathologies, ...
        cfg.basal_ganglia_circuits, cfg.hyperdirect_firing, CM_history, current_index);

    desired_input = 0.04 + 0.32*desired(i);
    muscle_input = 0.04 + 0.32*muscleposition(i);
    desired_input_rev = 0.36-desired_input;
    muscle_input_rev = 0.36-muscle_input;
    desired_input = put_in_range(desired_input + p(10)*cfg.bg_cerebellum_connectivity* ...
        rand*bg_output.gpi_dep_probability, 0, 0.36);
    muscle_input = put_in_range(muscle_input + p(9)*cfg.stn_brainstem_connectivity* ...
        rand*bg_output.stn_dep_probability + p(10)*0.1*cfg.bg_cerebellum_connectivity* ...
        rand*bg_output.gpi_dep_probability, 0, 0.36);
    desired_input_rev = put_in_range(desired_input_rev + p(10)*cfg.bg_cerebellum_connectivity* ...
        rand*bg_output.gpi_dep_probability, 0, 0.36);
    muscle_input_rev = put_in_range(muscle_input_rev + p(9)*cfg.stn_brainstem_connectivity* ...
        rand*bg_output.stn_dep_probability + p(10)*0.1*cfg.bg_cerebellum_connectivity* ...
        rand*bg_output.gpi_dep_probability, 0, 0.36);
    if dopamine(i) == 0
        desired_input = 0;
        desired_input_rev = 0;
    end
    for y = 1:cfg.cerebellar_circuits
        [cerebellar_agonist(y), golgi_flexor(y)] = hbm_vnr_dbs( ...
            desired_input, muscle_input_rev, golgi_flexor(y), p(8), p(7), cfg.neuron_firing, dcn);
        [cerebellar_antagonist(y), golgi_extensor(y)] = hbm_vnr_dbs( ...
            desired_input_rev, muscle_input, golgi_extensor(y), p(8), p(7), cfg.neuron_firing, dcn);
    end

    % Preserve the original shared BG motor program: the same GPi population
    % supplies agonist and antagonist inputs required by the new M1 interface.
    [L5_flexor, L5_extensor, ~, ~, CM] = M1_fcn_pts_dbs_corrected( ...
        cerebellar_agonist, cerebellar_antagonist, bg_output.gpi_firings, ...
        bg_output.gpi_firings, cfg.motor_drive, circuits, ...
        L5_flexor_prev, L5_extensor_prev, cortex);
    L5_flexor_prev = L5_flexor;
    L5_extensor_prev = L5_extensor;
    CM_history(current_index) = (CM.CM_agonist + CM.CM_antagonist)/2;

    % Pass the structure containing BOTH spinal modifiers.
    spinal_output = spinal_circuits_fcn_pts(L5_flexor, L5_extensor, muscleposition(i), spinal);
    velocity(i) = (sum(spinal_output.flexorpower(:)) - ...
        sum(spinal_output.extensorpower(:)))*gain/10;
    muscleposition(i+1) = put_in_range( ...
        muscleposition(i)+velocity(i)-cfg.resistance, 0, 1);
end

% Preserve original analysis indices, including the transition sample at the
% end of the resting window. The final N+1 position is outside these metrics.
rest_indices = 1:movement_start;
move_indices = movement_analysis_start:N;
joint = muscleposition(1:N);
covariance_matrix = cov(joint(move_indices), desired(move_indices));
metrics = [mean(abs(desired(move_indices)-joint(move_indices))), ...
    covariance_matrix(1, 2), mean(abs(velocity(rest_indices))), ...
    mean(movvar(joint(move_indices), cfg.variance_window)), ...
    mean(movvar(joint(rest_indices), cfg.variance_window))];
end
