# Corrected objective-guided MATLAB training

This bundle updates the supplied training script, connects the supplied DBS
interfaces, and adds DBS to the DCN in the supplied cerebellar circuit. Original
uploaded files are preserved. The new entry point is
`for_network_training_corrected.m`.

## Run

Extract the bundle and put its top-level folder on the MATLAB path, together
with the folder containing your existing project helpers:

```matlab
addpath('/path/to/objective_guided_dbs');
addpath('/path/to/your/project/helpers');

cfg = for_network_training_corrected('defaults');
cfg.use_parallel = true;  % optional; requires Parallel Computing Toolbox
results = for_network_training_corrected(cfg);
disp(results.parameter_table);
save('objective_guided_results.mat', 'results');
```

The default run uses the original 5,000 samples and 1,000 cortical circuits per
trial, so it is computationally expensive. Serial execution is the default.
Use the regression tests below before a full model run.

Three original project dependencies are still required and were not supplied:
`spinal_circuits_fcn_pts.m`, `direction_dbs_logistic_curve.m`, and
`put_in_range.m`. The driver checks that they can be found before running.
The scientific logistic curve and spinal model have not been reconstructed.

Add only the bundle's top-level folder for a model run; do not use
`addpath(genpath(...))`. The `tests/fixtures` folder contains synthetic test
doubles. They deliberately raise an error if invoked outside the tests.

## Included model and search files

| File | Purpose |
| --- | --- |
| `for_network_training_corrected.m` | Main function; defaults, optimization, model calls, scoring and history. |
| `hbm_vnr_dbs.m` | Supplied cerebellar equations with DBS confined to the DCN. |
| `M1_fcn_pts_dbs_corrected.m` | Supplied cortical function with the CM DBS multiplication parentheses corrected. |
| `basal_ganglia_fcn_dbs.m` | Exact supplied basal-ganglia file, with its canonical MATLAB filename. |
| `propose_parameter_change.m` | Normal increments with rejection sampling for DBS amplitudes. |
| `update_distributions_learning_corrected.m` | Score-proportional ticket, mean and spread updates. |

The corrected M1 function multiplies the entire CM input sum by the thalamic
DBS factor. In the supplied function, that factor multiplied only the final
opposite-side layer-5 input. The remaining equations are preserved, including
the **additive** L4/L23/L5 DBS terms in your supplied M1 function. Those cortical
terms differ from the multiplicative BG/VL/DCN input rule; they have not been
silently converted to a different formulation.

## DCN extension

The first six inputs retain their original meaning. A seventh, optional
structure supplies the DBS settings:

```matlab
dbs = struct('dbs_DCN_f', 130, 'dbs_DCN_a', 0.3, ...
             'dbs_DCN_activation_scale', 0.08);
[dcn, golgi, purkinje, granule, interneuron] = hbm_vnr_dbs( ...
    inferior_olive_probability, pontine_probability, previous_golgi, ...
    inhibition_modifier, excitation_modifier, 0.8, dbs);
```

For positive amplitude `a`, frequency `f`, and summed DCN input `I`, the update is:

```matlab
Q = 1 + 2.5*a*(rand + direction_dbs_logistic_curve(f, 0) - 0.5);
synaptic_firing = I*Q > (1-neuron_firing);
independent_firing = rand*a*activation_scale >= direction_dbs_logistic_curve(f, 1);
dcn_firing = synaptic_firing || independent_firing;
```

The two random draws are independent. `activation_scale = 0.08` is an explicit,
adjustable assumption matching most supplied BG targets and thalamus; it was
not fitted specifically for DCN. The original Golgi, granule, interneuron and
Purkinje equations are preserved. At zero amplitude the first five outputs and
random-number consumption match the original cerebellar function exactly.
The original declared sixth output was never assigned; the revised function
returns `[]` for that unused output instead of inventing a new quantity.

## Corrections to learning

- Eight baseline trials use identical initial parameters and distinct seeds.
- Each exploratory batch has 24 single-parameter changes from the retained state.
- All candidates, including the first, update the learned distributions.
- Ticket replacements increase with score magnitude: `floor(200*abs(score))`,
  with score clipped to [-1, 1] for learning only. Mean/spread updates remain
  proportional to the signed score; the default spread limits are 0.01 to 0.2.
- Eight validation trials use the identical selected candidate. The retained
  parameters change only if their mean score exceeds 0.05. Rejection retains
  the previous validated state. A pending final validation can add one batch
  beyond the nominal 50-iteration budget.
- Higher covariance improves its reference; lower error, resting movement and
  both variance metrics improve theirs. Variance is compared as a scalar mean.
  Trade-offs between metrics remain permitted through the combined score.
- Scores use the intentional divisor 15. Reference denominators are protected
  against zero; covariance uses its absolute reference magnitude in the
  denominator so negative covariance does not reverse the improvement direction.
- For the default resting-tremor search, a trial failing the 2.5% resting-variance
  improvement gate receives -0.05. This gate is applied to each trial before
  validation scores are averaged, as in the original algorithm.
- Both spinal modifiers are assigned to the correct fields and their structure
  is passed to the spinal function.
- DBS frequency/amplitude fields now use the names expected by the new modules.
  The old post-threshold DCN inhibition/excitation gates have been removed.
- Batch histories use their actual sizes, so old scores cannot remain in an
  eight-trial validation batch after a 24-trial exploration batch.

The normal distributions describe **increments**, preserving your original
algorithm. A DBS proposal is redrawn until `current_amplitude + increment > 0`.
A negative increment is valid if the resulting amplitude remains positive.
Initial zeros mean no stimulation; untouched targets can remain zero. No
absolute-value transformation, clipping or upper amplitude bound is applied.
A configurable retry limit raises an error instead of hanging indefinitely
if a learned distribution assigns negligible probability to admissible values.

## Settings commonly changed

Use `cfg = for_network_training_corrected('defaults')` to inspect all settings.

| Setting | Default | Meaning |
| --- | --- | --- |
| `iterations` | 50 | Batch budget, including baseline and validation. |
| `baseline_trials`, `validation_trials` | 8, 8 | Repeated simulations at identical parameter values. |
| `candidates_per_iteration` | 24 | Single-parameter proposals per exploratory batch. |
| `acceptance_threshold` | 0.05 | Strict score threshold for selection and validation. |
| `score_weights` | [5 5 5 2] | Error, covariance, resting movement, selected variance. |
| `score_divisor` | 15 | Score scaling; intentionally independent of weight sum. |
| `tremor_metric` | 'rest' | Use 'movement' to score moving variance instead. |
| `require_rest_tremor_improvement` | true | Enable the resting-variance gate. |
| `minimum_rest_tremor_improvement` | 0.025 | Required relative resting-variance improvement. |
| `initial_mu`, `initial_sigma` | 0, 0.05 | Initial normal distribution of increments. |
| `sigma_bounds` | [0.01 0.2] | Lower/upper spread limits. |
| `learning_rate`, `learning_rate_mu`, `learning_rate_sigma` | 200, 0.01, 0.005 | Ticket, mean and spread adaptation. |
| `dbs_frequency_hz` | 130 | One shared frequency or nine target-specific frequencies. |
| `dcn_activation_scale` | 0.08 | DCN independent-activation coefficient. |
| `random_seed` | 1 | Repeatable proposals and distinct seeded trial simulations. |
| `use_parallel` | false | Optional parallel evaluation of each trial batch. |

Parameter order is unchanged: BG excitation/inhibition, cortex
excitation/inhibition, spinal excitation/inhibition, cerebellum
excitation/inhibition, brainstem excitation/inhibition, followed by DBS amplitudes
for striatum, GPe, STN, GPi, thalamus, L4, L2/3, L5 and DCN. The nine frequencies
use that same DBS target order. Initial modifiers are [0.5 -0.5] for each pair;
all nine DBS amplitudes initially equal zero.

For example, to emphasize movement and disable the resting-tremor gate:

```matlab
cfg.score_weights = [5, 5, 1, 2];
cfg.tremor_metric = 'movement';
cfg.require_rest_tremor_improvement = false;
```

## Feedback and timing choices

The original driver has one BG motor program. Its GPi output is therefore
passed to both GPi arguments required by the new cortical interface. Previous
layer-5 activity is initialized to zero and updated every sample. CM agonist and
antagonist probabilities are averaged into a 50-sample history; the oldest
entry is read by the BG input sampler before replacement, while the history's
mean supplies its cholinergic modifier. `hyperdirect_firing` defaults to zero,
preserving the absence of this input in the original driver. These integration
choices are explicit because the old training script did not contain these
new feedback inputs.

The original metric intervals are retained: samples 1:2000 for rest and
2500:5000 for movement. Dopamine begins at sample 2000, so the resting window
includes that transition sample, just as in the supplied script. Changing
`sampling_frequency` alone does not rescale probabilities in the supplied BG
module, whose intrinsic rates use 250 Hz.

## Verification

All 12 MATLAB files passed the MISS_HIT syntax/lint check in this environment.
MATLAB and Octave were unavailable, so the MATLAB tests and full model have
**not** been executed here. The missing original spinal and curve dependencies
also prevent a scientific end-to-end validation in this environment.

Run the included 20 regression tests in MATLAB:

```matlab
test_results = runtests('tests/test_objective_guided_learning.m');
assert(all([test_results.Passed]));
```

Tests cover baseline/validation identity, rejection, final pending validation,
metric references and trade-offs, all-candidate updates, proportional learning,
amplitude rejection/decreases, score guards, reproducibility, exact no-DBS
cerebellar regression, unchanged upstream cells, both DCN DBS components, and
module/spinal interface wiring with synthetic test doubles. Follow these with
a small run using the actual project helpers before the full search.

The results include the retained parameter table, initial/final metrics, best
per-metric references, learned distributions, selection probabilities, and full
per-batch settings, seeds, scores, candidate changes and acceptance decisions.
