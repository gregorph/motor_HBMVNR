# Sensitivity analysis: `fig_constraints_check.m`

This script evaluates motor trajectories while varying intrinsic pacemaker parameters or the transformation of intended-position and proprioceptive inputs. The simulation and analysis blocks below generate and process new trajectories for the two active input-transformation sweeps.

## Dependencies and execution blocks

Use the environment specified in the [main README](README.md). The main simulation depends on `basal_ganglia_fcn_pts_rtl`, `M1_fcn_pts_rtl`, `hbm_vnr`, `spinal_circuits_fcn_pts` and `put_in_range`. Wavelet analysis uses `cwt`; spectral-peak estimation uses `findpeaks`.

| Main block | Purpose |
| --- | --- |
| `%% calculate motor trajectory changing pacemaker freq` | Generate trajectories across the selected parameter sweeps |
| `%% calculate freqs from motor trajectory` | Calculate spectral and trajectory measures |
| `%% test` | Independent exploratory simulation |
| `%% figure generator` | Assemble a composite plot across sweep families |

Execute the main blocks separately. Several blocks clear the workspace, so save generated outputs before proceeding. Section labels inside simulation loops are part of their parent block.

## Simulation protocol

The active trajectory configuration uses 1,000 thalamocortical circuits, 10,000 time steps, 250 Hz sampling and 48 repeats per parameter setting. Each simulation lasts 40 s. The target is:

```matlab
desired = 0.5 + 0.3*sin(2*pi*t/3);
```

Movement drive is 0.4 at indices 2501:7500 and 0.02 elsewhere. Dopamine is 0.8 during those movement indices and 0.75 elsewhere. This gives approximately 10 s rest, 20 s movement and 10 s rest.

The source defines six parameter branches:

| Condition index | Parameter varied | Values / definition |
| --- | --- | --- |
| 1 | STN pacemaker parameter | Derived from `freqrange_check` |
| 2 | GPe pacemaker parameter | Derived from `freqrange_check` |
| 3 | GPi pacemaker parameter | Derived from `freqrange_check` |
| 4 | Input transformation shift | `vertical_shift = -0.02:0.02:0.6` |
| 5 | Input transformation coefficient | `max_freqprop = 0.34:0.02:0.96` |
| 6 | Input transformation lower coefficient | `min_freqprop = -0.3:0.02:0.32`, with shift 0.6 |

The active loop selects conditions **4 and 5**, with 32 settings and 48 repeats each: 3,072 simulations. Condition 4 uses `min_freqprop=0.04` and `max_freqprop=0.36`; condition 5 uses `min_freqprop=0.04` and zero shift. Both intended-position and proprioceptive inputs are transformed, then bounded to `[0,0.36]`. Intended-position inputs are suppressed during rest.

The structure is indexed as `results(condition,setting,repeat)`. Preallocated entries outside the selected conditions, setting range and repeats are not observations.

## Generate trajectories

Open the script in the MATLAB Editor. Run the complete first main block, stopping before `%% calculate freqs from motor trajectory`. Save the generated trajectories and configuration:

```matlab
save('constraints_trajectories.mat', 'results', 'num_repetitions', ...
    'freqrange_check', 'repetitions', 'fs', 'desiredposition', '-v7.3');
```

Each populated `results(condition,setting,repeat).muscleposition_keep` contains 10,000 position samples. Spectral analysis is performed in the next block.

## Analyse the generated trajectories

In a working analysis copy, set the `matfile` assignment in the second main block to:

```matlab
matObj = matfile(fullfile(pwd, 'constraints_trajectories.mat'));
```

Keep this assignment after that block's workspace initialisation. Run the complete block, stopping before `%% test`, then save:

```matlab
save('constraints_summary.mat', 'results_avg', ...
    'repetitions', 'fs', 'desiredposition', '-v7.3');
```

The analysis loops match the active simulation selection: conditions 4:5, settings 1:32 and repeats 1:48. If the simulation selection is changed, update the analysis selection and repeat-count denominators consistently.

## Measures and outputs

| Quantity | Definition |
| --- | --- |
| `error` | Mean absolute target-position error at indices 2500:7500 |
| `peak_frequency_rest` | Magnitude-weighted centroid of the resting spectrum |
| `tallest_frequency_rest` | Mean across repeats of the tallest detected resting spectral peak |
| `power_rest` | Summed wavelet magnitude over indices 1:2500 and 7500:10000 |
| `power_movement` | Summed wavelet magnitude over indices 2500:7500 |
| `tremor_mov_ratio` | Mean across repeats of the rest/movement magnitude-sum ratio |

The transform uses frequency limits 0.5–20 Hz, 20 voices per octave and time-bandwidth 120. Complex coefficients are interpolated onto `0.5:0.25:20` using `makima` before taking their magnitude. Resting frequency estimates use the mean spectrum over indices 1000:2000 and 8000:9000. Peak detection uses a 10-point moving mean and minimum peak width 0.5 Hz.

The `power_*` names denote wavelet-magnitude sums, rather than squared power or PSD. The summarised fields have corresponding `_std` and `_se` fields. For the specified 48-repeat workflow, SE is calculated using `sqrt(48)` as the denominator. The mean spectral-peak field requires a detected peak for each included realisation.

## Plot the active sweeps

```matlab
D = load('constraints_summary.mat', 'results_avg');
S = D.results_avg;
x4 = -0.02:0.02:0.6;
x5 = 0.34:0.02:0.96;

figure;
subplot(1,2,1);
errorbar(x4, [S(4,1:32).tallest_frequency_rest], ...
    1.96*[S(4,1:32).tallest_frequency_rest_se], 'o-');
xlabel('Input shift (model coefficient)');
ylabel('Resting peak frequency (Hz)');

subplot(1,2,2);
errorbar(x5, [S(5,1:32).tallest_frequency_rest], ...
    1.96*[S(5,1:32).tallest_frequency_rest_se], 'o-');
xlabel('Input range coefficient');
ylabel('Resting peak frequency (Hz)');

savefig(gcf, 'constraints_input_sweeps.fig');
print(gcf, 'constraints_input_sweeps.png', '-dpng', '-r300');
```

The x-axes use the coefficients actually varied in the two active sweeps. Error bars represent mean ± 1.96 × SE.

The final `%% figure generator` block is a separate composite layout combining input-transformation and pacemaker sweeps. Its panels use condition rows 4/5 and 1/3 respectively. Use the condition-specific parameter vectors when assembling outputs from different sweep lengths. The two-panel workflow above covers the currently selected 4/5 simulations.

## Exploratory test and reproducibility

The `%% test` block runs a separate 500-circuit simulation with `min_freqprop=-0.15`, `max_freqprop=0.36` and `vertical_shift=0.6`, then plots a normalised wavelet-magnitude spectrum. It uses `basal_ganglia_fcn_dbs` and `direction_dbs_logistic_curve`. It is independent of the trajectory-sweep analysis and clears its own workspace.

Archive the generated trajectories, summaries, source files and parameter selection together. Preserve the condition/setting/repeat mapping when comparing sweeps. Parallel worker random streams are used without an explicit saved seed/substream per repeat.
