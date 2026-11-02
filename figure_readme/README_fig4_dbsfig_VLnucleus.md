# DBS simulations: `fig4_dbsfig_VLnucleus.m`

This script defines basal-ganglia response sweeps and motor simulations under deep brain stimulation. The motor block supports a fixed-amplitude switch-on protocol and a stepped-amplitude protocol. Newly generated outputs can be used for trajectory plots and statistical analysis.

## Dependencies

Use the environment specified in the [main README](README.md), together with:

- `basal_ganglia_fcn_dbs.m`
- `M1_fcn_pts_dbs.m` and `M1_fcn_pts_rtl.m`
- `hbm_vnr.m`
- `spinal_circuits_fcn_pts.m`
- `put_in_range.m`
- `direction_dbs_logistic_curve.m`

The DBS helper supplies frequency-dependent terms for the synaptic-modulation and added-activation components of stimulation.

## Execution structure

| Block | Purpose |
| --- | --- |
| `%% part 1` | BG responses across dopamine, motor drive and DBS amplitude |
| `%% part 2` | Repeated motor simulations under a selected DBS protocol |
| `%% ploting` | Composite figure assembly from the corresponding simulation outputs |

Run the main blocks separately in the MATLAB Editor. Begin each simulation stage with a clean workspace. Internal section labels within loops are not independent execution units. Save a stage's outputs before clearing its workspace or beginning figure assembly.

For independent target configurations, initialise the parameter structure at each sweep/condition setup with:

```matlab
weight_params = struct('keepdefaults', true);
```

For a working simulation copy, use this initialisation in place of the single-field `weight_params.keepdefaults = true;` assignment within the Part 1 `dbs_stim` loop and the Part 2 `condition_loop`. Apply the subsequent condition-specific assignments after initialisation.

## BG response sweep

Part 1 defines the following grids:

| Parameter | Values |
| --- | --- |
| Baseline dopamine | `0:0.025:3` |
| Motor drive | `0.01:0.005:0.4` |
| DBS amplitude | `0:0.025:1` |
| Stimulation frequency | 130 Hz |
| BG circuits per call | 10,000 |

The programmed condition indices are 2 for the unstimulated reference, 3 for STN stimulation and 4 for GPi stimulation. The dopamine input is half the baseline sweep value. The outputs `condition_checked(k).gpi_dp` and `.stn_dp` have amplitude × dopamine × motor-drive layout and represent firing probability multiplied by 250.

After running the complete Part 1 block, save:

```matlab
save('fig4_bg_sweep.mat', 'condition_checked', 'dbs_rate_simulate', ...
    'dopamine_levels', 'motor_drive_frequency', 'dop_acth', '-v7.3');
```

## Fixed-amplitude switch-on simulation

Run the complete Part 2 block with `repetitions = 10000`, stopping before `%% ploting`. This stage can be executed independently of the BG sweep.

The protocol uses 250 Hz sampling, 1,000 thalamocortical circuits and 48 repeats per condition. Each simulation lasts 40 s. Odd repeats represent rest and even repeats represent movement, giving 24 runs of each state per condition.

| Condition index | Configuration | Amplitude after onset | Frequency |
| --- | --- | --- | --- |
| 1 | Unstimulated reference | 0 | — |
| 2 | STN DBS | 0.8 | 130 Hz |
| 3 | GPi DBS | 0.9 | 130 Hz |
| 4 | Thalamic DBS | 0.5 | 130 Hz |

All four conditions use the adapted-parkinsonian network profile. The analysis label `Control` denotes the unstimulated reference. The thalamic DBS implementation affects both VL outputs and CM circuitry.

Rest uses motor drive 0.05 and dopamine 0.75. Movement uses motor drive 0.4 and dopamine 0.8. Programme 1 is simulated with target:

```matlab
desired = 0.5 + 0.3*sin(2*pi*0.1*t);
```

Stimulation is enabled from sample index 5000 by `i >= repetitions/2`, corresponding to approximately 20 s. The plotting convention marks onset at 20 s.

After completing the simulation:

```matlab
save('fig4_switchon.mat', '-v7.3');
```

The saved workspace retains the actual target and parameter settings. The principal outputs have condition × repeat × sample layout, normally 4 × 48 × 10000:

| Variable | Meaning |
| --- | --- |
| `fig3b_muscleposition` | Achieved position at each simulation time step |
| `fig3b_error_observed` | Desired minus achieved position |
| `fig3b_velocity` | Signed per-step motor output; multiply by `fs` for units per second |
| `fig3b_thalamus_net` | Flexor minus extensor thalamic activity |
| `fig3b_thalamus_summative` | Summed thalamic activity |

The `fig3b_` variable prefix is retained from the source naming convention.

## Plot a newly generated trajectory

Run this immediately after the switch-on simulation, or load `fig4_switchon.mat` first:

```matlab
t = (0:size(fig3b_muscleposition,3)-1)/250;
condition = 2;                       % STN DBS
restRun = 1;
movementRun = 2;

figure;
subplot(1,2,1);
plot(t, squeeze(fig3b_muscleposition(condition,restRun,:)));
xline(20, '--', 'DBS onset');
xlabel('Time (s)'); ylabel('Joint position'); title('Rest');

subplot(1,2,2);
plot(t, squeeze(fig3b_muscleposition(condition,movementRun,:)));
hold on;
plot(t, desiredposition(1,1:numel(t)), '--');
xline(20, '--', 'DBS onset');
xlabel('Time (s)'); ylabel('Joint position'); title('Movement');
legend('Simulated position', 'Target', 'DBS onset');
```

## Stepped-amplitude simulation

Use a separate working copy of Part 2 and a fresh workspace. Set `repetitions = 27500`, retaining the existing amplitude-array definitions. Replace the fixed-amplitude/onset condition block inside the time loop with:

```matlab
weight_internal = weight_params;
if condition_loop == 2
    weight_internal.dbs_stn_f = 130;
    weight_internal.dbs_stn_a = dbs_stn_inhibition(1,i);
elseif condition_loop == 3
    weight_internal.dbs_gpi_f = 130;
    weight_internal.dbs_gpi_a = dbs_gpi_inhibition(1,i);
elseif condition_loop == 4
    weight_internal.dbs_thalamus_f = 130;
    weight_internal.dbs_thalamus_a = dbs_gpi_inhibition(1,i);
end
```

The protocol lasts 110 s and contains eleven 2,500-sample blocks with amplitudes 0, 0.1, …, 1.0. Network state evolves continuously between blocks. The thalamic condition uses the same amplitude schedule as the GPi condition. Amplitude is a dimensionless model parameter.

Save this simulation separately:

```matlab
save('fig4_amplitude_steps.mat', '-v7.3');
```

## Composite figure and export

In a working copy of the `%% ploting` block, point its three sequential dataset loads to the newly generated outputs, keeping each load at its original position:

| Load position | Generated file |
| --- | --- |
| Before the top-row maps | `fig4_bg_sweep.mat` |
| Before the rest/movement example trajectories | `fig4_switchon.mat` |
| Before the amplitude-response panels | `fig4_amplitude_steps.mat` |

For example, use `load(fullfile(pwd,'fig4_switchon.mat'));` for the switch-on dataset. Preserve its actual `desiredposition` variable when plotting movement traces. The example run indices assume 48 repeats; the amplitude panels use eleven 2,500-sample blocks.

After figure generation:

```matlab
savefig(gcf, 'fig4_dbs.fig');
print(gcf, 'fig4_dbs.png', '-dpng', '-r300');
```

## Before/after statistical analysis

Apply the analysis to the 40-second switch-on dataset:

```matlab
D = load('fig4_switchon.mat', 'fig3b_muscleposition', 'desiredposition');
nTime = size(D.fig3b_muscleposition,3);
R = analyse_fig3b_dbs(D.fig3b_muscleposition, ...
    D.desiredposition(1,1:nTime), 'Fs', 250, 'TargetFrequency', 0.1, ...
    'MaxLag', 500, 'NumPermutations', 19999, ...
    'NumBootstraps', 20000, 'RandomSeed', 1);
save('fig4_statistics.mat', 'R', '-v7.3');
```

The [statistical analysis guide](README_statistics.md) describes paired comparisons and circular phase-latency inference. These fixed before/after windows apply to the switch-on protocol. The stepped-amplitude figure uses a separate moving-standard-deviation tremor measure; the statistical function uses moving variance.

Archive generated outputs, source files and settings together. The simulation script uses parallel worker random streams without a recorded per-realisation seed allocation; the analysis seed controls resampling only.
