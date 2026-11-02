# Figure 2: baseline and parkinsonian simulations

`fig2_48_realisations.m` simulates four network conditions, saves each realisation, calculates state-specific metrics and joint-position spectra, and exports the Figure 2 and metric-summary panels. Installation and a reduced demonstration are described in the [main README](README.md).

## Entry point

```matlab
summary = fig2_48_realisations(outputDir, baseSeed, plotOnly, plotRealisations);
```

| Argument | Default | Description |
| --- | --- | --- |
| `outputDir` | Timestamped directory under the current folder | Destination for a new simulation, or existing generated dataset for replotting |
| `baseSeed` | `20260928` | Integer simulation seed in `[0,2^32-1]` |
| `plotOnly` | `false` | Logical flag; `true` regenerates figures from saved simulation outputs |
| `plotRealisations` | `1` | One run index, three run indices, or `'random'` for illustrative trajectory selection |

The three trajectory rows correspond to Healthy, PD and PD severe, in that order. A scalar selects the same realisation for each row. A three-element vector selects one realisation per row. This choice affects illustrative traces; aggregate panels continue to use all completed realisations.

## Run the simulation

```matlab
addpath(pwd);
outDir = fullfile(pwd, 'fig2_results_01');
summary = fig2_48_realisations(outDir, 20260928, false, 1);
```

The standard configuration runs 48 realisations × 4 conditions. Each simulation lasts 84 s at 250 Hz. Use a new or empty output directory and a local process pool. Each completed realisation is saved before aggregate processing; a progress message identifies the saved run.

The model dependencies are `basal_ganglia_fcn_pts_rtl.m`, `hbm_vnr.m`, `M1_fcn_pts_rtl.m`, `spinal_circuits_fcn_pts.m` and `put_in_range.m`.

## Configuration

Model parameters are defined in `cfg` under `%% Settings`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `cfg.n_realisations` | 48 | Independent repeats |
| `cfg.circuits` | 500 | Circuit count passed to the thalamocortical module |
| BG circuit argument | 10 | Basal-ganglia circuits per call |
| `cfg.cerebellar_circuits` | 100 | Cerebellar circuits per agonist/antagonist pathway |
| `cfg.fs` | 250 Hz | Sampling frequency |
| State durations | 6, 24, 48, 6 s | Rest, single-programme movement, two-programme movement, hyperdirect activation |
| `cfg.sma_rest` | 0.05 | Resting motor drive |
| `cfg.sma_active` | 0.4 | Selected-programme drive |
| `cfg.sma_inactive` | 0.15 | Other enabled-programme drive |
| `cfg.target_period` | 10 s | Target period |
| `cfg.target_centre` | 0.5 | Target centre |
| `cfg.target_amplitude` | 0.15 | Programme 1 amplitude; programme 2 uses 1.5 times this amplitude |
| `cfg.save_full_cwt` | `false` | Retain complete wavelet coefficient arrays |
| `cfg.save_neuron_thalamus` | `false` | Retain neuron-level thalamic activity |
| `cfg.make_figure` | `true` | Export figures after simulation |

Programme 2 has a quarter-cycle target phase offset. It is disabled during 6–30 s, enabled at 30 s, and becomes the selected programme at 54 s. State variables are initialised separately for every condition and realisation.

## Output structure

The output directory contains `manifest.mat`, one `realisation_*.mat` file per run, `summary.mat`, CSV tables, and FIG/PNG images. Each raw file contains `runData`, whose `conditions` cell array follows the condition order below.

| Condition index | Label |
| --- | --- |
| 1 | Healthy |
| 2 | PD |
| 3 | PD mild adaptations |
| 4 | PD severe |

| State index | Interval |
| --- | --- |
| 1 | Rest, 0–6 s |
| 2 | Single-programme movement, 6–30 s |
| 3 | Two-programme movement, 30–78 s |
| 4 | Hyperdirect activation, 78–84 s |

| Field | Array layout |
| --- | --- |
| `summary.metrics.<metric>` | Condition × state × realisation |
| `summary.stats.<metric>.<statistic>` | Condition × state |
| `summary.movement_all.metrics.<metric>` | Condition × realisation, over 6–78 s |
| `summary.latency_single.values_s` | Condition × realisation |
| `summary.joint_spectra.<quantity>` | Condition × state × realisation × frequency |
| `summary.joint_spectra.stats.<quantity>.<statistic>` | Condition × state × frequency |

For the full run, scalar metric arrays are 4 × 4 × 48. Raw position has 21,001 samples, including the final endpoint; per-step metrics use 21,000 time-aligned samples.

```matlab
% PD single-programme error for each realisation:
pdError = reshape(summary.metrics.error(2,2,:), [], 1);

% Mean, SE and valid sample counts by condition/state:
meanError = summary.stats.error.mean;
seError = summary.stats.error.se;
nValid = summary.stats.error.n_valid;

% Frequency-weighted spectrum for PD at rest:
f = summary.joint_spectra.frequency_hz;
p = squeeze(summary.joint_spectra.stats.frequency_weighted_power.mean(2,1,:));
figure; plot(f, p);
xlabel('Frequency (Hz)');
ylabel('Frequency-weighted wavelet power');
```

## Measures

| Metric | Definition |
| --- | --- |
| `error` | Mean absolute same-time target-position error |
| `velocity` | Mean absolute within-state 10-sample moving mean of motor output, multiplied by `fs` |
| `tracking_covariance` | Sample covariance of scheduled target and position |
| `tracking_correlation` | Pearson correlation of scheduled target and position |
| `rigidity_covariance` | Sample covariance of summed unsmoothed flexor and extensor outputs |
| `rigidity_correlation` | Pearson correlation of those flexor and extensor outputs |

Muscle outputs are gain-scaled model signals. A constant signal gives undefined Pearson correlation (`NaN`) and zero covariance. Summaries report valid and missing counts. Figure intervals are pointwise mean ± 1.96 × SE across realisations; SE uses the number of valid values.

Single-programme latency is the peak lag of mean-centred `xcorr(position,target)` over 6–30 s, converted to seconds. Positive latency denotes achieved position lagging the target. Periodic targets can produce more than one correlation peak.

Joint-position spectra use a Morse continuous wavelet transform with gamma 3, time-bandwidth 60 and 10 voices per octave. Original power is the within-state mean of squared coefficient magnitude. Displayed spectra multiply this measure by frequency. These wavelet measures are not calibrated PSD per Hz. An additional `power_coi` output excludes recording-edge effects; the standard figure uses all state samples.

## Figure regeneration

After generating the dataset:

```matlab
% Select one realisation for all illustrative rows:
summary = fig2_48_realisations(outDir, [], true, 7);

% Select individual realisations for Healthy, PD and PD severe:
summary = fig2_48_realisations(outDir, [], true, [3 12 42]);

% Reproducible random selection of illustrative traces:
rng(123);
summary = fig2_48_realisations(outDir, [], true, 'random');
```

Use run indices present in the generated dataset. The plot-only branch writes `summary_replotted.mat`, updated summary/latency tables and the standard figure files while preserving the original summary and raw realisations. A newer `summary_replotted.mat` is used in preference to `summary.mat`. Preserve the full output directory so selected trajectories and any derived metrics can be recovered.

## Statistical analysis

```matlab
results = analyse_fig2_statistics(outDir);
disp(results.publication_table);
```

See [Statistical analysis](README_statistics.md) for the comparison families, adjusted intervals and correlation convention. Conditions use independent simulation substreams; repeated states within a condition share the same realisation. Simulation settings, seed, model versions and the full results directory should accompany an archived analysis.
