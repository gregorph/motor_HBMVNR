# Statistical analysis

The analysis functions operate on generated simulation outputs. Independent realisations or runs are the statistical units. Before/after measurements are paired within a run; individual time samples and overlapping moving windows are not treated as independent replicates.

| Entry point | Input | Purpose |
| --- | --- | --- |
| `analyse_fig2_statistics` | Figure 2 output directory, summary MAT file or summary structure | Between-condition contrasts and paired state comparisons |
| `analyse_fig3b_dbs` | DBS position array and matching target | Before/after DBS metrics and inference |
| `statistics_dbs.m` | Generated switch-on dataset | Driver and caption-table output for DBS analysis |
| `extract_pharmacological_stats.m` | Generated pharmacology FIG file | Pharmacological comparisons of plotted distributions |

Statistical inference requires Statistics and Machine Learning Toolbox. The DBS function also supports metric calculation without inferential tests through `RunStatistics=false`.

## Figure 2

After generating the full dataset:

```matlab
outDir = fullfile(pwd, 'fig2_results_01');
results = analyse_fig2_statistics(outDir);
disp(results.publication_table);
```

The function also accepts a specific summary file or an existing structure:

```matlab
results = analyse_fig2_statistics(fullfile(outDir, 'summary.mat'));
results = analyse_fig2_statistics(summary, 'SaveFiles', false);
```

| Option | Default | Description |
| --- | --- | --- |
| `PowerMeasure` | `'frequency_weighted_power'` | Use the displayed frequency-weighted wavelet measure; `'power'` selects unweighted power |
| `Alpha` | 0.05 | Family-level significance threshold |
| `OutputDir` | `anova_bonferroni` under the source directory | Analysis destination |
| `SaveFiles` | `true` | Export analysis files |

The required schema has four conditions, four states and at least two realisations. The analysis includes conditions 1, 2 and 4: Healthy, PD and PD severe. The mild-adaptation condition is excluded. A newer `summary_replotted.mat` is preferred when the input is a directory; supplying an explicit MAT filename fixes the source selection.

### Comparisons

| Metric | Comparisons | Bonferroni family size |
| --- | --- | --- |
| Velocity | Three condition pairs in each of four states | 12 |
| Flexor/extensor correlation | Three pairs in rest, single-programme and two-programme states | 9 |
| Target/position correlation | Three pairs in the two movement states | 6 |
| Single-programme latency | Three condition pairs | 3 |
| 3–7 Hz spectral measure | Healthy versus each PD condition in four states, plus rest versus hyperdirect within each of the three conditions | 11 |

Between-condition comparisons use independent-condition ANOVA contrasts. Rest/hyperdirect spectral comparisons are paired within condition and realisation. For a comparison labelled A versus B, the reported difference is mean(B) − mean(A). Percent change uses mean(A) as the reference.

The publication-table F statistic is the contrast's `t^2`; omnibus ANOVA results are reported separately. Difference intervals use Student-t critical values with the same Bonferroni family as the P values. These intervals differ from the figures' mean ± 1.96 × SE intervals.

The spectral measure is the average over 3–7 Hz of frequency-weighted wavelet power, calculated by trapezoidal integration with interpolated band endpoints and division by the band width. Selecting `'power'` applies the same operation to unweighted squared wavelet magnitude. These are wavelet measures rather than calibrated PSD bandpower.

At rest, undefined flexor/extensor correlations in Healthy and PD are assigned zero only when stored covariance is exactly zero. This explicit constant-signal convention is recorded in an audit table. Original simulation summaries are preserved.

### Exported outputs

Analysis files use the prefix `fig2_anova_bonferroni_weighted`, or `fig2_anova_bonferroni_unweighted` for the alternative measure:

| Output | Contents |
| --- | --- |
| `_table.csv` | Formatted publication table |
| `_pairwise.csv` | Contrast estimates, statistics, adjusted P values and status |
| `_anova.csv` | ANOVA results |
| `_descriptives.csv` | Descriptive statistics |
| `_band_power_per_realisation.csv` | Per-realisation spectral measures |
| `_rest_correlation_convention.csv` | Constant-signal convention audit |
| `_methods.txt` | Analysis-method description |
| `.xlsx` | Corresponding worksheets |
| `.mat` | Complete `results` structure |

Use a different `OutputDir` to retain multiple analysis variants. Valid counts and test status accompany the estimates; numerical values are retained independently of publication formatting.

## DBS before/after analysis

`analyse_fig3b_dbs` uses the switch-on data produced by the [DBS workflow](README_fig4_dbsfig_VLnucleus.md). The position input has layout **condition × run × sample**, with four conditions, an even number of runs and at least 9,000 samples.

| Index | Analysis label |
| --- | --- |
| 1 | Control: unstimulated reference |
| 2 | STN |
| 3 | GPi |
| 4 | VP/VL |

Odd runs are rest and even runs are movement. With 48 runs, resting-tremor measures use 24 rest runs and tracking measures use 24 movement runs per condition. The VP/VL label identifies the thalamic condition; the simulation also applies its thalamic DBS mechanism to CM circuitry.

```matlab
D = load('fig4_switchon.mat', 'fig3b_muscleposition', 'desiredposition');
nTime = size(D.fig3b_muscleposition,3);
R = analyse_fig3b_dbs(D.fig3b_muscleposition, ...
    D.desiredposition(1,1:nTime), ...
    'Fs', 250, 'TargetFrequency', 0.1, 'MaxLag', 500, ...
    'NumPermutations', 19999, 'NumBootstraps', 20000, 'RandomSeed', 1);

disp(R.tables.latencyPhaseSeconds);
disp(R.tables.covariance);
disp(R.tables.tremor);
```

The original target and position must share time alignment. A common target vector or appropriately indexed per-run/per-condition targets are supported, with time as the last dimension. One extra final target sample can be trimmed by the function.

### Windows and options

Before/after windows are fixed inclusive sample indices **1000:4000** and **6000:9000**. At 250 Hz they correspond approximately to 4–16 s and 24–36 s. These windows are intended for the 40-second switch-on protocol. Changing `Fs` does not change the sample-index windows.

| Option | Default | Meaning |
| --- | --- | --- |
| `Fs` | 250 | Sampling rate, Hz |
| `TargetFrequency` | 0.1 | Sinusoidal target frequency, Hz |
| `MaxLag` | 500 | Forward covariance-lag limit in samples; 2 s at 250 Hz |
| `RunStatistics` | `true` | Calculate inferential statistics |
| `NumPermutations` | 19999 | Phase-label permutation budget |
| `NumBootstraps` | 20000 | Circular-change bootstrap resamples |
| `RandomSeed` | 1 | Local resampling seed |

The resampling seed does not control simulation random numbers. `MaxLag` must be an integer from 1 to 2999.

### Linear metrics

`R.summary` contains covariance, maximum covariance, covariance-based latency, resting tremor, fitted tracking amplitude and fitted tracking R-squared. Change is **after − before**. Before/after tests are paired within each condition, with `F=t^2` for the two-level repeated-measures comparison. P values and change intervals are Bonferroni-adjusted across the four conditions separately within each metric. The function does not perform between-target tests.

Resting tremor is the mean of a 100-sample moving variance of position within each analysis window. Covariance-based latency searches forward lags on equal-length supports. Its maximum covariance should be compared with `lagZeroCovariance`, which uses the same support, when quantifying the gain from lag alignment.

### Circular phase latency

`R.phaseSummary` and `R.tables.latencyPhaseSeconds` contain circular phase results. Sinusoidal fits at the target frequency estimate phase delay in each window. Positive delay indicates that the achieved trajectory lags the target. At 0.1 Hz, the principal range is `[-5,5)` seconds, modulo the 10-second period; `MaxLag` does not bound this estimate.

| Field | Interpretation |
| --- | --- |
| `CircularChangeSeconds` | Circular mean of paired phase-delay changes |
| `RBefore`, `RAfter`, `RChange` | Resultant lengths describing phase concentration |
| `D` | Distance between mean before/after phase vectors |
| `pBonferroni` | Adjusted paired-permutation P value for the phase-vector statistic |
| `CI95BonferroniLower`, `CI95BonferroniUpper` | Approximate paired-bootstrap interval for circular change |

The phase-vector test is sensitive to changes in direction or concentration; its P value is not solely a test of the displayed mean change. Circular phase results have no t statistic, F statistic or percent change. Retain their separate table rather than combining them with the linear latency statistics. Fitted amplitude, fit quality, pair counts and interval status support interpretation of the phase estimates.

### Save the analysis

```matlab
save('fig4_statistics.mat', 'R', '-v7.3');
writetable(R.summary, 'fig4_linear_statistics.csv');
writetable(R.phaseSummary, 'fig4_phase_statistics.csv');
```

The driver `statistics_dbs.m` reconstructs the standard 0.3-amplitude, 0.1-Hz target and displays selected tables. To use newly generated data, set its initial load to `load(fullfile(pwd,'fig4_switchon.mat'));`. For other target configurations, call the function directly with the actual saved target.

## Pharmacological comparisons

First generate and export `fig3_pharmacology.fig` using the [pharmacology workflow](README_fig3.md). In a working copy of `extract_pharmacological_stats.m`, set:

```matlab
filename1 = fullfile(pwd, 'fig3_pharmacology.fig');
```

Run the extraction script. It expects the original panel/object ordering of 15 violin objects, grouped into tremor, correlation and latency across Healthy, PD, Levodopa, Anticholinergic and Amantadine. The inputs are the plotted distributions, including the figure's correlation clipping and absolute-latency transformation.

The script applies one-way ANOVA and Bonferroni pairwise comparisons. Its formatted comparisons report second-listed minus first-listed group, and pairwise `F=t^2`; the caption block reports omnibus ANOVA statistics separately. Percentage changes use Healthy as reference when present in the pair and PD otherwise.

Preserve the extracted arrays and group labels:

```matlab
save('fig3_statistics_inputs.mat', ...
    'tremor', 'correlation', 'latency', 'groupNames', '-v7.3');
```

Archive analysis outputs with the generated simulation data, plotting configuration and source files so the relationship between numerical measures and displayed panels remains explicit.
