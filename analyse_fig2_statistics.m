function results = analyse_fig2_statistics(source, varargin)
% ANALYSE_FIG2_STATISTICS  ANOVA and Bonferroni tables from saved Figure 2 data.
%
% Run after fig2_48_realisations; no simulation is repeated:
%   results = analyse_fig2_statistics('path/to/fig2_results');
%   results = analyse_fig2_statistics('path/to/summary_replotted.mat');
%   results = analyse_fig2_statistics(summary);
%
% The default spectral measure matches the figure's f * wavelet power:
%   results = analyse_fig2_statistics('path/to/fig2_results', ...
%       'PowerMeasure', 'power'); % optional UNWEIGHTED wavelet power instead
%
% Other options: 'Alpha', 0.05; 'OutputDir', 'path/to/statistics';
%                'SaveFiles', true (false returns tables without exporting).
%
% Conditions: indices [1 2 4] = Healthy, PD, PD severe. Mild is excluded.
% Unit of analysis: one scalar per independent realisation (normally n=48).
% Scalar comparisons: all three condition pairs, separately in each state:
%   velocity: rest, single programme, two programmes, hyperdirect;
%   flexor-extensor r: rest, single programme, two programmes;
%   target-position r: single programme and two programmes;
%   latency: single programme only, using the existing xcorr lag measure.
% Power comparisons: Healthy vs PD and Healthy vs severe in all four states,
%   plus REST vs HYPERDIRECT paired by realisation within each condition.
% At rest only, undefined flexor-extensor r in Healthy/PD is defined as zero
% when the stored covariance is exactly zero. In the simulation's metric
% helper, NaN r with zero covariance identifies a constant muscle signal;
% invalid/missing signals instead have NaN covariance and are NOT replaced.
% This is an explicit zero convention, not a calculated Pearson correlation.
% Original values and every substitution are exported in a separate audit.
%
% Bonferroni correction is WITHIN EACH METRIC ACROSS ITS REQUESTED STATES:
%   velocity 12, flexor-extensor r 9, target-position r 6, latency 3;
%   power 11 (eight vs Healthy plus three paired state comparisons).
% Planned family sizes are retained even if a comparison cannot be estimated.
% All requested comparisons are reported, regardless of the omnibus P value.
%
% For A vs B: difference = mean(B)-mean(A); percent = 100*difference/mean(A).
% Paired differences and means use complete pairs only. A zero reference mean
% gives an undefined percent change. Negative reference means retain their sign.
% CIs use Student-t critical values with the SAME Bonferroni family as P values.
% They are CIs for differences, not the plot's mean +/- 1.96*SE error bars.
% Pairwise t uses the pooled three-condition ANOVA residual variance.
% The publication-table F is t^2 for that contrast, not the omnibus ANOVA F.
% The two-state repeated-measures ANOVA is computed exactly as F=t_paired^2.
%
% Band measure per run = integral_3^7 S(f) df / (7-3), with endpoint
% interpolation and trapezoidal integration on the actual frequency grid.
% S(f) is f*mean_t(abs(CWT).^2) by default, or mean_t(abs(CWT).^2) for 'power'.
% This is a band-average wavelet measure, NOT calibrated PSD bandpower.
% All-state-sample spectra are used, matching the current figure (no COI mask).
% Correlations are analysed on the original Pearson-r scale, without Fisher z.
% Conditions are independent; states within a condition share realisations.
% Classical between-condition ANOVA assumes normal errors and equal variances.
%
% Output: Excel workbook, publication/detailed CSVs, methods TXT and results MAT
% under <source folder>/anova_bonferroni. Original summaries are never modified.
% Missing single-programme latency is recovered from saved realisation files.
% MATLAB R2022b-compatible; Statistics and Machine Learning Toolbox required.

if nargin < 1 || isempty(source)
    error('Fig2Stats:MissingSource', 'Provide a results folder, summary MAT file, or summary structure.');
end
p = inputParser;
addParameter(p, 'PowerMeasure', 'frequency_weighted_power', @(x) ischar(x) || isstring(x));
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x>0 && x<1);
addParameter(p, 'OutputDir', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'SaveFiles', true, @(x) islogical(x) && isscalar(x));
parse(p, varargin{:});
opt = p.Results;
opt.PowerMeasure = validatestring(opt.PowerMeasure, {'frequency_weighted_power','power'});
for dependency = {'anova1','tcdf','tinv'}
    if isempty(which(dependency{1}))
        error('Fig2Stats:MissingToolbox', '%s requires the Statistics and Machine Learning Toolbox.', dependency{1});
    end
end
[S, sourceDir, sourceName] = load_summary(source);
[conditions, stateIndex, R] = validate_summary(S);
conditionIndex = [1 2 4];
states = {'Rest','Single programme','Two programmes','Hyperdirect'};

metricData.velocity = double(S.metrics.velocity(conditionIndex,stateIndex,:));
metricData.rigidity_correlation = double(S.metrics.rigidity_correlation(conditionIndex,stateIndex,:));
[metricData.rigidity_correlation, restCorrelationConvention] = ...
    rest_zero_convention(S,metricData.rigidity_correlation,conditionIndex,stateIndex(1),R);
metricData.tracking_correlation = double(S.metrics.tracking_correlation(conditionIndex,stateIndex,:));
latency = single_latency(S, sourceDir, R, stateIndex(2));
metricData.latency_single = nan(3,4,R);
metricData.latency_single(:,2,:) = reshape(latency(conditionIndex,:),3,1,R);
[rawBand, weightedBand, bandGrid] = band_average(S, conditionIndex, stateIndex, R);
if strcmp(opt.PowerMeasure, 'frequency_weighted_power')
    metricData.band_power_3_7 = weightedBand;
    bandLabel = 'Mean frequency-weighted wavelet power 3-7 Hz';
    bandUnit = 'position^2 Hz';
    suffix = 'weighted';
else
    metricData.band_power_3_7 = rawBand;
    bandLabel = 'Mean wavelet power 3-7 Hz';
    bandUnit = 'position^2';
    suffix = 'unweighted';
end

spec = struct('key', {'velocity','rigidity_correlation','tracking_correlation', ...
    'latency_single','band_power_3_7'}, ...
    'label', {'Velocity','Flexor-extensor correlation','Target-position correlation', ...
    'Single-programme latency',bandLabel}, ...
    'unit', {'units/s','r','r','s',bandUnit}, ...
    'states', {1:4,1:3,2:3,2,1:4});
pairRows = repmat(pair_template(),0,1);
anovaRows = repmat(anova_template(),0,1);
descriptiveRows = repmat(descriptive_template(),0,1);

for k = 1:numel(spec)
    key = spec(k).key;
    for s = spec(k).states
        % Y rows are independent realisations; columns are the three conditions.
        Y = reshape(metricData.(key)(:,s,:),3,R).';
        [fit, description] = independent_anova(Y, spec(k), states{s}, conditions);
        anovaRows(end+1,1) = fit; %#ok<AGROW>
        descriptiveRows = [descriptiveRows; description]; %#ok<AGROW>
        if strcmp(key,'band_power_3_7')
            pairs = [1 2; 1 3];
        else
            pairs = [1 2; 1 3; 2 3];
        end
        for j = 1:size(pairs,1)
            a = pairs(j,1); b = pairs(j,2);
            row = pair_template();
            row.MetricKey = string(key); row.Metric = string(spec(k).label);
            row.Unit = string(spec(k).unit); row.State = string(states{s});
            row.Test = "Independent ANOVA contrast";
            row.A = string(conditions{a}); row.B = string(conditions{b});
            row.Comparison = row.A + " vs " + row.B;
            row.N_A = description(a).N; row.N_B = description(b).N;
            row.N_expected = R; row.Mean_A = description(a).Mean;
            row.Mean_B = description(b).Mean;
            row.Mean_difference = row.Mean_B-row.Mean_A;
            row.Change_percent = percent_change(row.Mean_A,row.Mean_difference);
            row.df = fit.df_den; row.Status = fit.Status;
            if fit.Status == "OK"
                row.SE_difference = sqrt(fit.ErrorMS*(1/row.N_A+1/row.N_B));
                row.t = row.Mean_difference/row.SE_difference;
                row.F = row.t^2;
                row.P_raw = min(1,2*tcdf(-abs(row.t),row.df));
            end
            pairRows(end+1,1) = row; %#ok<AGROW>
        end
    end
end

% Same run, same condition: the two states MUST be treated as paired.
for c = 1:3
    a = reshape(metricData.band_power_3_7(c,1,:),R,1);
    b = reshape(metricData.band_power_3_7(c,4,:),R,1);
    [row, fit] = paired_state_anova(a,b,spec(end),conditions{c});
    pairRows(end+1,1) = row; %#ok<AGROW>
    anovaRows(end+1,1) = fit; %#ok<AGROW>
end
pairRows = apply_bonferroni(pairRows,opt.Alpha);

results.created = datestr(now,30);
results.analysis_revision = '2026-10-07_rest_hyperdirect_zero_rest_correlation';
results.source = sourceName;
results.options = opt;
results.condition_indices = conditionIndex;
results.condition_labels = conditions;
results.state_names = S.state_names(stateIndex);
results.state_labels = states;
results.expected_realisations = R;
results.metric_values = metricData;
results.rest_correlation_zero_convention = restCorrelationConvention;
results.pairwise = struct2table(pairRows);
results.anova = struct2table(anovaRows);
results.descriptives = struct2table(descriptiveRows);
results.publication_table = publication_table(pairRows);
results.band_power.raw = rawBand;
results.band_power.frequency_weighted = weightedBand;
results.band_power.integration_grid_hz = bandGrid;
results.band_power.layout = 'condition x state x realisation';
results.band_power.per_realisation = band_table(rawBand,weightedBand,conditions,states,R);
results.methods = method_notes(opt,R);
results.files = struct();

if opt.SaveFiles
    outputDir = char(opt.OutputDir);
    if isempty(outputDir), outputDir = fullfile(sourceDir,'anova_bonferroni'); end
    results = export_results(results,outputDir,suffix);
end
fprintf('Completed %d contrasts and %d ANOVA analyses using up to %d runs per condition.\n', ...
    height(results.pairwise),height(results.anova),R);
fprintf('Bonferroni family sizes: velocity 12; flexor-extensor r 9; target-position r 6; latency 3; power 11.\n');
bad = sum(results.pairwise.Status ~= "OK");
if bad>0
    warning('Fig2Stats:UnestimableTests', ...
        '%d comparisons could not be estimated; see Status and sample counts in the detailed table.',bad);
end
end

function [S, folder, name] = load_summary(source)
if isstruct(source)
    S = source; name = 'summary structure supplied directly';
    if isfield(S,'output_dir') && isfolder(S.output_dir)
        folder = char(S.output_dir);
    else
        folder = pwd;
    end
    return;
end
source = char(source);
if isfolder(source)
    folder = source;
    original = dir(fullfile(folder,'summary.mat'));
    updated = dir(fullfile(folder,'summary_replotted.mat'));
    if ~isempty(updated) && (isempty(original) || updated.datenum>=original.datenum)
        name = fullfile(folder,'summary_replotted.mat');
    else
        name = fullfile(folder,'summary.mat');
    end
else
    name = source; folder = fileparts(source);
    if isempty(folder), folder = pwd; end
end
if ~isfile(name), error('Fig2Stats:MissingSummary','Cannot find %s.',name); end
loaded = load(name,'summary');
if ~isfield(loaded,'summary')
    error('Fig2Stats:MissingVariable','%s must contain a variable named summary.',name);
end
S = loaded.summary;
end

function [labels, stateIndex, R] = validate_summary(S)
needed = {'condition_labels','state_names','metrics','joint_spectra'};
for j = 1:numel(needed)
    assert(isfield(S,needed{j}),'Fig2Stats:Schema','Missing summary.%s.',needed{j});
end
assert(numel(S.condition_labels)==4 && numel(S.state_names)==4, ...
    'Fig2Stats:Schema','Use the four-condition, four-state Figure 2 summary.');
labels = cellstr(string(S.condition_labels([1 2 4])));
wanted = {'rest','movement_single','movement_two','hyperdirect'};
[found,stateIndex] = ismember(wanted,cellstr(string(S.state_names)));
assert(all(found),'Fig2Stats:Schema','The four required state names were not found.');
assert(isfield(S.metrics,'velocity'),'Fig2Stats:Schema','Missing per-realisation velocity.');
R = size(S.metrics.velocity,3);
assert(R>=2,'Fig2Stats:Schema','At least two realisations are needed.');
if isfield(S,'completed_realisations')
    assert(S.completed_realisations==R,'Fig2Stats:Schema','Realisation count does not match the metric arrays.');
end
for field = {'velocity','rigidity_correlation','rigidity_covariance','tracking_correlation'}
    assert(isfield(S.metrics,field{1}),'Fig2Stats:Schema','Missing summary.metrics.%s.',field{1});
    A = S.metrics.(field{1});
    assert(isnumeric(A) && isreal(A) && isequal(size(A),[4 4 R]), ...
        'Fig2Stats:Schema','%s must be condition x state x realisation (4 x 4 x N).',field{1});
end
assert(isfield(S.joint_spectra,'frequency_hz') && isfield(S.joint_spectra,'power'), ...
    'Fig2Stats:Schema','Per-realisation wavelet power and its frequency grid are required.');
if isfield(S,'config') && isfield(S.config,'random_substreams')
    streams = S.config.random_substreams(:,[1 2 4]);
    assert(numel(unique(streams(:)))==numel(streams), ...
        'Fig2Stats:Pairing','Condition substreams are not independent; this dataset needs a different design.');
end
end

function [values,convention] = rest_zero_convention(S,values,conditions,restState,R)
% The generator returns covariance=0 and r=NaN for constant muscle signals.
% For nonfinite input signals it returns NaN for BOTH quantities instead.
original = reshape(values(:,1,:),3,R);
covariance = reshape(double(S.metrics.rigidity_covariance(conditions,restState,:)),3,R);
velocity = reshape(double(S.metrics.velocity(conditions,restState,:)),3,R);
mask = false(3,R);
for c = 1:3
    if ismember(conditions(c),[1 2])
        mask(c,:) = isnan(original(c,:)) & isfinite(covariance(c,:)) & covariance(c,:)==0;
    end
end
analysed = original;
analysed(mask) = 0;
values(:,1,:) = reshape(analysed,3,1,R);
[c,r] = ndgrid(1:3,1:R);
labels = S.condition_labels(conditions);
Condition = reshape(string(labels(c(:))),[],1);
Realisation = r(:);
Original_rest_correlation = original(:);
Rest_covariance = covariance(:);
Rest_velocity = velocity(:);
Analysed_rest_correlation = analysed(:);
Zero_by_convention = mask(:);
convention.applied_mask = mask;
convention.n_replaced_by_condition = sum(mask,2);
convention.condition_indices = conditions;
convention.definition = ['Healthy/PD at rest: NaN flexor-extensor correlation with finite zero ' ...
    'covariance is set to zero by convention. This identifies a constant muscle signal ' ...
    'in the original metric helper. Other missing correlations remain missing.'];
convention.audit = table(Condition,Realisation,Original_rest_correlation,Rest_covariance, ...
    Rest_velocity,Analysed_rest_correlation,Zero_by_convention);
fprintf('Rest correlation zero convention: Healthy %d runs; PD %d runs.\n', ...
    convention.n_replaced_by_condition(1),convention.n_replaced_by_condition(2));
end

function values = single_latency(S,folder,R,singleState)
values = [];
if isfield(S,'latency_single') && isfield(S.latency_single,'values_s')
    L = S.latency_single;
    compatible = isequal(size(L.values_s),[4 R]) && ...
        isfield(L,'method_version') && L.method_version==1;
    if isfield(S,'state_times_s')
        compatible = compatible && isfield(L,'window_s') && ...
            isequal(L.window_s,S.state_times_s(singleState,:));
    end
    if compatible, values = double(L.values_s); return; end
end
assert(~isempty(which('xcorr')),'Fig2Stats:MissingLatencyDependency', ...
    'xcorr is needed to recover single-programme latency from saved traces.');
values = nan(4,R);
fprintf('Recovering single-programme latency from saved traces (no simulation).\n');
for r = 1:R
    file = fullfile(folder,sprintf('realisation_%03d.mat',r));
    if ~isfile(file) && isfield(S,'files') && numel(S.files)>=r && isfile(S.files{r})
        file = S.files{r};
    end
    assert(isfile(file),'Fig2Stats:MissingLatency', ...
        'Single-programme latency is absent and the saved run is missing: %s.',file);
    saved = load(file,'runData'); D = saved.runData;
    assert(D.realisation==r && numel(D.conditions)==4, ...
        'Fig2Stats:RunMismatch','Saved realisation does not match the summary.');
    assert(isequal(D.config.condition_labels,S.condition_labels) && ...
        D.config.base_seed==S.config.base_seed && D.config.fs==S.config.fs && ...
        isequal(D.inputs.phase_times_s(singleState,:),S.state_times_s(singleState,:)), ...
        'Fig2Stats:RunMismatch','Saved run settings or state windows do not match the summary.');
    idx = D.inputs.phase_indices.movement_single;
    target = D.inputs.active_target(idx); target = target(:);
    for c = [1 2 4]
        actual = D.conditions{c}.muscleposition(idx); actual = actual(:);
        if numel(actual)<2 || numel(actual)~=numel(target) || ...
                any(~isfinite(actual)) || any(~isfinite(target)) || ...
                all(actual==actual(1)) || all(target==target(1))
            continue;
        end
        [crossCov,lags] = xcorr(actual-mean(actual),target-mean(target));
        [~,peak] = max(crossCov);
        values(c,r) = lags(peak)/D.config.fs;
    end
end
end

function [raw,weighted,grid] = band_average(S,conditions,states,R)
f = double(S.joint_spectra.frequency_hz(:));
assert(numel(f)>=2 && all(isfinite(f) & f>0), ...
    'Fig2Stats:FrequencyGrid','Frequency grid must contain finite positive frequencies.');
[f,order] = sort(f);
assert(all(diff(f)>0) && f(1)<=3 && f(end)>=7, ...
    'Fig2Stats:FrequencyGrid','A unique frequency grid covering the full 3-7 Hz band is required.');
P = S.joint_spectra.power;
assert(isnumeric(P) && isreal(P) && isequal(size(P),[4 4 R numel(f)]), ...
    'Fig2Stats:SpectrumShape','Power must be 4 x 4 x N x frequency, retaining individual runs.');
P = double(P(conditions,states,:,order));
curves = reshape(P,[],numel(f));
assert(~any(curves(isfinite(curves))<0),'Fig2Stats:NegativePower','Wavelet power cannot be negative.');
lo = find(f<=3,1,'last'); hi = find(f>=7,1,'first');
f = f(lo:hi); curves = curves(:,lo:hi);
grid = unique([3; f(f>3 & f<7); 7]);
valid = all(isfinite(curves),2);
rawFlat = nan(size(curves,1),1); weightedFlat = rawFlat;
if any(valid)
    % Integrate on the real frequency grid; CWT bins are not equally spaced.
    yi = interp1(f,curves(valid,:).',grid,'linear');
    weightedCurves = bsxfun(@times,curves(valid,:),f.');
    ywi = interp1(f,weightedCurves.',grid,'linear');
    rawFlat(valid) = (trapz(grid,yi,1)/4).';
    weightedFlat(valid) = (trapz(grid,ywi,1)/4).';
end
raw = reshape(rawFlat,3,4,R);
weighted = reshape(weightedFlat,3,4,R);
end

function [fit, descriptions] = independent_anova(Y,spec,state,conditions)
fit = anova_template();
fit.MetricKey = string(spec.key); fit.Metric = string(spec.label);
fit.State = string(state); fit.Test = "One-way independent ANOVA";
fit.Condition = "All three plotted conditions";
descriptions = repmat(descriptive_template(),3,1);
counts = zeros(1,3);
for c = 1:3
    v = Y(isfinite(Y(:,c)),c); counts(c) = numel(v);
    d = descriptive_template();
    d.Metric = string(spec.label); d.State = string(state);
    d.Condition = string(conditions{c}); d.Unit = string(spec.unit);
    d.N = numel(v); d.Missing = size(Y,1)-d.N;
    if ~isempty(v), d.Mean = mean(v); end
    if numel(v)>=2
        d.STD = std(v,0); d.SE = d.STD/sqrt(d.N);
        d.Mean_CI95_lower = d.Mean-1.96*d.SE;
        d.Mean_CI95_upper = d.Mean+1.96*d.SE;
    end
    descriptions(c) = d;
end
fit.N = sum(counts); fit.Missing = numel(Y)-fit.N;
fit.df_num = 2; fit.df_den = fit.N-3;
if any(counts<2)
    fit.Status = "Fewer than two valid runs in at least one condition"; return;
end
y = Y(:); group = repelem((1:3).',size(Y,1));
valid = isfinite(y);
[p,T] = anova1(y(valid),group(valid),'off');
fit.F = T{2,5}; fit.P_raw = p; fit.ErrorMS = T{3,4};
fit.df_num = T{2,3}; fit.df_den = T{3,3};
if ~isfinite(fit.ErrorMS) || fit.ErrorMS<=0
    fit.F = NaN; fit.P_raw = NaN;
    fit.Status = "Zero or invalid residual variance; inference unavailable";
else
    fit.Status = "OK";
end
end

function [row,fit] = paired_state_anova(a,b,spec,condition)
row = pair_template();
row.MetricKey = string(spec.key); row.Metric = string(spec.label);
row.Unit = string(spec.unit); row.State = "Rest vs Hyperdirect";
row.Test = "Paired two-level repeated-measures ANOVA";
row.A = string(condition)+": Rest";
row.B = string(condition)+": Hyperdirect";
row.Comparison = string(condition)+": Rest vs Hyperdirect";
row.N_expected = numel(a);
valid = isfinite(a) & isfinite(b); a = a(valid); b = b(valid);
n = numel(a); row.N_A = n; row.N_B = n; row.df = n-1;
fit = anova_template();
fit.MetricKey = row.MetricKey; fit.Metric = row.Metric;
fit.State = row.State; fit.Test = row.Test; fit.Condition = string(condition);
fit.N = n; fit.Missing = row.N_expected-n; fit.df_num = 1; fit.df_den = n-1;
if n>0
    row.Mean_A = mean(a); row.Mean_B = mean(b);
    row.Mean_difference = mean(b-a);
    row.Change_percent = percent_change(row.Mean_A,row.Mean_difference);
end
if n<2
    row.Status = "Fewer than two complete pairs";
elseif ~isfinite(var(b-a,0)) || var(b-a,0)<=0
    row.Status = "Zero or invalid variance of paired differences; inference unavailable";
else
    % For two repeated levels: SS_state=n*mean(d)^2/2, MS_error=var(d)/2.
    % Hence F=n*mean(d)^2/var(d)=t_paired^2, with df=(1,n-1).
    row.SE_difference = std(b-a,0)/sqrt(n);
    row.t = row.Mean_difference/row.SE_difference;
    row.F = row.t^2; row.P_raw = min(1,2*tcdf(-abs(row.t),n-1));
    row.Status = "OK"; fit.ErrorMS = var(b-a,0)/2;
    fit.F = row.F; fit.P_raw = row.P_raw;
end
fit.Status = row.Status;
end

function rows = apply_bonferroni(rows,alpha)
keys = [rows.MetricKey];
for key = unique(keys,'stable')
    indices = find(keys==key);
    m = numel(indices);
    for j = indices
        rows(j).Bonferroni_family = key;
        rows(j).Family_comparisons = m;
        rows(j).Alpha_family = alpha;
        if rows(j).Status ~= "OK", continue; end
        rows(j).P_bonferroni = min(1,m*rows(j).P_raw);
        critical = tinv(1-alpha/(2*m),rows(j).df);
        margin = critical*rows(j).SE_difference;
        rows(j).CI_lower = rows(j).Mean_difference-margin;
        rows(j).CI_upper = rows(j).Mean_difference+margin;
        rows(j).Significant = rows(j).P_bonferroni<alpha;
    end
end
end

function value = percent_change(reference,difference)
if ~isfinite(reference) || reference==0 || ~isfinite(difference)
    value = NaN;
else
    value = 100*difference/reference;
end
end

function T = publication_table(rows)
cells = cell(0,6); previous = "";
for j = 1:numel(rows)
    r = rows(j);
    heading = r.Metric+" ("+r.Unit+") | "+r.State;
    if heading~=previous
        cells(end+1,:) = {char(heading),'','','','',''}; %#ok<AGROW>
        previous = heading;
    end
    effect = sprintf('%s (%s, %s)',number_text(r.Mean_difference), ...
        number_text(r.CI_lower),number_text(r.CI_upper));
    cells(end+1,:) = {char(r.Comparison),effect,number_text(r.Change_percent), ...
        number_text(r.t),number_text(r.F),p_text(r.P_bonferroni)}; %#ok<AGROW>
end
T = cell2table(cells,'VariableNames',{'Comparison / condition', ...
    'Mean difference or change (lower, upper)','Change (%)','t','F','P value'});
end

function str = number_text(value)
if isnan(value), str = 'NA'; else, str = sprintf('%.4g',value); end
end

function str = p_text(value)
if isnan(value)
    str = 'NA';
elseif value==0
    str = sprintf('<%.3g',realmin);
else
    str = sprintf('%.3g',value);
end
end

function T = band_table(raw,weighted,conditions,states,R)
[c,s,r] = ndgrid(1:3,1:4,1:R);
originalConditions = [1 2 4];
ConditionIndex = reshape(originalConditions(c(:)),[],1);
Condition = reshape(string(conditions(c(:))),[],1);
State = reshape(string(states(s(:))),[],1);
Realisation = r(:);
WaveletPower_3_7Hz = raw(:);
FrequencyWeightedPower_3_7Hz = weighted(:);
T = table(ConditionIndex,Condition,State,Realisation, ...
    WaveletPower_3_7Hz,FrequencyWeightedPower_3_7Hz);
end

function notes = method_notes(opt,R)
notes = [
    "One observation per realisation; expected runs per condition: "+string(R)+". Valid sample counts are reported."
    "Conditions: Healthy, PD and PD severe (indices 1,2,4); mild adaptations excluded."
    "Scalar comparisons: all three condition pairs for each state displayed in the main figure."
    "Between conditions: separate classical one-way ANOVA for each metric/state; post hoc t uses that ANOVA's pooled residual variance."
    "Power: PD and severe vs Healthy in all four states; rest vs hyperdirect paired within each of the three conditions (change=hyperdirect-rest)."
    "Two repeated states: exact repeated-measures ANOVA F=t_paired^2, df=(1,n_complete_pairs-1)."
    "Bonferroni families span all requested comparisons within each metric, including its states."
    "Family sizes: velocity=12; flexor-extensor r=9; target-position r=6; latency=3; power=11."
    "Family alpha="+string(opt.Alpha)+". Adjusted P=min(1,m*P_raw). CI=difference +/- tinv(1-alpha/(2*m),df)*SE_difference."
    "Publication-table P values and difference CIs are Bonferroni-adjusted. Omnibus ANOVA P values in the ANOVA sheet are unadjusted."
    "Publication-table F is the pairwise contrast t^2, not the overall three-group ANOVA F. The overall statistics are exported separately."
    "Difference for A vs B is B-A; percent change=100*(B-A)/A. Paired estimates use complete pairs for both means and the difference."
    "Percent change is undefined at zero reference mean; near-zero or negative references, including correlations, need interpretation alongside the means."
    "Correlations retain the raw Pearson-r scale, with no Fisher transformation. All requested contrasts are reported regardless of omnibus significance."
    "Resting flexor-extensor correlation: Healthy/PD NaN values with finite zero covariance are assigned zero by convention for the comparisons."
    "The original metric helper gives zero covariance and NaN r for a constant muscle signal; missing/nonfinite signals give NaN covariance and are not converted."
    "These assigned zeros are not measured Pearson correlations. Each replacement is recorded in the Rest correlation sheet and CSV; original summaries are preserved."
    "Zero reference means produce NA percentage changes. If zero-convention groups are constant while another group varies, the equal-variance ANOVA assumption is not met; interpret its inferential results accordingly."
    "Nonfinite values are excluded separately by group; paired analyses exclude an entire pair when either value is nonfinite."
    "Missing or zero-variance tests are flagged; planned Bonferroni family sizes are not reduced."
    "Band measure: frequency-average over 3-7 Hz, trapz on the actual grid divided by 4 Hz; exact band endpoints are linearly interpolated."
    "Selected spectral measure: "+string(opt.PowerMeasure)+". Both unweighted and frequency-weighted band values are exported per realisation."
    "Wavelet power is mean_t(abs(CWT).^2), not PSD per Hz; frequency-weighted power multiplies this by f. These are wavelet measures, not calibrated physical bandpower."
    "Spectra include all state samples, matching the figure; no cone-of-influence exclusion or duration pooling is added."
    "Single-programme latency uses the saved mean-centred cross-correlation lag in seconds. No circular-phase analysis or latency method change is introduced."
    "Descriptive mean intervals match the figure: mean +/- 1.96*SE. They differ from the Bonferroni-adjusted intervals for differences."
    "Classical between-condition ANOVA assumes normal errors and equal variances. For the paired test, normality concerns within-run differences."
    "The original summary and simulation files are not modified."
    ];
end

function results = export_results(results,folder,suffix)
if ~isfolder(folder)
    [ok,message] = mkdir(folder);
    assert(ok,'Fig2Stats:OutputDirectory','Could not create output folder: %s',message);
end
prefix = fullfile(folder,['fig2_anova_bonferroni_' suffix]);
results.files.publication_csv = [prefix '_table.csv'];
results.files.pairwise_csv = [prefix '_pairwise.csv'];
results.files.anova_csv = [prefix '_anova.csv'];
results.files.descriptives_csv = [prefix '_descriptives.csv'];
results.files.band_power_csv = [prefix '_band_power_per_realisation.csv'];
results.files.rest_correlation_csv = [prefix '_rest_correlation_convention.csv'];
results.files.methods_txt = [prefix '_methods.txt'];
results.files.excel = [prefix '.xlsx'];
results.files.mat = [prefix '.mat'];
writetable(results.publication_table,results.files.publication_csv);
writetable(results.pairwise,results.files.pairwise_csv);
writetable(results.anova,results.files.anova_csv);
writetable(results.descriptives,results.files.descriptives_csv);
writetable(results.band_power.per_realisation,results.files.band_power_csv);
writetable(results.rest_correlation_zero_convention.audit,results.files.rest_correlation_csv);
fid = fopen(results.files.methods_txt,'wt');
assert(fid>=0,'Fig2Stats:OutputFile','Cannot write %s.',results.files.methods_txt);
closer = onCleanup(@() fclose(fid));
for j = 1:numel(results.methods)
    fprintf(fid,'%s\n',char(results.methods(j)));
end
clear closer;
results.excel_saved = false;
try
    if isfile(results.files.excel), delete(results.files.excel); end
    writetable(results.publication_table,results.files.excel,'Sheet','Publication','UseExcel',false);
    writetable(results.pairwise,results.files.excel,'Sheet','Pairwise','UseExcel',false);
    writetable(results.anova,results.files.excel,'Sheet','ANOVA','UseExcel',false);
    writetable(results.descriptives,results.files.excel,'Sheet','Descriptives','UseExcel',false);
    writetable(results.band_power.per_realisation,results.files.excel,'Sheet','Band power per run','UseExcel',false);
    writetable(results.rest_correlation_zero_convention.audit,results.files.excel,'Sheet','Rest correlation','UseExcel',false);
    writetable(table(results.methods,'VariableNames',{'Method'}), ...
        results.files.excel,'Sheet','Methods','UseExcel',false);
    results.excel_saved = true;
catch ME
    warning('Fig2Stats:ExcelExport','Excel export failed: %s. CSV and MAT outputs are retained.',ME.message);
end
save(results.files.mat,'results','-v7.3');
fprintf('Analysis files saved in %s\n',folder);
end

function row = pair_template()
row = struct('MetricKey',"",'Metric',"",'Unit',"",'State',"",'Test',"", ...
    'A',"",'B',"",'Comparison',"",'N_A',NaN,'N_B',NaN,'N_expected',NaN, ...
    'Mean_A',NaN,'Mean_B',NaN,'Mean_difference',NaN,'Change_percent',NaN, ...
    'SE_difference',NaN,'df',NaN,'t',NaN,'F',NaN,'P_raw',NaN, ...
    'Bonferroni_family',"",'Family_comparisons',NaN,'Alpha_family',NaN, ...
    'CI_lower',NaN,'CI_upper',NaN,'P_bonferroni',NaN,'Significant',false,'Status',"");
end

function row = anova_template()
row = struct('MetricKey',"",'Metric',"",'State',"",'Test',"",'Condition',"", ...
    'N',NaN,'Missing',NaN,'df_num',NaN,'df_den',NaN,'F',NaN,'P_raw',NaN, ...
    'ErrorMS',NaN,'Status',"");
end

function row = descriptive_template()
row = struct('Metric',"",'State',"",'Condition',"",'Unit',"", ...
    'N',NaN,'Missing',NaN,'Mean',NaN,'STD',NaN,'SE',NaN, ...
    'Mean_CI95_lower',NaN,'Mean_CI95_upper',NaN);
end
