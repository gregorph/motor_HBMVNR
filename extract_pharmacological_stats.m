filename1 = "C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\in_preparation\PD_ssc4\figures\fig_pharmacological.fig";
fig = openfig(filename1, 'invisible');
v = findobj(fig, 'Type', 'violinplot');

for i = 1:15

    % Extract data from each violin plot and perform analysis
    data(i,:) = v(i).YData;
    % figure
    % violinplot(data(i,:))
    % (Add your analysis code here, e.g., calculate mean or standard deviation)

end

tremor(1:5,:) = flip(data(11:15,:));
correlation(1:5,:) = flip(data(6:10,:));
latency(1:5,:) = flip(data(1:5,:));

tremor = tremor.';
correlation = correlation.';
latency = latency.';


% Create sample data and run ANOVA
[p, tbl, stats] = anova1(correlation, [], 'off');

% Perform post hoc multiple comparisons (Tukey-Kramer by default)
if p < 0.05
    results = multcompare(stats, 'CType','bonferroni');
end


%%
%% =========================================================
%  One-way ANOVA + pairwise multiple comparisons
% ==========================================================

% Your array should be:
% data(:,1) = Healthy
% data(:,2) = PD
% data(:,3) = Levodopa
% data(:,4) = Anticholinergic
% data(:,5) = Amantadine

%% Group names
groupNames = {'Healthy', ...
              'PD', ...
              'Levodopa', ...
              'Anticholinergic', ...
              'Amantadine'};

data = tremor;

%% Convert array into a table
T = array2table(data, 'VariableNames', groupNames);

% Display the table
disp('Data table:')
disp(T)

%% =========================================================
% One-way ANOVA
% ==========================================================

[p, anovaTable, stats] = anova1(data, groupNames);

fprintf('\n================ ANOVA RESULT ================\n');
fprintf('ANOVA p-value = %.6f\n', p);

if p < 0.05
    fprintf('There is a statistically significant difference between groups.\n');
else
    fprintf('No statistically significant overall difference was detected.\n');
end

%% =========================================================
% Post-hoc pairwise comparisons
% ==========================================================

% Tukey-Kramer multiple comparison test
results = multcompare(stats, 'CType', 'bonferroni');

% multcompare returns:
% Column 1 = first group
% Column 2 = second group
% Column 3 = lower confidence interval
% Column 4 = estimated mean difference
% Column 5 = upper confidence interval
% Column 6 = adjusted p-value

%% Convert numeric group numbers to group names

nComparisons = size(results,1);

Group1 = strings(nComparisons,1);
Group2 = strings(nComparisons,1);

for i = 1:nComparisons
    Group1(i) = groupNames{results(i,1)};
    Group2(i) = groupNames{results(i,2)};
end

%% Create readable results table

pairwiseResults = table( ...
    Group1, ...
    Group2, ...
    results(:,3), ...
    results(:,4), ...
    results(:,5), ...
    results(:,6), ...
    'VariableNames', { ...
    'Group1', ...
    'Group2', ...
    'LowerCI', ...
    'MeanDifference', ...
    'UpperCI', ...
    'pValue'});

%% Add significance column

pairwiseResults.Significant = pairwiseResults.pValue < 0.05;

%% Display all pairwise comparisons

fprintf('\n================ ALL PAIRWISE COMPARISONS ================\n');
disp(pairwiseResults)

%% Display only significant comparisons

significantResults = pairwiseResults(pairwiseResults.Significant, :);

fprintf('\n================ SIGNIFICANT COMPARISONS ================\n');

if isempty(significantResults)
    fprintf('No pairwise comparisons were significant at p < 0.05.\n');
else
    disp(significantResults)
end

%% =========================================================
% Optional: print significant comparisons clearly
% ==========================================================

if ~isempty(significantResults)

    fprintf('\nSignificant differences:\n');

    for i = 1:height(significantResults)

        fprintf('%s vs %s: p = %.6f\n', ...
            significantResults.Group1(i), ...
            significantResults.Group2(i), ...
            significantResults.pValue(i));

    end
end








%%
groupNames = {'Healthy', 'PD', 'Levodopa', ...
    'Anticholinergic', 'Amantadine'};

measureNames = {'Correlation', 'Tremor', 'Latency'};
measureData = {correlation, tremor, latency};

for k = 1:numel(measureNames)

    data = measureData{k};

    [~, ~, stats] = anova1(data, groupNames, 'off');

    results = multcompare(stats, ...
        'CType', 'bonferroni', ...
        'Display', 'off');

    healthyIndex = find(strcmp(stats.gnames, 'Healthy'), 1);
    pdIndex = find(strcmp(stats.gnames, 'PD'), 1);

    againstHealthy = ...
        results(:,1) == healthyIndex | ...
        results(:,2) == healthyIndex;

    againstPD = ...
        results(:,1) == pdIndex | ...
        results(:,2) == pdIndex;

    % Include all Healthy/PD comparisons and other significant pairs.
    keep = againstHealthy | againstPD | results(:,6) < 0.05;
    comparisonsToPrint = results(keep, :);

    fprintf('\n%s\n', upper(measureNames{k}));
    fprintf('%s\n', repmat('-', 1, 129));

    fprintf('%-32s %-43s %12s %12s %12s %13s\n', ...
        'Comparison / condition', ...
        'Mean difference or change (lower, upper)', ...
        'Change (%)', 't', 'F', 'P value');

    fprintf('%s\n', repmat('-', 1, 129));

    for i = 1:size(comparisonsToPrint,1)

        g1 = comparisonsToPrint(i,1);
        g2 = comparisonsToPrint(i,2);

        % Difference = second-listed group minus first-listed group.
        meanDifference = -comparisonsToPrint(i,4);

        % Reverse confidence limits to match the difference.
        lowerCI = -comparisonsToPrint(i,5);
        upperCI = -comparisonsToPrint(i,3);

        pValue = comparisonsToPrint(i,6);

        SE = stats.s * sqrt(1/stats.n(g1) + 1/stats.n(g2));

        tValue = meanDifference / SE;

        % Pairwise F statistic, not the overall ANOVA F.
        FValue = tValue^2;

        % Percentage reference: Healthy if in the pair; otherwise PD.
        if g1 == healthyIndex || g2 == healthyIndex
            referenceMean = stats.means(healthyIndex);
        else
            referenceMean = stats.means(pdIndex);
        end

        if referenceMean ~= 0
            changePercent = ...
                100 * meanDifference / referenceMean;
            changeText = sprintf('%.2f%%', changePercent);
        else
            changeText = 'Undefined';
        end

        % Keep comparison names in their original order.
        comparisonText = sprintf('%s vs %s', ...
            stats.gnames{g1}, stats.gnames{g2});

        differenceText = sprintf('%.5g (%.5g, %.5g)', ...
            meanDifference, lowerCI, upperCI);

        fprintf('%-32s %-43s %12s %12.5g %12.5g %13.5g\n', ...
            comparisonText, ...
            differenceText, ...
            changeText, ...
            tValue, ...
            FValue, ...
            pValue);
    end

    fprintf('%s\n', repmat('-', 1, 129));
end

%%
captionData = {tremor, correlation, latency};
captionNames = {'Tremor', 'Correlation', 'Latency'};

for k = 1:numel(captionData)
    [captionP, captionANOVA] = anova1(captionData{k}, [], 'off');

    fprintf('%s: F(%d, %d) = %.10g; P = %.10g\n', ...
        captionNames{k}, ...
        captionANOVA{2,3}, captionANOVA{3,3}, ...
        captionANOVA{2,5}, captionP);
end