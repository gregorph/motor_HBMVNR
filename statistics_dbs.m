load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\dbs_statistics_variable.mat")

repetitions = 10000; %Define number of timeslices run

A = 0.3;          % Amplitude of sine wave
T = 10;        % Length of one cycle (period) in seconds
C = 0.5;          % Centre (vertical offset)
fs = 250;      % Sampling frequency (points per second)
L = repetitions/fs;          % Total signal length in seconds

% Derived values
f = 1/T;                    % Frequency (Hz)
ttime = 0:1/fs:L;               % Time vector

% Generate sine wave with offset
desiredposition(1,:) = C + A * sin(2*pi*f*ttime);%ones(1,repetitions+1)*0.5;

R = analyse_fig3b_dbs(fig3b_muscleposition, desiredposition(1,1:10000), ...
    'Fs', 250, ...
    'TargetFrequency', 0.1, ...
    'MaxLag', 500, ...
    'NumPermutations', 19999, ...
    'NumBootstraps', 20000, ...
    'RandomSeed', 1);

disp(R.tables.latencyPhaseSeconds)
disp(R.tables.trackingAmplitude)
disp(R.tables.trackingR2)

disp(R.tables.covariance)
% disp(R.tables.maxCovariance)
disp(R.tables.tremor)
%%

captionMetrics = {'tremor', 'covariance', ...
    'trackingAmplitude', 'trackingR2'};

captionRows = ismember(R.summary.Metric, captionMetrics);

disp(R.summary(captionRows, ...
    {'Metric', 'Condition', 'NPairs', ...
    'tStatistic', 'F', 'DF1', 'DF2', 'pBonferroni'}))

disp(R.phaseSummary(:, ...
    {'Condition', 'NPairs', 'D', 'pBonferroni'}))

