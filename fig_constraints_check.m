%% calculate motor trajectory changing pacemaker freq
clear all

circuits = 1000; %Number of circuits to run in parallel %1000

gravity_anti = 0; %Define external forces acting on joint. See eq. 20

neuron_fire = 0.80; %Define the default neuron excitability. In this case neuron has 85% chance firing given a single default excitatory release

resistance = gravity_anti;

repetitions = 10000;
A = 0.3;          % Amplitude of sine wave
T = 3;        % Length of one cycle (period) in seconds
C = 0.5;          % Centre (vertical offset)
fs = 250;      % Sampling frequency (points per second)
L = repetitions/fs;          % Total signal length in seconds

% Derived values
f = 1/T;                    % Frequency (Hz)
ttime = 0:1/fs:L;               % Time vector

% Generate sine wave with offset
desiredposition(1,:) = C + A * sin(2*pi*f*ttime);

mod_inhibition = -0.5; %Set strength of GABA synapses
mod_excitation = 0.5; %Set strength of Glutamate synapses

golgi_cellfire_flexor = zeros(100,1);
golgi_cellfire_extensor = zeros(100,1);
cerebellar_agonist = zeros(100,1);
cerebellar_antagonist = zeros(100,1);

muscleposition_str = 0.5;
gain = 1.5; %Define Strength of muscle
gain = gain/circuits;

n = repetitions;

idx1 = 1:5000; %0.1

resting_md = 0.02;
motor_drive = ones(1, n)*resting_md;
% motor_drive(idx1) = 0.4;
motor_drive(idx1+2500) = 0.4;%resting_md+0.1;

% Define segment indices

dopamine_levels = ones(1,n)*0.75;
dopamine_levels(idx1+2500) = 0.8;

hyperdirect_firing = zeros(1,n);
cerebellum_VL_weight = 0.025; %0.025 for PD
cerebellum_CM_weight = 0.05; %0.06 for PD
stn_bg_connectivity = 2; %2
bg_cer_connectivity = 1.8;
CM_striatum_w = 0.24;
CM_acth_modifier = 0.85; %1
ach_levels = 1.5;
L5_CM_same_weight = 0.5;
L5_CM_opp_weight = 0.2;
stn_fire_rate = 0.24;
cortex_acth_modifier = 0.75;
BG_CM_weight = 1.2;
BG_VL_weight = 3;
CM_diffuse = zeros(size(motor_drive,1),fix(0.200/(1/fs)));
% CM_cummulative = zeros(1,fix(0.200/(1/fs)));
CM_nucleus.agonist = zeros(1,fix(0.200/(1/fs)));
CM_nucleus.antagonist = zeros(1,fix(0.200/(1/fs)));
CM_agonist_keep = zeros(1,repetitions);
CM_antagonist_keep = zeros(1,repetitions);
CM_index = 1;
acth_rr_agonist = zeros(1,repetitions);
acth_rr_antagonist = zeros(1,repetitions);
gpi_agonist_keep = zeros(1,repetitions);
gpi_antagonist_keep = zeros(1,repetitions);

L5_flexor_prev = zeros(size(desiredposition,1),circuits);
L5_extensor_prev = zeros(size(desiredposition,1),circuits);

freqrange_check = 0:2:100;
num_repetitions = 48;

CM_nucleus_template.agonist = zeros(1,fix(0.200/(1/fs)));
CM_nucleus_template.antagonist = zeros(1,fix(0.200/(1/fs)));

results = struct( ...
    'muscleposition_keep', zeros(1, repetitions), ...
    'cwt_power', zeros(67, repetitions), ...
    'freqs', zeros(67, 1), ...
    'covariance', 0, ...
    'error', 0, ...
    'power_rest', 0, ...
    'peak_frequency_rest', 0 ...
    );

results = repmat(results, 6, length(freqrange_check), 100);

for condition_loop = 4:5%3
    weight_params.keepdefaults = true;
    pathologies.keepdefaults = true;
    weight_params.cerebellum_VL_weight = cerebellum_VL_weight;
    weight_params.cerebellum_CM_weight = cerebellum_CM_weight;
    weight_params.CM_acth_modifier = CM_acth_modifier;
    weight_params.CM_striatum_w = CM_striatum_w;
    weight_params.L5_CM_same_weight = L5_CM_same_weight;
    weight_params.L5_CM_opp_weight = L5_CM_opp_weight;
    weight_params.cortex_acth_modifier = cortex_acth_modifier;
    weight_params.BG_CM_weight = BG_CM_weight;
    weight_params.BG_VL_weight = BG_VL_weight;
    weight_params.stn_fire_rate = stn_fire_rate;

    parfor freq_check_index = 1:32%length(freqrange_check) %parfor
tic
        for tremor_repeat = 1:num_repetitions%parfor tremor_repeat = 1:10

            muscleposition = ones(1,repetitions)*muscleposition_str;
            golgi_cellfire_flexor = zeros(100,1);
            golgi_cellfire_extensor = zeros(100,1);
            cerebellar_agonist = zeros(1,100);
            cerebellar_antagonist = zeros(1,100);
            output = zeros(1,repetitions);
            output_mp = zeros(1,2);
            CM_index = 1;

            thalamus_flexor = zeros(1,circuits);
            thalamus_extensor = zeros(1,circuits);
            thalamus_net = zeros(1,repetitions);
            thalamus_summative = zeros(1,repetitions);
            velocity = zeros(1,repetitions);
            error_observed = zeros(1,repetitions);
            muscleposition_keep = zeros(1,repetitions);

            CM_nucleus = CM_nucleus_template;

            L5_flexor_prev = zeros(1,circuits);
            L5_extensor_prev = zeros(1,circuits);

            for i = 1:length(ttime)-1
                weight_internal = weight_params;
                if(condition_loop == 1)
                    weight_internal.stn_fire_rate = (1-0.8)/(1-(freqrange_check(freq_check_index)/250));
                elseif(condition_loop == 2)
                    weight_internal.gpe_fire_rate = (1-0.8)/(1-(freqrange_check(freq_check_index)/250));
                elseif(condition_loop == 3)%&&i>=repetitions/2)
                    weight_internal.gpi_fire_rate = (1-0.8)/(1-(freqrange_check(freq_check_index)/250));
                elseif(condition_loop == 4)
                    vertical_shift = -0.02:0.02:0.6;
                    min_freqprop = zeros(1,32);
                    max_freqprop = zeros(1,32);
                    min_freqprop(1:32) = 0.04;
                    max_freqprop(1:32) = 0.36
                elseif(condition_loop == 5)
                    vertical_shift = zeros(1,32);
                    min_freqprop = zeros(1,32);
                    max_freqprop = zeros(1,32);
                    min_freqprop(1:32) = 0.04;
                    max_freqprop(1:32) = 0.34:0.02:0.96;
                elseif(condition_loop == 6)
                    vertical_shift = zeros(1,32);
                    vertical_shift(1,1:32) = 0.6;
                    min_freqprop = zeros(1,32);
                    max_freqprop = zeros(1,32);
                    min_freqprop(1:32) = -0.3:0.02:0.32;
                    max_freqprop(1:32) = 0.36;
                end

                %% BG
                output = zeros(1,circuits);
                %[b_ganglia_dep_probabilities] = basal_ganglia_fcn_pts_rtl(motor_drive(motor_programs,i),2,dopamine_levels(i),weight_internal,pathologies,10,hyperdirect_firing(i),CM_diffuse(motor_programs,:),CM_index);
                [b_ganglia_dep_probabilities_agonist] = basal_ganglia_fcn_pts_rtl(motor_drive(1,i),ach_levels,dopamine_levels(i),weight_internal,pathologies,10,hyperdirect_firing(i),CM_nucleus.agonist,CM_index);
                [b_ganglia_dep_probabilities_antagonist] = basal_ganglia_fcn_pts_rtl(motor_drive(1,i),ach_levels,dopamine_levels(i),weight_internal,pathologies,10,hyperdirect_firing(i),CM_nucleus.antagonist,CM_index);

                %acth_rr_agonist(i) = b_ganglia_dep_probabilities_agonist.acth_release_rate;
                %acth_rr_antagonist(i) = b_ganglia_dep_probabilities_antagonist.acth_release_rate;

                %% Cerebellum
                % desiredposition(i) = centre+sin(i/(cycle/(2*pi)))*amplitude;

                if(condition_loop <=3)
                    desired_input = (0.04 + 0.32*desiredposition(1,i));
                    muscleposition_input = 0.04 + 0.32*muscleposition(i);
                    desired_input_rev = (0.04 + 0.32) - desired_input + 0.04; %corrected by adding + 0.04
                    muscleposition_input_rev = (0.04 + 0.32) - muscleposition_input + 0.04; %corrected by adding + 0.04
                else

                    desired_input = ((min_freqprop(freq_check_index) + vertical_shift(freq_check_index)) + (max_freqprop(freq_check_index) + vertical_shift(freq_check_index))*desiredposition(1,i));
                    muscleposition_input = (min_freqprop(freq_check_index) + vertical_shift(freq_check_index)) + (max_freqprop(freq_check_index) + vertical_shift(freq_check_index))*muscleposition(i);
                    desired_input_rev = ((min_freqprop(freq_check_index) + vertical_shift(freq_check_index)) + (max_freqprop(freq_check_index) + vertical_shift(freq_check_index))) - desired_input + (min_freqprop(freq_check_index) + vertical_shift(freq_check_index)); %corrected by adding + 0.04
                    muscleposition_input_rev = ((min_freqprop(freq_check_index) + vertical_shift(freq_check_index)) + (max_freqprop(freq_check_index) + vertical_shift(freq_check_index))) - muscleposition_input + (min_freqprop(freq_check_index) + vertical_shift(freq_check_index)); %corrected by adding + 0.04
                    
                    maxfreq_limit = max_freqprop(freq_check_index) + vertical_shift(freq_check_index);
                end
                desired_input = put_in_range(desired_input - bg_cer_connectivity*rand*b_ganglia_dep_probabilities_agonist.gpi_dep_probability,0,0.36);
                muscleposition_input = put_in_range(muscleposition_input + stn_bg_connectivity*rand*b_ganglia_dep_probabilities_antagonist.stn_dep_probability ...
                    -0.075*bg_cer_connectivity*rand*b_ganglia_dep_probabilities_antagonist.gpi_dep_probability,0,0.36);
                desired_input_rev = put_in_range(desired_input_rev - bg_cer_connectivity*rand*b_ganglia_dep_probabilities_antagonist.gpi_dep_probability,0,0.36);
                muscleposition_input_rev = put_in_range(muscleposition_input_rev + stn_bg_connectivity*rand*b_ganglia_dep_probabilities_agonist.stn_dep_probability...
                    -0.075*bg_cer_connectivity*rand*b_ganglia_dep_probabilities_agonist.gpi_dep_probability,0,0.36);

                if(motor_drive(i) == resting_md)
                    desired_input = 0;
                    desired_input_rev = 0;
                end

                for y = 1:100
                    [cerebellar_agonist(y),golgi_cellfire_flexor(y),~] = hbm_vnr(desired_input,muscleposition_input_rev,golgi_cellfire_flexor(y),mod_inhibition,mod_excitation,neuron_fire);

                    [cerebellar_antagonist(y),golgi_cellfire_extensor(y),~] = hbm_vnr(desired_input_rev,muscleposition_input,golgi_cellfire_extensor(y),mod_inhibition,mod_excitation,neuron_fire);
                end

                %% M1
                CM_index = mod(i,length(CM_diffuse)) + ((mod(i,length(CM_diffuse)) == 0)*length(CM_diffuse));
                % [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_diffuse(1,CM_index)] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,b_ganglia_dep_probabilities.gpi_firings,motor_drive(1,i),circuits,weight_params);
                [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_tmp] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,...
                    b_ganglia_dep_probabilities_agonist.gpi_firings,b_ganglia_dep_probabilities_antagonist.gpi_firings,...
                    motor_drive(1,i),circuits,L5_flexor_prev(1,:),L5_extensor_prev(1,:),weight_internal);

                L5_flexor_prev(1,:) = L5_flexor;
                L5_extensor_prev(1,:) = L5_extensor;
                % gpi_agonist_keep(i) = mean(b_ganglia_dep_probabilities_agonist.gpi_firings);
                % gpi_antagonist_keep(i) = mean(b_ganglia_dep_probabilities_antagonist.gpi_firings);

                CM_nucleus.agonist(CM_index) = CM_tmp.CM_agonist;
                % CM_agonist_keep(i) = CM_nucleus.agonist(CM_index);
                CM_nucleus.antagonist(CM_index) = CM_tmp.CM_antagonist;
                % CM_antagonist_keep(i) = CM_nucleus.antagonist(CM_index);
                CM_nucleus.greatest(CM_index) = max(CM_nucleus.agonist(CM_index),CM_nucleus.antagonist(CM_index));
                CM_nucleus.mean(CM_index) = (CM_nucleus.agonist(CM_index)+CM_nucleus.antagonist(CM_index))/2;
                %% Spinal Circuits

                spinal_circuit_output = spinal_circuits_fcn_pts(L5_flexor,L5_extensor,muscleposition(i),weight_internal);

                output(i) = sum(spinal_circuit_output.flexorpower*gain/10,"all")-sum(spinal_circuit_output.extensorpower*gain/10,"all"); %See Eq. 20

                % CM_activity_keep(condition_loop,i) = CM_nucleus.greatest(CM_index);

                muscleposition(i+1) = muscleposition(i) + output(i)-resistance; %See Eq. 21

                %Define boundaries of joint range of motion
                if(muscleposition(i+1)<0)
                    muscleposition(i+1) = 0;
                elseif(muscleposition(i+1)>1)
                    muscleposition(i+1) = 1;
                end
            end

            results(condition_loop,freq_check_index,tremor_repeat).muscleposition_keep = muscleposition(1,1:10000);

            % results(condition_loop,freq_check_index,tremor_repeat).covariance = min(cov(muscleposition(1,2500:7500),desiredposition(1,2500:7500)),[],'all');
            %
            % results(condition_loop,freq_check_index,tremor_repeat).error = mean(abs(muscleposition(1,2500:7500) - desiredposition(1,2500:7500)));
            %
            % [full_res,freqs] = cwt(muscleposition(1,1:10000),250,FrequencyLimits=[1 20],VoicesPerOctave=20,TimeBandwidth=120);
            % power_movement = abs(full_res);
            %
            % power_rest = mean([power_movement(:,1:2500) power_movement(:,7500:10000)],2);
            %
            % results(condition_loop,freq_check_index,tremor_repeat).peak_frequency_rest = sum(freqs .* power_rest) / sum(power_rest);
            %
            % results(condition_loop,freq_check_index,tremor_repeat).power_rest = sum([power_movement(:,1:2500) power_movement(:,7500:10000)],"all");
            %
            % results(condition_loop,freq_check_index,tremor_repeat).power_movement = sum(power_movement(:,2500:7500),"all");
            %
            % results(condition_loop,freq_check_index,tremor_repeat).cwt_power = power_movement;
        end
        toc
    end

    % output = movmean(output,10);
    % ax1 = subplot(3,1,1);
    % plot(ttime,muscleposition)
    % hold on
    % ylim([0.1 0.9])
    % ax2 = subplot(3,1,2);
    % % plot(movmean(acth_rr_agonist,100),'*')
    % % hold on
    % % plot(movmean(acth_rr_antagonist,100),'x')
    % tmp = max(CM_agonist_keep,CM_antagonist_keep);
    % plot(ttime(1,1:30000),movmean(tmp,100))
    % hold on
    % plot(ttime(1,1:30000),movmean(min(gpi_antagonist_keep,gpi_agonist_keep),100))
    % ax3 = subplot(3,1,3);
    % % plot(movmean(CM_agonist_keep,100),'*')
    % % hold on
    % % plot(movmean(CM_antagonist_keep,100),'x')
    % tmp = min(acth_rr_agonist,acth_rr_antagonist);
    % plot(ttime(1,1:30000),movmean(tmp,100))
    % linkaxes([ax1 ax2 ax3],'x')
    % xlim([0 ttime(1,10000)])
end

%% calculate freqs from motor trajectory
clear all
% load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\freq_tremor_pacemaker_27march.mat")

matObj = matfile("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\results_sensitivity_9september.mat");

repetitions = 10000;
A = 0.3;          % Amplitude of sine wave
T = 3;        % Length of one cycle (period) in seconds
C = 0.5;          % Centre (vertical offset)
fs = 250;      % Sampling frequency (points per second)
L = repetitions/fs;          % Total signal length in seconds

% Derived values
f = 1/T;                    % Frequency (Hz)
ttime = 0:1/fs:L;               % Time vector

% Generate sine wave with offset
desiredposition(1,:) = C + A * sin(2*pi*f*ttime);
for condition_loop = 4:5%3
    results = matObj.results(condition_loop,:,:);
    for freq_check_index = 1:32 %parfor
        for tremor_repeat = 1:48%parfor tremor_repeat = 1:10
            muscleposition = results(1,freq_check_index,tremor_repeat).muscleposition_keep;
            results(1,freq_check_index,tremor_repeat).covariance = min(cov(muscleposition(1,2500:7500),desiredposition(1,2500:7500)),[],'all');

            results(1,freq_check_index,tremor_repeat).error = mean(abs(muscleposition(1,2500:7500) - desiredposition(1,2500:7500)));

            [full_res,freqs] = cwt(muscleposition(1,1:10000),250,FrequencyLimits=[0.5 20],VoicesPerOctave=20,TimeBandwidth=120);
            
            f_new = 0.5:0.25:20;
            freqs = flipud(freqs);
            full_res = flipud(full_res);

            full_res_interp = interp1(freqs, full_res, f_new, 'makima');
            clear freqs
            freqs = f_new;

            power_movement = abs(full_res_interp);

            power_rest = mean([power_movement(:,1000:2000) power_movement(:,8000:9000)],2); %power_movement(:,7500:10000)],2);

            results(1,freq_check_index,tremor_repeat).peak_frequency_rest = sum(freqs.' .* power_rest) / sum(power_rest);

            [pks, locs] = findpeaks(movmean(power_rest,10), freqs,"MinPeakWidth",0.5);   % locs are frequencies
            [~, idx] = max(pks);              % tallest peak
            results(1,freq_check_index,tremor_repeat).tallestFreq = locs(idx);

            results(1,freq_check_index,tremor_repeat).power_rest = sum([power_movement(:,1:2500) power_movement(:,7500:10000)],"all");

            results(1,freq_check_index,tremor_repeat).power_movement = sum(power_movement(:,2500:7500),"all");



            % results(condition_loop,freq_check_index,tremor_repeat).cwt_power = power_movement;
        end
        results_avg(condition_loop,freq_check_index).covariance = mean( ...
            [results(1,freq_check_index,1:48).covariance], 'omitnan');
        results_avg(condition_loop,freq_check_index).covariance_std = std( ...
            [results(1,freq_check_index,1:48).covariance], 'omitnan');
        results_avg(condition_loop,freq_check_index).covariance_se = std( ...
            [results(1,freq_check_index,1:48).covariance], 'omitnan') / sqrt(48);

        results_avg(condition_loop,freq_check_index).error = mean( ...
            [results(1,freq_check_index,1:48).error], 'omitnan');
        results_avg(condition_loop,freq_check_index).error_std = std( ...
            [results(1,freq_check_index,1:48).error], 'omitnan');
        results_avg(condition_loop,freq_check_index).error_se = std( ...
            [results(1,freq_check_index,1:48).error], 'omitnan') / sqrt(48);

        results_avg(condition_loop,freq_check_index).peak_frequency_rest = mean( ...
            [results(1,freq_check_index,1:48).peak_frequency_rest], 'omitnan');
        results_avg(condition_loop,freq_check_index).peak_frequency_rest_std = std( ...
            [results(1,freq_check_index,1:48).peak_frequency_rest], 'omitnan');
        results_avg(condition_loop,freq_check_index).peak_frequency_rest_se = std( ...
            [results(1,freq_check_index,1:48).peak_frequency_rest], 'omitnan') / sqrt(48);

        results_avg(condition_loop,freq_check_index).tallest_frequency_rest = mean( ...
            [results(1,freq_check_index,1:48).tallestFreq], 'omitnan');
        results_avg(condition_loop,freq_check_index).tallest_frequency_rest_std = std( ...
            [results(1,freq_check_index,1:48).tallestFreq], 'omitnan');
        results_avg(condition_loop,freq_check_index).tallest_frequency_rest_se = std( ...
            [results(1,freq_check_index,1:48).tallestFreq], 'omitnan') / sqrt(48);

        results_avg(condition_loop,freq_check_index).power_rest = mean( ...
            [results(1,freq_check_index,1:48).power_rest], 'omitnan');
        results_avg(condition_loop,freq_check_index).power_rest_std = std( ...
            [results(1,freq_check_index,1:48).power_rest], 'omitnan');
        results_avg(condition_loop,freq_check_index).power_rest_se = std( ...
            [results(1,freq_check_index,1:48).power_rest], 'omitnan') / sqrt(48);

        results_avg(condition_loop,freq_check_index).power_movement = mean( ...
            [results(1,freq_check_index,1:48).power_movement], 'omitnan');
        results_avg(condition_loop,freq_check_index).power_movement_std = std( ...
            [results(1,freq_check_index,1:48).power_movement], 'omitnan');
        results_avg(condition_loop,freq_check_index).power_movement_se = std( ...
            [results(1,freq_check_index,1:48).power_movement], 'omitnan') / sqrt(48);

        results_avg(condition_loop,freq_check_index).tremor_mov_ratio = mean( ...
            ([results(1,freq_check_index,1:48).power_rest]./[results(1,freq_check_index,1:48).power_movement]), 'omitnan');
        results_avg(condition_loop,freq_check_index).tremor_mov_ratio_std = std( ...
            ([results(1,freq_check_index,1:48).power_rest]./[results(1,freq_check_index,1:48).power_movement]), 'omitnan');
        results_avg(condition_loop,freq_check_index).tremor_mov_ratio_se = std( ...
           ([results(1,freq_check_index,1:48).power_rest]./[results(1,freq_check_index,1:48).power_movement]), 'omitnan') / sqrt(48);

    end
    clear results
end

% figure
% plot([results_avg(1,:).tallest_frequency_rest])
% hold on
% plot([results_avg(2,:).tallest_frequency_rest])
% plot([results_avg(3,:).tallest_frequency_rest])


%% test
clear all
repetitions = 10000; %Define number of timeslices run

circuits = 500; %Number of circuits to run in parallel %1000

gravity_anti = 0; %Define external forces acting on joint. See eq. 20

neuron_fire = 0.80; %Define the default neuron excitability. In this case neuron has 85% chance firing given a single default excitatory release

resistance = gravity_anti;

A = 0.3;          % Amplitude of sine wave
T = 3;        % Length of one cycle (period) in seconds
C = 0.5;          % Centre (vertical offset)
fs = 250;      % Sampling frequency (points per second)
L = repetitions/fs;          % Total signal length in seconds

% Derived values
f = 1/T;                    % Frequency (Hz)
ttime = 0:1/fs:L;               % Time vector

% Generate sine wave with offset
desiredposition(1,:) = C + A * sin(2*pi*f*ttime);

mod_inhibition = -0.5; %Set strength of GABA synapses
mod_excitation = 0.5; %Set strength of Glutamate synapses

golgi_cellfire_flexor = zeros(100,1);
golgi_cellfire_extensor = zeros(100,1);
cerebellar_agonist = zeros(100,1);
cerebellar_antagonist = zeros(100,1);

muscleposition_str = 0.5;
gain = 1.5; %Define Strength of muscle
gain = gain/circuits;

n = repetitions;

idx1 = 1:5000; %0.1

resting_md = 0.02;
motor_drive = ones(1, n)*resting_md;
% motor_drive(idx1) = 0.4;
motor_drive(idx1+2500) = 0.4;%resting_md+0.1;

% Define segment indices

dopamine_levels = ones(1,n)*0.75;
dopamine_levels(idx1+2500) = 0.8;

hyperdirect_firing = zeros(1,n);
cerebellum_VL_weight = 0.025; %0.025 for PD
cerebellum_CM_weight = 0.05; %0.06 for PD
stn_bg_connectivity = 2; %2
bg_cer_connectivity = 1.8;
CM_striatum_w = 0.24;
CM_acth_modifier = 0.85; %1
ach_levels = 1.5;
L5_CM_same_weight = 0.5;
L5_CM_opp_weight = 0.2;
stn_fire_rate = 0.24;
cortex_acth_modifier = 0.75;
BG_CM_weight = 1.2;
BG_VL_weight = 3;
CM_diffuse = zeros(size(motor_drive,1),fix(0.200/(1/fs)));
% CM_cummulative = zeros(1,fix(0.200/(1/fs)));
CM_nucleus.agonist = zeros(1,fix(0.200/(1/fs)));
CM_nucleus.antagonist = zeros(1,fix(0.200/(1/fs)));
CM_agonist_keep = zeros(1,repetitions);
CM_antagonist_keep = zeros(1,repetitions);
CM_index = 1;
acth_rr_agonist = zeros(1,repetitions);
acth_rr_antagonist = zeros(1,repetitions);
gpi_agonist_keep = zeros(1,repetitions);
gpi_antagonist_keep = zeros(1,repetitions);

L5_flexor_prev = zeros(size(desiredposition,1),circuits);
L5_extensor_prev = zeros(size(desiredposition,1),circuits);

freqrange_check = 0:2:100;
num_repetitions = 48;

CM_nucleus_template.agonist = zeros(1,fix(0.200/(1/fs)));
CM_nucleus_template.antagonist = zeros(1,fix(0.200/(1/fs)));

results = struct( ...
    'muscleposition_keep', zeros(1, repetitions), ...
    'cwt_power', zeros(67, repetitions), ...
    'freqs', zeros(67, 1), ...
    'covariance', 0, ...
    'error', 0, ...
    'power_rest', 0, ...
    'peak_frequency_rest', 0 ...
    );

results = repmat(results, 3, length(freqrange_check), 100);

min_freqprop = -0.15;
max_freqprop = 0.36;
vertical_shift = 0.6; %increase vertical shift == downshift in frequency

for condition_loop = 1:1%3
    weight_params.keepdefaults = true;
    pathologies.keepdefaults = true;
    weight_params.cerebellum_VL_weight = cerebellum_VL_weight;
    weight_params.cerebellum_CM_weight = cerebellum_CM_weight;
    weight_params.CM_acth_modifier = CM_acth_modifier;
    weight_params.CM_striatum_w = CM_striatum_w;
    weight_params.L5_CM_same_weight = L5_CM_same_weight;
    weight_params.L5_CM_opp_weight = L5_CM_opp_weight;
    weight_params.cortex_acth_modifier = cortex_acth_modifier;
    weight_params.BG_CM_weight = BG_CM_weight;
    weight_params.BG_VL_weight = BG_VL_weight;
    weight_params.stn_fire_rate = stn_fire_rate;

    for freq_check_index = 1:1 %parfor

        for tremor_repeat = 1:1%parfor tremor_repeat = 1:10

            muscleposition = ones(1,repetitions)*muscleposition_str;
            golgi_cellfire_flexor = zeros(100,1);
            golgi_cellfire_extensor = zeros(100,1);
            cerebellar_agonist = zeros(1,100);
            cerebellar_antagonist = zeros(1,100);
            output = zeros(1,repetitions);
            output_mp = zeros(1,2);
            CM_index = 1;

            thalamus_flexor = zeros(1,circuits);
            thalamus_extensor = zeros(1,circuits);
            thalamus_net = zeros(1,repetitions);
            thalamus_summative = zeros(1,repetitions);
            velocity = zeros(1,repetitions);
            error_observed = zeros(1,repetitions);
            muscleposition_keep = zeros(1,repetitions);

            CM_nucleus = CM_nucleus_template;

            L5_flexor_prev = zeros(1,circuits);
            L5_extensor_prev = zeros(1,circuits);

            for i = 1:length(ttime)-1
                weight_internal = weight_params;
                % if(condition_loop == 1)
                %     weight_internal.stn_fire_rate = (1-0.8)/(1-(freqrange_check(freq_check_index)/250));
                % elseif(condition_loop == 2)
                %     weight_internal.gpe_fire_rate = (1-0.8)/(1-(freqrange_check(freq_check_index)/250));
                % elseif(condition_loop == 3)%&&i>=repetitions/2)
                %     weight_internal.gpi_fire_rate = (1-0.8)/(1-(freqrange_check(freq_check_index)/250));
                % end

                %% BG
                output = zeros(1,circuits);
                %[b_ganglia_dep_probabilities] = basal_ganglia_fcn_pts_rtl(motor_drive(motor_programs,i),2,dopamine_levels(i),weight_internal,pathologies,10,hyperdirect_firing(i),CM_diffuse(motor_programs,:),CM_index);
                [b_ganglia_dep_probabilities_agonist] = basal_ganglia_fcn_dbs(motor_drive(1,i),ach_levels,dopamine_levels(i),weight_internal,pathologies,10,hyperdirect_firing(i),CM_nucleus.agonist,CM_index);
                [b_ganglia_dep_probabilities_antagonist] = basal_ganglia_fcn_dbs(motor_drive(1,i),ach_levels,dopamine_levels(i),weight_internal,pathologies,10,hyperdirect_firing(i),CM_nucleus.antagonist,CM_index);

                %acth_rr_agonist(i) = b_ganglia_dep_probabilities_agonist.acth_release_rate;
                %acth_rr_antagonist(i) = b_ganglia_dep_probabilities_antagonist.acth_release_rate;

                %% Cerebellum
                % desiredposition(i) = centre+sin(i/(cycle/(2*pi)))*amplitude;

                desired_input = ((min_freqprop + vertical_shift) + (max_freqprop + vertical_shift)*desiredposition(1,i));
                muscleposition_input = (min_freqprop + vertical_shift) + (max_freqprop + vertical_shift)*muscleposition(i);
                desired_input_rev = ((min_freqprop + vertical_shift) + (max_freqprop + vertical_shift)) - desired_input + (min_freqprop + vertical_shift); %corrected by adding + 0.04
                muscleposition_input_rev = ((min_freqprop + vertical_shift) + (max_freqprop + vertical_shift)) - muscleposition_input + (min_freqprop + vertical_shift); %corrected by adding + 0.04

                desired_input = put_in_range(desired_input - bg_cer_connectivity*rand*b_ganglia_dep_probabilities_agonist.gpi_dep_probability,0,0.36);
                muscleposition_input = put_in_range(muscleposition_input + stn_bg_connectivity*rand*b_ganglia_dep_probabilities_antagonist.stn_dep_probability ...
                    -0.075*bg_cer_connectivity*rand*b_ganglia_dep_probabilities_antagonist.gpi_dep_probability,0,0.36);
                desired_input_rev = put_in_range(desired_input_rev - bg_cer_connectivity*rand*b_ganglia_dep_probabilities_antagonist.gpi_dep_probability,0,0.36);
                muscleposition_input_rev = put_in_range(muscleposition_input_rev + stn_bg_connectivity*rand*b_ganglia_dep_probabilities_agonist.stn_dep_probability...
                    -0.075*bg_cer_connectivity*rand*b_ganglia_dep_probabilities_agonist.gpi_dep_probability,0,0.36);

                if(motor_drive(i) == resting_md)
                    desired_input = 0;
                    desired_input_rev = 0;
                end

                for y = 1:100
                    [cerebellar_agonist(y),golgi_cellfire_flexor(y),~] = hbm_vnr(desired_input,muscleposition_input_rev,golgi_cellfire_flexor(y),mod_inhibition,mod_excitation,neuron_fire);

                    [cerebellar_antagonist(y),golgi_cellfire_extensor(y),~] = hbm_vnr(desired_input_rev,muscleposition_input,golgi_cellfire_extensor(y),mod_inhibition,mod_excitation,neuron_fire);
                end

                %% M1
                CM_index = mod(i,length(CM_diffuse)) + ((mod(i,length(CM_diffuse)) == 0)*length(CM_diffuse));
                % [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_diffuse(1,CM_index)] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,b_ganglia_dep_probabilities.gpi_firings,motor_drive(1,i),circuits,weight_params);
                [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_tmp] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,...
                    b_ganglia_dep_probabilities_agonist.gpi_firings,b_ganglia_dep_probabilities_antagonist.gpi_firings,...
                    motor_drive(1,i),circuits,L5_flexor_prev(1,:),L5_extensor_prev(1,:),weight_internal);

                L5_flexor_prev(1,:) = L5_flexor;
                L5_extensor_prev(1,:) = L5_extensor;
                % gpi_agonist_keep(i) = mean(b_ganglia_dep_probabilities_agonist.gpi_firings);
                % gpi_antagonist_keep(i) = mean(b_ganglia_dep_probabilities_antagonist.gpi_firings);

                CM_nucleus.agonist(CM_index) = CM_tmp.CM_agonist;
                % CM_agonist_keep(i) = CM_nucleus.agonist(CM_index);
                CM_nucleus.antagonist(CM_index) = CM_tmp.CM_antagonist;
                % CM_antagonist_keep(i) = CM_nucleus.antagonist(CM_index);
                CM_nucleus.greatest(CM_index) = max(CM_nucleus.agonist(CM_index),CM_nucleus.antagonist(CM_index));
                CM_nucleus.mean(CM_index) = (CM_nucleus.agonist(CM_index)+CM_nucleus.antagonist(CM_index))/2;
                %% Spinal Circuits

                spinal_circuit_output = spinal_circuits_fcn_pts(L5_flexor,L5_extensor,muscleposition(i),weight_internal);

                output(i) = sum(spinal_circuit_output.flexorpower*gain/10,"all")-sum(spinal_circuit_output.extensorpower*gain/10,"all"); %See Eq. 20

                % CM_activity_keep(condition_loop,i) = CM_nucleus.greatest(CM_index);

                muscleposition(i+1) = muscleposition(i) + output(i)-resistance; %See Eq. 21

                %Define boundaries of joint range of motion
                if(muscleposition(i+1)<0)
                    muscleposition(i+1) = 0;
                elseif(muscleposition(i+1)>1)
                    muscleposition(i+1) = 1;
                end
            end
        end
    end
end
    [full_res,freqs] = cwt(muscleposition(1,1:5000),250);
    tst = mean(abs(full_res(:,1000:2000)),2);
    plot(freqs,tst./max(tst))
    hold on

%% figure generator
%load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\results_fig_brain_pacemaker_9sep.mat")
load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\results_fig_brain_pacemaker_29march.mat")
freqs = 0:2:100;
vertical_shift = -0.02:0.02:0.6;


subplot(2,2,1)
hold on
errorbar(vertical_shift*250,movmean([results_avg(4,2:end).tallest_frequency_rest],1)...
    ,movmean(1.96*[results_avg(4,2:end).tallest_frequency_rest_se],1),"LineStyle","--","Marker","x","LineWidth",1 ...
    ,"Color","k")
ylabel("Frequency (Hz)")
xlabel("Horizontal shift (Hz)")
subtitle("Tremor frequency")
ylim([3 10])

subplot(2,2,2)
hold on
errorbar(vertical_shift*250,movmean([results_avg(5,:).tallest_frequency_rest],10)...
    ,movmean(1.96*[results_avg(5,:).tallest_frequency_rest_se],10),"LineStyle","--","Marker","x","LineWidth",1 ...
    ,"Color","k")
xlabel("Frequency range (Hz)")
subtitle("Tremor frequency")
title("Tremor frequency against proprioceptor frequency changes")
ylim([3 10])

subplot(2,3,4)
errorbar(freqs,movmean([results_avg(1,:).tallest_frequency_rest],10)...
    ,movmean(1.96*[results_avg(1,:).tallest_frequency_rest_se],10),"LineStyle","--","Marker","*","LineWidth",1)
hold on
errorbar(freqs,movmean([results_avg(3,:).tallest_frequency_rest],10)...
    ,movmean(1.96*[results_avg(3,:).tallest_frequency_rest_se],10),"LineStyle","--","Marker","x","LineWidth",1)
xlabel("Pacemaker frequency (Hz)")
ylabel("Frequency (Hz)")
title("Tremor frequency")

subplot(2,3,5)
errorbar(freqs,movmean([results_avg(1,:).error],10)...
    ,movmean(1.96*[results_avg(1,:).error_se],10),"LineStyle","--","Marker","*","LineWidth",1)
hold on
errorbar(freqs,movmean([results_avg(3,:).error],10)...
    ,movmean(1.96*[results_avg(3,:).error_se],10),"LineStyle","--","Marker","x","LineWidth",1)
xlabel("Pacemaker frequency (Hz)")
ylabel("Error (units)")
title("Error")

subplot(2,3,6)
errorbar(freqs,movmean([results_avg(1,:).tremor_mov_ratio],10),1.96*movmean([results_avg(1,:).tremor_mov_ratio_se],10) ...
    ,"LineStyle","--","Marker","*","LineWidth",1)
hold on
errorbar(freqs,movmean([results_avg(3,:).tremor_mov_ratio],10),1.96*movmean([results_avg(3,:).tremor_mov_ratio_se],10) ...
    ,"LineStyle","--","Marker","x","LineWidth",1)
xlabel("Pacemaker frequency (Hz)")
title("Strength tremor over voluntary movement")
ylabel("Tremor to movement ratio")



