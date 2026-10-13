%% part 1
clear all

num_circuits_bg = 10000;
motor_drive = 0.01:0.005:0.4;
motor_drive_frequency = motor_drive*250;
dopamine_levels = 0:0.025:3;
dop_acth = dopamine_levels/2;

parkinsonian_level = [1 0.5 0.5 0.5 0.5 0.5];
levodopa_addition = [0 0 0.35 0.65 0 0];
anticholinergic = [1 1 1 1 0.6 1];
amantadine = [0.5 0.5 0.5 0.5 0.5 0.25];
stn_gpi_w = [1 1 1 1 1 0.5];
titles_top_row = {"Control","PD","Levodopa","Levodopa High dose","Anticholinergic","Amantadine"};

for condition_index = 1:length(amantadine)
    tic
    weight_params.keepdefaults = true;
    weight_params.glutamine_mod = amantadine(condition_index);
    weight_params.stn_gpi_w = stn_gpi_w(condition_index);
    pathologies.keepdefaults = true;
    for dop_index = 1:length(dopamine_levels)
        
        parfor motor_index = 1:length(motor_drive)

            [b_ganglia_dep_probabilities] = basal_ganglia_fcn_pts(motor_drive(motor_index),2*anticholinergic(condition_index),(dopamine_levels(dop_index)*parkinsonian_level(condition_index))+levodopa_addition(condition_index),weight_params,pathologies,num_circuits_bg);

            gpi_dp(dop_index,motor_index) = b_ganglia_dep_probabilities.gpi_dep_probability*250;
        end
        condition_checked(condition_index).gpi_dp(dop_index,:) = gpi_dp(dop_index,:);
    end
    toc
end
% fig generation
figure
for i = 1:length(condition_checked)
    ax(i) = subplot(1,length(condition_checked),i);
    mxfreq(i) = max(condition_checked(i).gpi_dp,[],"all");
    minfreq(i) = min(condition_checked(i).gpi_dp,[],"all");
    contourf(dop_acth,motor_drive_frequency,movmean(condition_checked(i).gpi_dp.',10),20)
    hold on;
    if(i==3)
    xlabel('Dopamine Levels (scaled)');
    end
    if(i == 1)
    ylabel('Motor Drive Frequency (Hz)');
    end
    title(titles_top_row{i});
end

for i = 1:length(mxfreq)
    clim(ax(i),[min(minfreq),max(mxfreq)])
end

%% part 2
clear desiredposition muscleposition output
neuron_fire = 0.80; %Define the default neuron excitability. In this case neuron has 85% chance firing given a single default excitatory release

repetitions = 1500; %Define number of timeslices run

circuits = 1000; %Number of circuits to run in parallel

gravity_anti = 0; %Define external forces acting on joint. See eq. 20

extensor_excitability = 1; %Change strength of extensor muscle
flexor_excitability = 1; %Change strength of extensor muscle

resistance = gravity_anti;

avwindow = 5; %window for moving average

mod_inhibition = -0.5; %Set strength of GABA synapses
mod_excitation = 0.5; %Set strength of Glutamate synapses

muscleposition(1) = 0;
gain = 1.5; %Define Strength of muscle
gain = gain/circuits;

A = 0.2;          % Amplitude of sine wave
T = 10;        % Length of one cycle (period) in seconds
C = 0.5;          % Centre (vertical offset)
fs = 250;      % Sampling frequency (points per second)
L = repetitions/fs;          % Total signal length in seconds

% Derived values
f = 1/T;                    % Frequency (Hz)
ttime = 0:1/fs:L;               % Time vector

% Generate sine wave with offset
desiredposition(1,:) = C + A * sin(2*pi*f*ttime);
A = 0.15;
desiredposition(2,:) = C + A * sin(2*pi*f*ttime+(pi*0.5));

n = length(ttime);
motor_drive = zeros(2, n);
motor_drive(1,:) = 0.4; 

idx2 = floor(0.5*n)+1 : n;

% Assign values
motor_drive(2, idx2) = 0.1;   % Row 3, next 20%

% Define segment indices

% Assign values

% cerebellum_VL_weight = [0.0135 0.03 0.03 0.03 0.03];
% bg_cer_connectivity = [1.2 0.7 0.7 0.7 0.7];
cerebellum_VL_weight = [0.014 0.016 0.016 0.016 0.016];
bg_cer_connectivity = [1 0.6 0.6 0.6 0.6];
stn_bg_connectivity = [0.025 0.1 0.1 0.1 0.1];

parkinsonian_level = [1 0.4 0.4 0.4 0.4];
levodopa_addition = [0 0 0.3 0 0];
anticholinergic = [1 1 1 0.6 1];
amantadine = [0.5 0.5 0.5 0.5 0.25];
stn_gpi_w = [1 1 1 1 0.2];


for condition_loop = 1:length(amantadine)
    weight_params.keepdefaults = true;
    weight_params.glutamine_mod = amantadine(condition_loop);
    weight_params.stn_gpi_w = stn_gpi_w(condition_loop);
    pathologies.keepdefaults = true;
    weight_params.cerebellum_VL_weight = cerebellum_VL_weight(condition_loop);
    fig3b = struct;
    num_parts_graph = 48;
    parfor cor_index = 1:num_parts_graph
        golgi_cellfire_flexor = zeros(100,1);
        golgi_cellfire_extensor = zeros(100,1);
        cerebellar_agonist = zeros(1,100);
        cerebellar_antagonist = zeros(1,100);
        muscleposition = ones(1,repetitions+1)*0.5;
        output = zeros(1,repetitions);
        output_mp = zeros(1,2);
        
        thalamus_flexor_collective = zeros(1,repetitions);
        thalamus_extensor_collective = zeros(1,repetitions);
        thalamus_net = zeros(1,repetitions);
        thalamus_summative = zeros(1,repetitions);
        velocity = zeros(1,repetitions);
        error_observed = zeros(1,repetitions);

        tic
        for i = 1:length(ttime)-1
            %% BG
            thalamus_flexor = zeros(2,circuits);
            thalamus_extensor = zeros(2,circuits);
            for motor_programs = 1:2
                
                [b_ganglia_dep_probabilities] = basal_ganglia_fcn_pts(motor_drive(motor_programs,i),2*anticholinergic(condition_loop),2*parkinsonian_level(condition_loop)+levodopa_addition(condition_loop),weight_params,pathologies,10);

                %% Cerebellum
                % desiredposition(i) = centre+sin(i/(cycle/(2*pi)))*amplitude;
                % desiredposition_input = put_in_range(desiredposition(motor_programs,i) - rand*b_ganglia_dep_probabilities.gpi_dep_probability*bg_cer_connectivity(condition_loop),0,1); %normal weight 1.5 for resting tremor PD 1
                % 
                % muscleposition_tmp = put_in_range(muscleposition(i) + rand*b_ganglia_dep_probabilities.gpi_dep_probability*bg_cer_connectivity(condition_loop),0,1);


                desired_input = (0.04 + 0.32*desiredposition(motor_programs,i));
                muscleposition_input = 0.04 + 0.32*muscleposition(i);
                desired_input_rev = (0.04 + 0.32) - desired_input;
                muscleposition_input_rev = (0.04 + 0.32) - muscleposition_input;

                desired_input = put_in_range(desired_input - bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities.gpi_dep_probability,0,0.36);
                muscleposition_input = put_in_range(muscleposition_input + stn_bg_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities.stn_dep_probability - 0.1*bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities.gpi_dep_probability,0,0.36);
                desired_input_rev = put_in_range(desired_input_rev - bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities.gpi_dep_probability,0,0.36);
                muscleposition_input_rev = put_in_range(muscleposition_input_rev + stn_bg_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities.stn_dep_probability-0.1*bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities.gpi_dep_probability,0,0.36);


                for y = 1:100
                    [cerebellar_agonist(y),golgi_cellfire_flexor(y),~] = hbm_vnr(desired_input,muscleposition_input_rev,golgi_cellfire_flexor(y),mod_inhibition,mod_excitation,neuron_fire);

                    [cerebellar_antagonist(y),golgi_cellfire_extensor(y),~] = hbm_vnr(desired_input_rev,muscleposition_input,golgi_cellfire_extensor(y),mod_inhibition,mod_excitation,neuron_fire);
                end

                %% M1
                [L5_flexor,L5_extensor,thalamus_flexor(motor_programs,:),thalamus_extensor(motor_programs,:)] = M1_fcn_pts(cerebellar_agonist,cerebellar_antagonist,b_ganglia_dep_probabilities.gpi_firings,motor_drive(motor_programs,i),circuits,weight_params);
                %% Spinal Circuits

                spinal_circuit_output = spinal_circuits_fcn_pts(L5_flexor,L5_extensor,muscleposition(i),weight_params);

                output_mp(motor_programs) = sum(spinal_circuit_output.flexorpower*gain/10,"all")-sum(spinal_circuit_output.extensorpower*gain/10,"all"); %See Eq. 20

            end

            output(i) = sum(output_mp);
            muscleposition(i+1) = muscleposition(i) + output(i)-resistance; %See Eq. 21

            thalamus_flexor_collective(1,i) = sum(thalamus_flexor,"all");
            thalamus_extensor_collective(1,i) = sum(thalamus_extensor,"all");
            thalamus_net(1,i) = thalamus_flexor_collective(1,i) - thalamus_extensor_collective(1,i);
            thalamus_summative(1,i) = thalamus_flexor_collective(1,i) + thalamus_extensor_collective(1,i);
            velocity(1,i) = output(i);
            error_observed(1,i) = desiredposition(1,i) - muscleposition(i);

            %Define boundaries of joint range of motion
            if(muscleposition(i+1)<0)
                muscleposition(i1) = 0;
            elseif(muscleposition(i+1)>1)
                muscleposition(i+1) = 1;
            end
        end
        [tmp_power,freqs] = cwt(muscleposition,fs);
        tmp_power = abs(tmp_power);
        [~,index_freqlow] = min(abs(freqs-2.5));
        [~,index_freqhigh] = min(abs(freqs-8.5));
        power_motor_tremor1_fig3(condition_loop,cor_index) = mean(tmp_power(index_freqhigh:index_freqlow,1:repetitions/2),"all");
        power_motor_tremor2_fig3(condition_loop,cor_index) = mean(tmp_power(index_freqhigh:index_freqlow,1+(repetitions/2):end),"all");
        tremor_var(condition_loop,cor_index) = mean(movvar(muscleposition,100),"all");

        [tmp_power,freqs] = cwt(thalamus_summative,fs);
        tmp_power = abs(tmp_power);
        [~,index_freqlow] = min(abs(freqs-2.5));
        [~,index_freqhigh] = min(abs(freqs-8.5));

        power_thalamustremor1_fig3(condition_loop,cor_index) = mean(tmp_power(index_freqhigh:index_freqlow,1:repetitions/2),"all");
        power_thalamustremor2_fig3(condition_loop,cor_index) = mean(tmp_power(index_freqhigh:index_freqlow,1+(repetitions/2):end),"all");
        [~,index_freqlow] = min(abs(freqs-12.5));
        [~,index_freqhigh] = min(abs(freqs-30.5));
        power_thalamusbeta1_fig3(condition_loop,cor_index) = mean(tmp_power(index_freqhigh:index_freqlow,1:repetitions/2),"all");
        power_thalamusbeta2_fig3(condition_loop,cor_index) = mean(tmp_power(index_freqhigh:index_freqlow,1+(repetitions/2):end),"all");
        

        %correlation_keep(condition_loop,cor_index) = min(corrcoef(muscleposition(1+(repetitions/2):end),desiredposition(1,1+(repetitions/2):end)),[],'all');
        tmp_cov = cov(muscleposition,desiredposition(1,:));
        correlation_keep(condition_loop,cor_index) = tmp_cov(1,2);

        [r,lags] = xcorr(thalamus_net,error_observed(1,:));

        r(lags>=0) = 0;
        [max_corthalamus(condition_loop,cor_index),max_index_tmp] = max(r,[],"all");
        cor_thalamus_error_centre(condition_loop,cor_index) = r(repetitions/2);
        lag_keep(condition_loop,cor_index) = lags(max_index_tmp);
        corr_thalamus_error(condition_loop,cor_index) = min(corrcoef(thalamus_net(1+(repetitions/2):end),error_observed(1,1+(repetitions/2):end)),[],'all');
        corr_thalamus_desired(condition_loop,cor_index) = min(corrcoef(thalamus_net(1+(repetitions/2):end),desiredposition(1,1+(repetitions/2):end-1)),[],'all');

        [r,lags] = xcorr(muscleposition,desiredposition(1,:));
        r(lags>=0) = 0;
        [~,max_index_tmp] = max(r,[],"all");
        lag_keep_movement_desired(condition_loop,cor_index) = lags(max_index_tmp);

        fig3b_thalamus_flexor_collective(condition_loop,cor_index,:) = thalamus_flexor_collective;
        fig3b_thalamus_extensor_collective(condition_loop,cor_index,:) = thalamus_extensor_collective;
        fig3b_thalamus_net(condition_loop,cor_index,:) = thalamus_net;
        fig3b_thalamus_summative(condition_loop,cor_index,:) = thalamus_summative;
        fig3b_velocity(condition_loop,cor_index,:) = velocity;
        fig3b_error_observed(condition_loop,cor_index,:) = error_observed;


        toc
    end
end
%% ploting
% load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\fig_phamacological.mat")
load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\fig_phamacological_top_hd.mat")
num_circuits_bg = 1000;
motor_drive = 0.01:0.005:0.4;
motor_drive_frequency = motor_drive*250;
dopamine_levels = 0:0.025:3;
dop_acth = dopamine_levels/2;
titles_part1 = ["Healthy" "PD" "Levodopa" "Levodopa ++" "Anticholinergic" "Amantadine"];
fig = figure;
for i = 1:length(condition_checked)
    ax(i) = subplot(3,length(condition_checked),i);
    mxfreq(i) = max(condition_checked(i).gpi_dp,[],"all");
    minfreq(i) = min(condition_checked(i).gpi_dp,[],"all");
    contourf(dop_acth,motor_drive_frequency,movmean(condition_checked(i).gpi_dp.',40),20)
    xlim([0 1.5])
    hold on;
    xlabel('Dopamine - ACh release rate');
    if(i == 1)
        ylabel('Motor Drive Frequency (Hz)');
    end
    title(titles_part1(i));
end

for i = 1:length(mxfreq)
    clim(ax(i),[min(minfreq),max(mxfreq)])
end

subplot(3,2,3)
cmap = lines(5);
hold on
for g = 1:5
    velocity_plot(:,:) = abs(fig3b_velocity(g,:,:)*250);
    power_motor_tremor2_fig3(g,:) = (power_motor_tremor1_fig3(g,:))./(mean(velocity_plot,2)).';

    scatter(mean(velocity_plot,2), power_motor_tremor2_fig3(g,:), 30, ...
        'MarkerFaceColor', cmap(g,:), ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha',0.5);
    plotEllipse(mean(velocity_plot,2).', power_motor_tremor2_fig3(g,:), cmap(g,:));
    xlabel("Velocity (units/s)")
    ylabel("Power tremor / velocity")
    
end
title("Tremor against Velocity")
legend(["Healthy", '' "PD", '', "Levodopa", '', "Anticholinergic",'' , "Amantadine",''])

subplot(3,2,4)
cmap = lines(5);
hold on
for g = 1:5
    error_plot(:,:) = abs(fig3b_error_observed(g,:,:));
    scatter(abs(correlation_keep(g,:)), mean(error_plot,2), 30, ...
        'MarkerFaceColor', cmap(g,:), ...
        'MarkerEdgeColor', 'k', 'MarkerFaceAlpha',0.5);
    plotEllipse(abs(correlation_keep(g,:)), mean(error_plot,2).', cmap(g,:));
    xlabel("Covariance")
    ylabel("Error (units)")
end
ylim([0 1])
%xlim([0 11])
title("Error against Correlation")

subplot(3,3,7)
tremor_movement_ratio = power_thalamustremor2_fig3./power_thalamusbeta2_fig3;
violinplot(tremor_movement_ratio.')
hold on
boxchart(tremor_movement_ratio.')
title("Thalamus Tremor - Movement frequency ratio")
xticklabels(["Healthy" "PD" "Levodopa" "Anticholinergic" "Amantadine"])


subplot(3,3,8)
corr_thalamus_error(corr_thalamus_error<0) = 0;
violinplot(corr_thalamus_error.')
hold on
boxchart(corr_thalamus_error.')
ylabel("Correlation coeficient")
title("Correlation thalamus to error")
xticklabels(["Healthy" "PD" "Levodopa" "Anticholinergic" "Amantadine"])

subplot(3,3,9)
% lag_keep(lag_keep == 500) = 0;
violinplot(abs((lag_keep.')/250))
hold on
boxchart(abs((lag_keep.')/250))
title("Latency thalamus")
ylabel("Latency (s)");
xticklabels(["Healthy" "PD" "Levodopa" "Anticholinergic" "Amantadine"])