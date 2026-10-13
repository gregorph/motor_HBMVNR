%% part 1
clear all

num_circuits_bg = 10000;
motor_drive = 0.01:0.005:0.4;%0.01:0.005:0.4;
motor_drive_frequency = motor_drive*250;
dopamine_levels = 0:0.025:3;
dop_acth = dopamine_levels/2;

parkinsonian_level = [1 0.5 0.5 0.5];
levodopa_addition = [0 0 0 0];
anticholinergic = [1 1 1 1];
amantadine = [0.5 0.5 0.5 0.5];
CM_nucleus_agonist = zeros(1,fix(0.200/(1/250)));

stn_fire_rate = [0.3 0.3 0.3 0.3];

gpi_dbs = [0 0 0 0.125];
stn_dbs = [0 0 0.4 0];

dbs_rate_simulate = 0:0.025:1;


for condition_index = 2:length(amantadine)
    for dbs_stim = 1:length(dbs_rate_simulate)
        weight_params.keepdefaults = true;
        weight_params.stn_fire_rate = stn_fire_rate(condition_index);
        if(condition_index == 4)
            weight_params.dbs_gpi_a = dbs_rate_simulate(dbs_stim);
            weight_params.dbs_gpi_f = 130;
        elseif(condition_index == 3)
            weight_params.dbs_stn_a = dbs_rate_simulate(dbs_stim);
            weight_params.dbs_stn_f = 130;
        end
        pathologies.keepdefaults = true;
        for dop_index = 1:length(dopamine_levels)
            parfor motor_index = 1:length(motor_drive)

                [b_ganglia_dep_probabilities] = basal_ganglia_fcn_dbs(motor_drive(motor_index),2*anticholinergic(condition_index),(dopamine_levels(dop_index)*parkinsonian_level(condition_index)),weight_params,pathologies,num_circuits_bg,0,CM_nucleus_agonist,1);

                gpi_dp(dbs_stim,dop_index,motor_index) = b_ganglia_dep_probabilities.gpi_dep_probability*250;
                stn_dp(dbs_stim,dop_index,motor_index) = b_ganglia_dep_probabilities.stn_dep_probability*250;
                % gpe_dp(dbs_stim,dop_index,motor_index) = b_ganglia_dep_probabilities.gpe_dep_probability*250;

            end
            condition_checked(condition_index).gpi_dp(dbs_stim,dop_index,:) = gpi_dp(dbs_stim,dop_index,:);
            condition_checked(condition_index).stn_dp(dbs_stim,dop_index,:) = stn_dp(dbs_stim,dop_index,:);
        end
        % % if(condition_index<=2)
        % %     break;
        % % end
        disp(dbs_stim)
    end
end

figure
% tmp(:,:) = condition_checked(3).gpi_dp(:,1,:);
tmp(:,:) = gpi_dp(20,:,:);
contourf(movmean(tmp,10).',20)

for index_dif = 1:2
    condition_checked(4+index_dif).gpi_dp = movmean(condition_checked(2).gpi_dp(1,:,:),50,1) - movmean(condition_checked(2+index_dif).gpi_dp(:,:,:),50,1);
    mean_reduction(index_dif,:) = mean(condition_checked(4+index_dif).gpi_dp,[2 3]).';
    range_reduction(index_dif,:) = 1.96*std(condition_checked(4+index_dif).gpi_dp,0,[2 3]).';
    min_activation(index_dif,:) = min(condition_checked(2+index_dif).gpi_dp,[],[2 3]);
    max_activation(index_dif,:) = max(condition_checked(2+index_dif).gpi_dp,[],[2 3]);
    mean_activation(index_dif,:) = (min_activation(index_dif,:) + max_activation(index_dif,:))./2;
end

% fig generation
figure
plot_sequence = [3 4 0 5 6];
index_plot = [20 5 0 20 5];
counter1 = 1;
counter2 = 1;
for i = 1:length(plot_sequence)
    if(plot_sequence(i)~=0)
        ax(i) = subplot(1,length(plot_sequence),i);
        toplot(:,:) = condition_checked(plot_sequence(i)).gpi_dp(index_plot(i),:,:);
        if(plot_sequence(i)<=4)
            mxfreq(counter1) = max(condition_checked(plot_sequence(i)).gpi_dp,[],"all");
            minfreq(counter1) = min(condition_checked(plot_sequence(i)).gpi_dp,[],"all");
            counter1 = counter1 + 1;
            contourf(dop_acth,motor_drive_frequency,movmean(toplot.',20),20)
        else
            mxfreq2(counter2) = max(condition_checked(plot_sequence(i)).gpi_dp,[],"all");
            minfreq2(counter2) = min(condition_checked(plot_sequence(i)).gpi_dp,[],"all");
            counter2  = counter2 + 1;
            contourf(dop_acth,motor_drive_frequency,movmean(toplot.',50),10)
        end
        hold on;
        xlabel('Dop - ACh');
        ylabel('Motor Drive Frequency (Hz)');
    end
end

for i = 1:length(mxfreq)
    clim(ax(i),[min(minfreq),max(mxfreq)])
    %clim(ax(2+(i-1)*3),[min(minfreq2),max(mxfreq2)])
end



ax(3) = subplot(1,length(plot_sequence),3);
colors_plot = {[214,124,48],[62,150,81]};
markers_use = ['*' 'x'];
for i = 1:size(mean_reduction,1)
    % X=[dbs_rate_simulate,fliplr(dbs_rate_simulate)];                %#create continuous x value array for plotting
    % Y=[min_activation(i,:),fliplr(max_activation(i,:))];              %#create y values for out and then back
    % fill(X,Y,'b',"FaceAlpha",0.2);
    hold on
    yyaxis left
    plot(dbs_rate_simulate,movmean(mean_reduction(i,:),10),"Color",colors_plot{1}./250,Marker=markers_use(i)...
        ,LineStyle="--",MarkerIndices=1:10:length(dbs_rate_simulate),MarkerSize=15,LineWidth=1.5)
    ylabel("Mean GPi F reduction")
    yyaxis right
    plot(dbs_rate_simulate,movmean(max_activation(i,:)-min_activation(i,:),10),"Color",colors_plot{2}./250,Marker=markers_use(i)...
        ,LineStyle="--",MarkerIndices=1:10:length(dbs_rate_simulate),MarkerSize=15,LineWidth=1.5)
    ylabel("Mean GPi F range")
end
ax(3).YAxis(1).Color = colors_plot{1}./250;
ax(3).YAxis(2).Color = colors_plot{2}./250;

%% part 2
%this part can be used to generate data for 2 purposes firstly for
%continuous gradient of increasing dbs amplitude and secondly for a single
%dbs turn on. For first case have repetitions equal to 27500 and comment
%fixed dbs amplitude in lines 301-314 including the added and condition on
%the loop. Else change repetitions to 10000 and uncomment fixed switch
%with added && condition in the if statements lines 301-314.
clear desiredposition muscleposition output
neuron_fire = 0.80; %Define the default neuron excitability. In this case neuron has 85% chance firing given a single default excitatory release

repetitions = 10000; %Define number of timeslices run
%to run data for left hand side of figure change to 10000 else 27500

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
A = 0.15;
desiredposition(2,:) = C + A * sin(2*pi*f*ttime+(pi*0.5));

resting_md = 0.05;
n = length(ttime);
motor_drive = zeros(1, n);
motor_drive(1,:) = 0.4;

idx2 = 1 : n;

muscleposition_str = 0.5;

% Assign values
% motor_drive(2, idx2) = resting_md;

% Define segment indices

% Assign values

hyperdirect_firing = zeros(1,n);
cerebellum_VL_weight = [0.04 0.04 0.04 0.04]; %0.025 for PD
cerebellum_CM_weight = [0.05 0.05 0.05 0.05]; %0.06 for PD
stn_bg_connectivity = [2 2 2 2]; %0.8
bg_cer_connectivity = [1.2 1.2 1.2 1.2];
L5_CM_same_weight = [0.5 0.5 0.5 0.5];
L5_CM_opp_weight = [0.2 0.2 0.2 0.2];
CM_striatum_w = [0.24 0.24 0.24 0.24];
stn_fire_rate = [0.3 0.3 0.3 0.3];
cortex_acth_modifier = [0.75 0.75 0.75 0.75];
BG_CM_weight = [1 1 1 1];
BG_VL_weight = [2.5 2.5 2.5 2.5];
CM_acth_modifier = [1 1 1 1];
ach_levels = [1.5 1.5 1.5 1.5];

CM_diffuse = zeros(size(motor_drive,1),fix(0.200/(1/fs)));
% CM_cummulative = zeros(1,fix(0.200/(1/fs)));

gpi_agonist_keep = zeros(1,repetitions);
gpi_antagonist_keep = zeros(1,repetitions);

acth_rr_agonist = zeros(size(desiredposition,1),repetitions);
acth_rr_antagonist = zeros(size(desiredposition,1),repetitions);

L5_flexor_prev = zeros(size(desiredposition,1),circuits);
L5_extensor_prev = zeros(size(desiredposition,1),circuits);

parkinsonian_level = [0.4 0.4 0.4 0.4];
levodopa_addition = [0 0 0 0];
anticholinergic = [1 1 1 1];
amantadine = [0.5 0.5 0.5 0.5];
stn_gpi_w = [1 1 1 1];

%stn dbs settings
x = 0.1;          % increment size
y = repetitions/11;          % number of steps before increment
N = repetitions;         % total length of array

idx = 0:(N-1);                  % index vector
dbs_stn_inhibition = floor(idx / y) * x;       % compute increments

%GPi dbs settings
x = 0.1;          % increment size
y = repetitions/11;          % number of steps before increment
N = repetitions;         % total length of array

idx = 0:(N-1);                  % index vector
dbs_gpi_inhibition = floor(idx / y) * x;       % compute increments

% %GPi dbs settings
% x = 0.1;          % increment size
% y = repetitions/11;          % number of steps before increment
% N = repetitions;         % total length of array
% 
% idx = 0:(N-1);                  % index vector
% dbs_thalamus_inhibition = floor(idx / y) * x;       % compute increments

for condition_loop = 1:4
    muscleposition = ones(1,repetitions)*muscleposition_str;
    weight_params.keepdefaults = true;
    pathologies.keepdefaults = true;
    weight_params.cerebellum_VL_weight = cerebellum_VL_weight(condition_loop);
    weight_params.cerebellum_CM_weight = cerebellum_CM_weight(condition_loop);
    weight_params.CM_acth_modifier = CM_acth_modifier(condition_loop);
    weight_params.CM_striatum_w = CM_striatum_w(condition_loop);
    weight_params.L5_CM_same_weight = L5_CM_same_weight(condition_loop);
    weight_params.L5_CM_opp_weight = L5_CM_opp_weight(condition_loop);
    weight_params.stn_fire_rate = stn_fire_rate(condition_loop);
    weight_params.cortex_acth_modifier = cortex_acth_modifier(condition_loop);
    weight_params.BG_CM_weight = BG_CM_weight(condition_loop);
    weight_params.BG_VL_weight = BG_VL_weight(condition_loop);
    fig3b = struct;
    num_parts_graph = 48;

    CM_nucleus_template.agonist = zeros(1,fix(0.200/(1/fs)));
    CM_nucleus_template.antagonist = zeros(1,fix(0.200/(1/fs)));

    parfor cor_index = 1:num_parts_graph %parfor

        motor_drive_local = zeros(1,length(amantadine));
        motor_drive_local = motor_drive;

        golgi_cellfire_flexor = zeros(100,1);
        golgi_cellfire_extensor = zeros(100,1);
        cerebellar_agonist = zeros(1,100);
        cerebellar_antagonist = zeros(1,100);
        muscleposition = ones(1,repetitions+1)*0.5;
        output = zeros(1,repetitions);
        output_mp = zeros(1,2);
        CM_index = 1;

        thalamus_flexor_collective = zeros(1,repetitions);
        thalamus_extensor_collective = zeros(1,repetitions);
        thalamus_net = zeros(1,repetitions);
        thalamus_summative = zeros(1,repetitions);
        velocity = zeros(1,repetitions);
        error_observed = zeros(1,repetitions);
        muscleposition_keep = zeros(1,repetitions);

        CM_nucleus = CM_nucleus_template;

        L5_flexor_prev = zeros(1,circuits);
        L5_extensor_prev = zeros(1,circuits);

        if(mod(cor_index,2) == 0)
            dop_level_release = 0.8;
            % elseif(mod(cor_index,3)== 0)
            %     dop_level_release = 0;
        else
            dop_level_release = 0.75;
            motor_drive_local(1,:) = resting_md;
        end

        tic
        for i = 1:length(ttime)-1
            weight_internal = weight_params;
            if(condition_loop == 2&&i>=repetitions/2)
                %weight_internal.dbs_stn_a = dbs_stn_inhibition(1,i);
                weight_internal.dbs_stn_f = 130;
                weight_internal.dbs_stn_a = 0.8;
            elseif(condition_loop == 3&&i>=repetitions/2)
                %weight_internal.dbs_gpi_a = dbs_gpi_inhibition(1,i);
                weight_internal.dbs_gpi_f = 130;
                weight_internal.dbs_gpi_a = 0.9;
            elseif(condition_loop == 4&&i>=repetitions/2)
                %weight_internal.dbs_thalamus_a = dbs_gpi_inhibition(1,i);
                weight_internal.dbs_thalamus_f = 130;
                weight_internal.dbs_thalamus_a = 0.5;
            end
            %% BG
            thalamus_flexor = zeros(2,circuits);
            thalamus_extensor = zeros(2,circuits);
            for motor_programs = 1:1

                %[b_ganglia_dep_probabilities] = basal_ganglia_fcn_pts_rtl(motor_drive(motor_programs,i),2,dopamine_levels(i),weight_params,pathologies,10,hyperdirect_firing(i),CM_diffuse(motor_programs,:),CM_index);
                [b_ganglia_dep_probabilities_agonist] = basal_ganglia_fcn_dbs(motor_drive_local(motor_programs,i),ach_levels(condition_loop),dop_level_release,weight_internal,pathologies,10,hyperdirect_firing(i),CM_nucleus.agonist(motor_programs,:),CM_index);
                [b_ganglia_dep_probabilities_antagonist] = basal_ganglia_fcn_dbs(motor_drive_local(motor_programs,i),ach_levels(condition_loop),dop_level_release,weight_internal,pathologies,10,hyperdirect_firing(i),CM_nucleus.antagonist(motor_programs,:),CM_index);

                %% Cerebellum
                % desiredposition(i) = centre+sin(i/(cycle/(2*pi)))*amplitude;
                desired_input = (0.04 + 0.32*desiredposition(motor_programs,i));
                muscleposition_input = 0.04 + 0.32*muscleposition(i);
                desired_input_rev = (0.04 + 0.32) - desired_input + 0.04; %corrected by adding + 0.04
                muscleposition_input_rev = (0.04 + 0.32) - muscleposition_input + 0.04; %corrected by adding + 0.04

                desired_input = put_in_range(desired_input - bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities_agonist.gpi_dep_probability,0,0.36);
                muscleposition_input = put_in_range(muscleposition_input + stn_bg_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities_antagonist.stn_dep_probability ...
                    -0.075*bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities_antagonist.gpi_dep_probability,0,0.36);
                desired_input_rev = put_in_range(desired_input_rev - bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities_antagonist.gpi_dep_probability,0,0.36);
                muscleposition_input_rev = put_in_range(muscleposition_input_rev + stn_bg_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities_agonist.stn_dep_probability...
                    -0.075*bg_cer_connectivity(condition_loop)*rand*b_ganglia_dep_probabilities_agonist.gpi_dep_probability,0,0.36);

                if(dop_level_release == 0.75)
                    desired_input = 0;
                    desired_input_rev = 0;
                end
                % disp([desired_input,muscleposition_input,desired_input_rev,muscleposition_input_rev])
                for y = 1:100
                    [cerebellar_agonist(y),golgi_cellfire_flexor(y),~] = hbm_vnr(desired_input,muscleposition_input_rev,golgi_cellfire_flexor(y),mod_inhibition,mod_excitation,neuron_fire);

                    [cerebellar_antagonist(y),golgi_cellfire_extensor(y),~] = hbm_vnr(desired_input_rev,muscleposition_input,golgi_cellfire_extensor(y),mod_inhibition,mod_excitation,neuron_fire);
                end

                %% M1
                CM_index = mod(i,length(CM_nucleus.agonist)) + ((mod(i,length(CM_nucleus.agonist)) == 0)*length(CM_nucleus.agonist));
                % [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_diffuse(1,CM_index)] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,b_ganglia_dep_probabilities.gpi_firings,motor_drive(1,i),circuits,weight_params);
            
                if(condition_loop == 4)
                    [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_tmp] = M1_fcn_pts_dbs(cerebellar_agonist,cerebellar_antagonist,...
                    b_ganglia_dep_probabilities_agonist.gpi_firings,b_ganglia_dep_probabilities_antagonist.gpi_firings,...
                    motor_drive_local(1,i),circuits,L5_flexor_prev(motor_programs,:),L5_extensor_prev(motor_programs,:),weight_internal);
                else
                [L5_flexor,L5_extensor,thalamus_flexor(1,:),thalamus_extensor(1,:),CM_tmp] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,...
                    b_ganglia_dep_probabilities_agonist.gpi_firings,b_ganglia_dep_probabilities_antagonist.gpi_firings,...
                    motor_drive_local(1,i),circuits,L5_flexor_prev(motor_programs,:),L5_extensor_prev(motor_programs,:),weight_internal);
                end
                L5_flexor_prev(motor_programs,:) = L5_flexor;
                L5_extensor_prev(motor_programs,:) = L5_extensor;

                CM_nucleus.agonist(motor_programs,CM_index) = CM_tmp.CM_agonist;
                CM_nucleus.antagonist(motor_programs,CM_index) = CM_tmp.CM_antagonist;
                CM_nucleus.greatest(motor_programs,CM_index) = max(CM_nucleus.agonist(motor_programs,CM_index),CM_nucleus.antagonist(motor_programs,CM_index));
                CM_nucleus.mean(motor_programs,CM_index) = (CM_nucleus.agonist(motor_programs,CM_index)+CM_nucleus.antagonist(motor_programs,CM_index))/2;

                %% Spinal Circuits

                spinal_circuit_output = spinal_circuits_fcn_pts(L5_flexor,L5_extensor,muscleposition(i),weight_internal);

                output_mp(motor_programs) = sum(spinal_circuit_output.flexorpower*gain/10,"all")-sum(spinal_circuit_output.extensorpower*gain/10,"all"); %See Eq. 20

            end

            output(i) = sum(output_mp);
            muscleposition_keep(1,i) = muscleposition(i);
            muscleposition(i+1) = muscleposition(i) + output(i)-resistance; %See Eq. 21

            if(muscleposition(i+1)<0)
                muscleposition(i+1) = 0;
            elseif(muscleposition(i+1)>1)
                muscleposition(i+1) = 1;
            end

            thalamus_flexor_collective(1,i) = sum(thalamus_flexor,"all");
            thalamus_extensor_collective(1,i) = sum(thalamus_extensor,"all");
            thalamus_net(1,i) = thalamus_flexor_collective(1,i) - thalamus_extensor_collective(1,i);
            thalamus_summative(1,i) = thalamus_flexor_collective(1,i) + thalamus_extensor_collective(1,i);
            velocity(1,i) = output(i);
            error_observed(1,i) = desiredposition(1,i) - muscleposition(i);

            %Define boundaries of joint range of motion

        end
        [tmp_power,freqs] = cwt(velocity,fs);
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


        correlation_keep(condition_loop,cor_index) = min(corrcoef(muscleposition(1+(repetitions/2):end),desiredposition(1,1+(repetitions/2):end)),[],'all');
        [r,lags] = xcorr(thalamus_net(1,1+(repetitions/2):repetitions),desiredposition(1,1+(repetitions/2):repetitions));
        r(lags>=0) = 0;
        [max_corthalamus(condition_loop,cor_index),max_index_tmp] = max(r,[],"all");
        cor_thalamus_error_centre(condition_loop,cor_index) = r(repetitions/2);
        lag_keep(condition_loop,cor_index) = lags(max_index_tmp);
        corr_thalamus_error(condition_loop,cor_index) = min(corrcoef(thalamus_net(1+(repetitions/2):end),error_observed(1,1+(repetitions/2):end)),[],'all');
        corr_thalamus_desired(condition_loop,cor_index) = min(corrcoef(thalamus_net(1+(repetitions/2):end),desiredposition(1,1+(repetitions/2):end-1)),[],'all');


        fig3b_thalamus_flexor_collective(condition_loop,cor_index,:) = thalamus_flexor_collective;
        fig3b_thalamus_extensor_collective(condition_loop,cor_index,:) = thalamus_extensor_collective;
        fig3b_thalamus_net(condition_loop,cor_index,:) = thalamus_net;
        fig3b_thalamus_summative(condition_loop,cor_index,:) = thalamus_summative;
        fig3b_velocity(condition_loop,cor_index,:) = velocity;
        fig3b_error_observed(condition_loop,cor_index,:) = error_observed;
        fig3b_muscleposition(condition_loop,cor_index,:) = muscleposition_keep;

        toc
    end
end

% tmp(:,:) = fig3b_muscleposition(4,:,:);
% plot(tmp.')
% hold on
% plot(desiredposition(1,:))


% figure
% clear tmp
% ttime = 0:1/250:(repetitions-1)/250;
% for index_plotbelow = 1:3
%     tmp(:,:) = fig3b_muscleposition(index_plotbelow+1,:,:);
%     subplot(4,3,4+(index_plotbelow-1)*3)
%     plot(ttime,tmp(1,:))
%     hold on
%     % plot(ttime,desiredposition(1,1:end-1))
%     xlim([0 repetitions/250])
%     % ylim([0 1])
%     subplot(4,3,5+(index_plotbelow-1)*3)
%     plot(ttime,tmp(2,:))
%     hold on
%     plot(ttime,desiredposition(1,1:end-1))
%     xlim([0 repetitions/250])
%     ylim([0 1])
% end

%% ploting
clear all
load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\fig_dbs_top_3march.mat")
num_circuits_bg = 1000;
motor_drive = 0.01:0.005:0.4;%0.01:0.005:0.4;
motor_drive_frequency = motor_drive*250;
% dopamine_levels = 0:0.1:3;
dopamine_levels = 0:0.025:3;
dop_acth = dopamine_levels/2;
titles_part1 = ["Healthy" "PD" "Levodopa" "Anticholinergic" "Amantadine"];
figure
plot_sequence = [3 3 2 4 4];
index_plot = [10 25 1 16 25];
counter1 = 1;
counter2 = 1;
dbs_rate_simulate = 0:0.025:1;

repetitions = 27500; %Define number of timeslices run
A = 0.2;          % Amplitude of sine wave
T = 10;        % Length of one cycle (period) in seconds
C = 0.5;          % Centre (vertical offset)
fs = 250;      % Sampling frequency (points per second)
L = repetitions/fs;          % Total signal length in seconds

% Derived values
f = 1/T;                    % Frequency (Hz)
ttime = 0:1/fs:L;               % Time vector

% Generate sine wave with offset
desiredposition(1,:) = C + A * sin(2*pi*f*ttime);%ones(1,repetitions+1)*0.5;


for index_dif = 1:2
    condition_checked(4+index_dif).gpi_dp = movmean(condition_checked(2).gpi_dp(1,:,:),50,1) - movmean(condition_checked(2+index_dif).gpi_dp(:,:,:),50,1);
    mean_reduction(index_dif,:) = mean(condition_checked(4+index_dif).gpi_dp,[2 3]).';
    range_reduction(index_dif,:) = 1.96*std(condition_checked(4+index_dif).gpi_dp,0,[2 3]).';
    min_activation(index_dif,:) = min(condition_checked(2+index_dif).gpi_dp,[],[2 3]);
    max_activation(index_dif,:) = max(condition_checked(2+index_dif).gpi_dp,[],[2 3]);
    mean_activation(index_dif,:) = (min_activation(index_dif,:) + max_activation(index_dif,:))./2;
end

for i = 1:length(plot_sequence)
    if(plot_sequence(i)~=0)
        ax(i) = subplot(4,length(plot_sequence),i);
        toplot(:,:) = condition_checked(plot_sequence(i)).gpi_dp(index_plot(i),:,:);
        % if(plot_sequence(i)<=4)
        % if(i<=2)
            % mxfreq(counter1) = max(condition_checked(plot_sequence(i)).gpi_dp,[],"all");
            % minfreq(counter1) = min(condition_checked(plot_sequence(i)).gpi_dp,[],"all");
            % mxfreq(counter1) = max(toplot,[],"all");
            % minfreq(counter1) = min(toplot,[],"all");
            counter1 = counter1 + 1;
            % imagesc(dop_acth,motor_drive_frequency,movmean(toplot.',10))
            % 
            % hold on
            contourf(dop_acth,motor_drive_frequency,movmean(toplot.',25),20,"LineWidth",0.5,"LineColor","k")
            colormap(jet(100));

        % else
        %     mxfreq2(counter2) = max(toplot,[],"all");
        %     minfreq2(counter2) = min(toplot,[],"all");
        %     counter2  = counter2 + 1;
        %     contourf(dop_acth,motor_drive_frequency,movmean(toplot.',20),20)
        % end
        hold on;
        xlabel('Dop - ACh');
        ylabel('Motor Drive Frequency (Hz)');
    end
end

for i = 1:5
    clim(ax(i),[20,140])
    % clim(ax(2+(i-1)*3),[min(minfreq2),max(mxfreq2)])
end

% ax(3) = subplot(4,length(plot_sequence),3);
% colors_plot = {[214,124,48],[62,150,81]};
% markers_use = ['*' 'x'];
% for i = 1:size(mean_reduction,1)
%     % X=[dbs_rate_simulate,fliplr(dbs_rate_simulate)];                %#create continuous x value array for plotting
%     % Y=[min_activation(i,:),fliplr(max_activation(i,:))];              %#create y values for out and then back
%     % fill(X,Y,'b',"FaceAlpha",0.2);
%     hold on
%     yyaxis left
%     plot(dbs_rate_simulate,movmean(mean_reduction(i,:),10),"Color",colors_plot{1}./250,Marker=markers_use(i)...
%         ,LineStyle="--",MarkerIndices=1:10:length(dbs_rate_simulate),MarkerSize=10,LineWidth=1.5)
%     ylabel("Mean GPi F reduction")
%     yyaxis right
%     plot(dbs_rate_simulate,movmean(max_activation(i,:)-min_activation(i,:),10),"Color",colors_plot{2}./250,Marker=markers_use(i)...
%         ,LineStyle="--",MarkerIndices=1:10:length(dbs_rate_simulate),MarkerSize=10,LineWidth=1.5)
%     ylabel("Mean GPi F range")
% end
% ax(3).YAxis(1).Color = colors_plot{1}./250;
% ax(3).YAxis(2).Color = colors_plot{2}./250;


load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\fig_dbs_bottom_left_3march.mat")
clear tmp
ttime = 0:1/250:(10000-1)/250;
for index_plotbelow = 1:3
    tmp(:,:) = fig3b_muscleposition(index_plotbelow+1,:,:);
    subplot(4,3,4+(index_plotbelow-1)*3)
    if(index_plotbelow==1)
        plot(ttime,tmp(13,:),"LineWidth",2)
    elseif(index_plotbelow==2)
        plot(ttime,tmp(1,:),"LineWidth",2)
    else
        plot(ttime,tmp(9,:),"LineWidth",2)
    end
    hold on
    plot([20 20],[0 1],Color="k",LineStyle="--",LineWidth=1)
    if(index_plotbelow == 1)
        title("Resting state")
    elseif(index_plotbelow == 3)
        xlabel("Time(s)")
    elseif(index_plotbelow == 2)
        ylabel("Joint Position")
    end
    % plot(ttime,(movvar(tmp(1,:),1000)*10^3)+0.4)
    % plot(ttime,desiredposition(1,1:end-1))
    xlim([0 40])
    ylim([0.1 0.9])
    subplot(4,3,5+(index_plotbelow-1)*3)
    if(index_plotbelow==1)
        plot(ttime,tmp(6,:),"LineWidth",2)
    elseif(index_plotbelow==2)
        plot(ttime,tmp(22,:),"LineWidth",2) %22 is good no dyskinesia 34 is good 36 like 34 more diskinesia
    else
        plot(ttime,tmp(8,:),"LineWidth",2)
    end
    hold on
    plot(ttime,desiredposition(1,1:10000),"LineWidth",2)
    plot([20 20],[0 1],Color="k",LineStyle="--",LineWidth=1)
    if(index_plotbelow == 1)
        title("Position seeking state")
    end
    xlim([0 40])
    ylim([0.1 0.9])
end
markers_use = {'+','*','x'};
load("C:\Users\grego\OneDrive - Queen Mary, University of London\manuscripts\accepted\cerebellum_paper\code\basal_ganglia\fig_dbs_bottom_right_3march.mat")
clear tmp
dbs_inhibition = 0:0.1:1;
subplot(4,3,6)
cmap = lines(3);
hold on
for i = 1:size(fig3b_muscleposition,1)
    counterpd = 1;
    counterspd = 1;
    for y = 1:size(fig3b_muscleposition,2)
        tmp(:,:) = fig3b_error_observed(i,y,:);
        error_tmp = abs(movmean(tmp,1));
        tremor_tmp = movvar(tmp,200);
        if(mod(y,2) == 0)
            for z = 1:11
                error_pd(counterpd,z) = mean(abs(error_tmp(1+(z-1)*2500:z*2500)),"all");
            end
            counterpd = counterpd + 1;
        else
            for z = 1:11
                error_spd(counterspd,z) = mean(abs(error_tmp(1+(z-1)*2500:z*2500)),"all");
            end
            counterspd = counterspd + 1;
        end
    end
    if(i == 1)
        error_control_pd = error_pd;
        error_control_spd = error_spd;
    else

        if(i~=4)
            yyaxis left
            plot(dbs_inhibition*100,movmean((mean(error_pd,1)-mean(error_control_pd,1))./mean(error_control_pd,1)*100,3),Color=cmap(2-1,:),Marker=markers_use{i-1},LineWidth=1.5,LineStyle="--",MarkerSize=10)

        else
            yyaxis right
            plot(dbs_inhibition*100,movmean((mean(error_pd,1)-mean(error_control_pd,1))./mean(error_control_pd,1)*100,3),Color=cmap(3-1,:),Marker=markers_use{i-1},LineWidth=1.5,LineStyle="--",MarkerSize=10)

        end

        %plot(dbs_inhibition,(mean(velolicity_pd,1)-mean(velocity_control_pd,1))./mean(velocity_control_pd,1),Color=cmap(i-1,:),Marker="*")
        ylabel("% Inhibition")
        title("Error % change")
        % yyaxis right
        % plot(dbs_inhibition,mean(tremor_pd,1),Color=cmap(i))
        % plot(dbs_inhibition,mean(tremor_spd,1),Color=cmap(i))
    end
end

clear velolicity_pd velolicity_spd velocity_control_pd velocity_control_spd
dbs_inhibition = 0:0.1:1;
subplot(4,3,12)
cmap = lines(3);
markers_use = {'+','*','x'};
hold on
for i = 1:size(fig3b_muscleposition,1)
    counterpd = 1;
    counterspd = 1;
    counterpd_nodbs = 1;
    for y = 1:size(fig3b_muscleposition,2)
        tmp(:,:) = fig3b_muscleposition(i,y,:);
        movement_tmp = movmean(tmp,2000);
        tremor_tmp = movstd(tmp,200);
        if(mod(y,2) == 0)
            for z = 1:11
                velolicity_pd(counterpd,z) = mean(movement_tmp(1+(z-1)*2500:z*2500),"all");
                tremor_pd(counterpd,z) = mean(tremor_tmp(1+(z-1)*2500:z*2500),"all");
                tremor_velocity_pd(counterpd,z) = tremor_pd(counterpd,z)/velolicity_pd(counterpd,z);
            end
            counterpd = counterpd + 1;
        elseif(mod(y,1) == 0)
            for z = 1:11
                velolicity_spd(counterspd,z) = mean(movement_tmp(1+(z-1)*2500:z*2500),"all");
                tremor_spd(counterspd,z) = mean(tremor_tmp(1+(z-1)*2500:z*2500),"all");
                tremor_velocity_spd(counterspd,z) = tremor_spd(counterspd,z)/velolicity_spd(counterspd,z);
            end
            counterspd = counterspd + 1;
        end
    end
    if(i == 1)
        velocity_control_pd = velolicity_pd;
        velocity_control_spd = velolicity_spd;
        tremor_control_pd = tremor_pd;
        tremor_control_spd = tremor_spd;
        tremor_velocity_control_pd = tremor_velocity_pd;
        tremor_velocity_control_spd = tremor_velocity_spd;
    else
        % if(i~=3)
        %     % yyaxis left
        %     plot(dbs_inhibition*100,movmean((mean(tremor_spd,1)-mean(tremor_control_spd,1))./mean(tremor_control_spd,1)*100,3),Color=cmap(2-1,:),Marker=markers_use{i-1},LineWidth=1.5,LineStyle="--")
        % 
        % else
            % yyaxis right
            plot(dbs_inhibition*100,movmean((mean(tremor_spd,1)-mean(tremor_control_spd,1))./mean(tremor_control_spd,1)*100,3),Color='k',Marker=markers_use{i-1},LineWidth=1.5,LineStyle="--",MarkerSize=10)

        % end
        %plot(dbs_inhibition,(mean(tremor_velocity_pd,1)-mean(tremor_velocity_control_pd,1))./mean(tremor_velocity_control_pd,1),Color=cmap(i-1,:),Marker="*")
        xlabel("% Inhibition")
        title("Resting tremor % change")
        % yyaxis right
        % plot(dbs_inhibition,mean(tremor_pd,1),Color=cmap(i))
        % plot(dbs_inhibition,mean(tremor_spd,1),Color=cmap(i))
    end
end

% % dbs_inhibition = 0:0.1:1;
% % subplot(4,3,12)
% % cmap = lines(3);
% % hold on
% % for i = 1:size(fig3b_muscleposition,1)
% %     counterpd = 1;
% %     counterspd = 1;
% %     counterpd_nodbs = 1;
% %     for y = 1:size(fig3b_muscleposition,2)
% %         tmp(:,:) = fig3b_muscleposition(i,y,:);
% %         movement_tmp = movmean(tmp,200);
% %         tremor_tmp = movvar(tmp,200);
% %         if(~(mod(y,2) == 0)&&~(mod(y,3) == 0))
% %             for z = 1:11
% %                velolicity_pd(counterpd,z) = mean(movement_tmp(1+(z-1)*2500:z*2500),"all");
% %                tremor_pd(counterpd,z) = std(mean(tremor_tmp(1+(z-1)*2500:z*2500),"all"));
% %                tremor_velocity_pd(counterpd,z) = tremor_pd(counterpd,z)/velolicity_pd(counterpd,z);
% %             end
% %             counterpd = counterpd + 1;
% %         end
% %     end
% %     if(i == 1)
% %         velocity_control_pd = velolicity_pd;
% %         tremor_control_pd = tremor_pd;
% %         tremor_velocity_control_pd = tremor_velocity_pd;
% %     else
% %
% %         % yyaxis left
% %         plot(dbs_inhibition,(mean(velolicity_pd,1)-mean(velocity_control_pd,1))./mean(velocity_control_pd,1),Color=cmap(i-1,:),Marker="*")
% %         xlabel("DBS inhibition")
% %         ylabel("Change Rest movement(%)")
% %         % plot(dbs_inhibition,(mean(tremor_velocity_spd,1)-mean(tremor_velocity_control_spd,1))./mean(tremor_velocity_control_spd,1),Color=cmap(i-1,:),Marker="x")
% %         % yyaxis right
% %         % plot(dbs_inhibition,mean(tremor_pd,1),Color=cmap(i))
% %         % plot(dbs_inhibition,mean(tremor_spd,1),Color=cmap(i))
% %     end
% % end


%% find percentage differce

