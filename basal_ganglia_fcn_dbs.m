function [b_ganglia_dep_probabilities] = basal_ganglia_fcn_dbs(cortex_dep_probability,acth_release_rate,dopamine_release_rate,weight_params,pathologies,circuits,hyperdirect_firing,CM_diffuse,current_index)
%BASAL_GANGLIA_FCN Summary of this function goes here
%   Detailed explanation goes here
%%  Settings
if(~exist("circuits"))
    circuits = 10;
end
if(pathologies.keepdefaults)
    huntigton_d = 0;
    parkinsons_dop_reduction = 1;
end

neuron_firing = 0.8;

%Neurotransmitter settings
gaba_mod = -0.5;
glutamine_mod = 0.5;

%Synaptic weights
cortex_striatum_w = 1;
CM_striatum_w = 0.25; %5 at last 0
d2_gpe_w = 1.5;
gpe_stn_w = 1;
stn_gpi_w = 0.9;
d1_gpi_w = 1.5;
hyperdirect_stn_weight = 10;

%Intrinsic firing rates
gpe_fire_rate = (1-neuron_firing)/(1-(60/250));
stn_fire_rate = (1-neuron_firing)/(1-(20/250));
gpi_fire_rate = (1-neuron_firing)/(1-(60/250));

dbs_striatum_f = 0;
dbs_stn_f = 0;
dbs_gpe_f = 0;
dbs_gpi_f = 0;

dbs_striatum_a = 0;
dbs_stn_a = 0;
dbs_gpe_a = 0;
dbs_gpi_a = 0;

%Neuromodulators
acth_modifier_d1 = 0.4;
acth_modifier_d2 = 0.5;
dopamine_modifier_d1 = 1 * parkinsons_dop_reduction;
dopamine_modifier_d2 = 0.8 * parkinsons_dop_reduction;
acth_dopamine_modifier = 0.01;
CM_acth_modifier = 0.5;
cortex_acth_modifier = 0.3;%0.4;
CM_fire = zeros(1,10);


var_list = { ...
    'neuron_firing','gaba_mod','glutamine_mod', ...
    'cortex_striatum_w','d2_gpe_w','gpe_stn_w','stn_gpi_w','d1_gpi_w','hyperdirect_stn_weight', ...
    'gpe_fire_rate','stn_fire_rate','gpi_fire_rate','thalamus_fire_rate', ...
    'acth_modifier_d1','acth_modifier_d2','dopamine_modifier_d1','dopamine_modifier_d2','acth_dopamine_modifier', ...
    'dbs_striatum_f','dbs_stn_f','dbs_gpe_f','dbs_gpi_f',...
    'dbs_striatum_a','dbs_stn_a','dbs_gpe_a','dbs_gpi_a',...
    'CM_acth_modifier','cortex_acth_modifier','CM_striatum_w'};

for i = 1:numel(var_list)
    v = var_list{i};
    if isfield(weight_params, v)
        eval([v ' = weight_params.(v);']);
    end
end

CM_acth_modifier = put_in_range(CM_acth_modifier*acth_release_rate*(mean(CM_diffuse,"all")/0.4),0,(0.5*acth_release_rate)); %0.3
cortex_acth_modifier = cortex_acth_modifier*acth_release_rate*(cortex_dep_probability/0.4);


% acth_release_rate_keep = zeros(1,circuits);

acth_release_rate = acth_release_rate - CM_acth_modifier...
    + cortex_acth_modifier;

% acth_release_rate_keep(i) = acth_release_rate;

dopamine_released = dopamine_release_rate*(1+(acth_release_rate*acth_dopamine_modifier));

d2_dep_threshold = (0.5) - (acth_release_rate*acth_modifier_d2) + (dopamine_released*dopamine_modifier_d2);
if(d2_dep_threshold < 0)
    d2_dep_threshold = 0;
end

d1_dep_threshold = (0.5) + (acth_release_rate*acth_modifier_d1) - (dopamine_released*dopamine_modifier_d1);
if(d1_dep_threshold < 0)
    d1_dep_threshold = 0;
end

for i = 1:circuits

    for g = 1:10
        CM_fire(g) = rand > (1-CM_diffuse(current_index));
    end

    % cortex_fire = rand > (1-cortex_dep_probability);
    % cortex_striatum_effect = rand * glutamine_mod * cortex_fire*cortex_striatum_w;

    for g = 1:4
        cortex_fire = rand > (1-cortex_dep_probability);
        cortex_striatum_effect = rand * glutamine_mod * cortex_fire * cortex_striatum_w;
        % CM_fire = rand > (1-CM_diffuse(current_index));
        %CM_striatum_effect = rand * glutamine_mod * CM_fire * CM_striatum_w;
        CM_striatum_effect = rand * glutamine_mod * sum(CM_fire) * CM_striatum_w;
        dbs_striatum_effect = 1+((rand + direction_dbs_logistic_curve(dbs_striatum_f,0) - 0.5)*dbs_striatum_a*2.5);
        MSN_D1(i,g) = (((cortex_striatum_effect + CM_striatum_effect) * dbs_striatum_effect > d1_dep_threshold));
        %MSN_D1(i,g) = (((cortex_striatum_effect) > d1_dep_threshold)&&(rand>=striatum_dbs_i));
        MSN_D1(i,g) = MSN_D1(i,g)||((rand*dbs_striatum_a*0.075)>=(direction_dbs_logistic_curve(dbs_striatum_f,1)));
    end

    for g = 1:6
        cortex_fire = rand > (1-cortex_dep_probability);
        cortex_striatum_effect = rand * glutamine_mod * cortex_fire*cortex_striatum_w;
        % CM_fire = rand > (1-CM_diffuse(current_index));
        %CM_striatum_effect = rand * glutamine_mod * CM_fire * CM_striatum_w;
        dbs_striatum_effect = 1+((rand + direction_dbs_logistic_curve(dbs_striatum_f,0) - 0.5)*dbs_striatum_a*2.5);
        CM_striatum_effect = rand * glutamine_mod * sum(CM_fire) * CM_striatum_w;
        MSN_D2(i,g) = ((cortex_striatum_effect + CM_striatum_effect) * dbs_striatum_effect > d2_dep_threshold);
        %MSN_D2(i,g) = ((cortex_striatum_effect) > d2_dep_threshold)&&(rand>=striatum_dbs_i);
        MSN_D2(i,g) = MSN_D2(i,g)||((rand*dbs_striatum_a*0.08)>=(direction_dbs_logistic_curve(dbs_striatum_f,1)));
    end


    for g = 1:3
        d2_gpe_effect = rand*gaba_mod*sum(MSN_D2(i,:))*d2_gpe_w;
        gpe_gpe_effect = ((rand*gpe_fire_rate)>(1-neuron_firing)) * gaba_mod;
        stn_gpe_effect = ((rand*stn_fire_rate)>(1-neuron_firing)) * glutamine_mod;
        dbs_gpe_effect = 1+((rand + direction_dbs_logistic_curve(dbs_gpe_f,0) - 0.5)*dbs_gpe_a*2.5);
        gpe(i,g) = ((rand*gpe_fire_rate) +  d2_gpe_effect + gpe_gpe_effect + stn_gpe_effect) * dbs_gpe_effect > (1-neuron_firing);
        gpe(i,g) = gpe(i,g) || ((rand*dbs_gpe_a*0.08)>=(direction_dbs_logistic_curve(dbs_gpe_f,1)));
    end

    if(exist("hyperdirect_firing"))
        for g = 1:6
            hyperdirect_fire(g) = rand > (1-hyperdirect_firing);
        end
        hyperdirect_stn_effect = rand*glutamine_mod*sum(hyperdirect_fire)*hyperdirect_stn_weight;
    else
        hyperdirect_stn_effect = 0;
    end
    for g = 1:2
        dbs_stn_effect = 1+((rand + direction_dbs_logistic_curve(dbs_stn_f,0) - 0.5)*dbs_stn_a*2.5);
        gpe_stn_effect = rand*gaba_mod*sum(gpe(i,:))*gpe_stn_w;
        stn(i,g) = ((rand*stn_fire_rate) + gpe_stn_effect + hyperdirect_stn_effect) * dbs_stn_effect> (1-neuron_firing);
        stn(i,g) = stn(i,g) || ((rand*dbs_stn_a*0.08)>=(direction_dbs_logistic_curve(dbs_stn_f,1)));
    end
    d1_gpi_effect = rand*gaba_mod*sum(MSN_D1(i,:))*d1_gpi_w;
    stn_gpi_effect = rand*glutamine_mod*sum(stn(i,:))*stn_gpi_w;
    dbs_gpi_effect = 1+((rand + direction_dbs_logistic_curve(dbs_gpi_f,0) - 0.5)*dbs_gpi_a*2.5);
    gpi(i) = ((rand*gpi_fire_rate) + d1_gpi_effect + stn_gpi_effect)*dbs_gpi_effect > (1-neuron_firing);
    gpi(i) = gpi(i) || ((rand*dbs_gpi_a*0.08)>=(direction_dbs_logistic_curve(dbs_gpi_f,1)));

    % gpi_thalamus_effect = gpi(i)*rand*gaba_mod;
    % thalamus(i) = (rand*thalamus_fire_rate) + gpi_thalamus_effect > (1-neuron_firing);

end
b_ganglia_dep_probabilities.gpi_dep_probability = mean(gpi);

b_ganglia_dep_probabilities.MSN_D1_dep_probability = mean(MSN_D1,"all");

b_ganglia_dep_probabilities.MSN_D2_dep_probability = mean(MSN_D2,"all");

% b_ganglia_dep_probabilities.thalamus_dep_probability = mean(thalamus);

b_ganglia_dep_probabilities.gpe_dep_probability = mean(gpe,"all");

b_ganglia_dep_probabilities.stn_dep_probability = mean(stn,"all");

b_ganglia_dep_probabilities.gpi_firings = gpi;

b_ganglia_dep_probabilities.acth_release_rate = acth_release_rate;%mean(acth_release_rate_keep,"all");

end

