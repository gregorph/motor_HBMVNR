function [L5_flexor,L5_extensor,VL_flexor,VL_extensor,CM_nucleus] = M1_fcn_pts_rtl(cerebellar_agonist,cerebellar_antagonist,bg_gpi_fire_agonist,bg_gpi_fire_antagonist,motor_drive,circuits,L5_flexor_prev,L5_extensor_prev,weight_params)
%M1_FCN Summary of this function goes here
%   Detailed explanation goes here
%%  Settings
    if(~exist("circuits"))
        circuits = 1000;
    end
    
        neuron_fire = 0.8;
    
        %Neurotransmitter settings
        gaba_mod = -0.5;
        glutamine_mod = 0.5;
    
        %Synaptic weights
        cerebellum_VL_weight = 0.0135;%0.0135 physiology 0.1 resting tremor PD
        BG_VL_weight = 3;
        cerebellum_CM_weight = 0.005;%0.005 physiology 
        BG_CM_weight = 1.2;
        VL_L4_weight = 0.6;
        L4_PV_weight = 0.6;
        motor_L23_weight = 0.3;
        L4_L23_weight = 0.24;
        PV_L23_weight = 0.4;
        L23_L5_weight = 0.16;
        L5_CM_same_weight = 0.1;
        L5_CM_opp_weight = 0.04;
        thalamus_dbs = 0;
        L4_dbs = 0;
        L23_dbs = 0;
        L5_dbs = 0;
        
    
    var_list = { ...
    'neuron_fire','gaba_mod','glutamine_mod', ...
    'cerebellum_VL_weight','BG_VL_weight','VL_L4_weight','L4_PV_weight','motor_L23_weight','L4_L23_weight', ...
    'PV_L23_weight','L23_L5_weight','thalamus_dbs','L4_dbs','L23_dbs','L5_dbs'...
     'cerebellum_CM_weight','BG_CM_weight','L5_CM_same_weight','L5_CM_opp_weight'};

for i = 1:numel(var_list)
    v = var_list{i};
    if isfield(weight_params, v)
        eval([v ' = weight_params.(v);']);
    end
end
    
thalamus_dbs_i = abs(thalamus_dbs) * (thalamus_dbs>=0);
thalamus_dbs_e = abs(thalamus_dbs) * (thalamus_dbs<=0);
L4_dbs_i = abs(L4_dbs) * (L4_dbs>=0);
L4_dbs_e = abs(L4_dbs) * (L4_dbs<=0);
L23_dbs_i = abs(L23_dbs) * (L23_dbs>=0);
L23_dbs_e = abs(L23_dbs) * (L23_dbs<=0);
L5_dbs_i = abs(L5_dbs) * (L5_dbs>=0);
L5_dbs_e = abs(L5_dbs) * (L5_dbs<=0);

        L4_flexor = zeros(1,10);
        L4_extensor = zeros(1,10);

        CM_diffuse_agonist_tmp = zeros(1,10);
        CM_diffuse_antagonist_tmp = zeros(1,10);

        VL_flexor = zeros(1,circuits);
        VL_extensor = zeros(1,circuits);
        CM_diffuse_agonist = zeros(1,circuits);
        CM_diffuse_antagonist = zeros(1,circuits);

        L5_CM_same_weight = L5_CM_same_weight/circuits;
        L5_CM_opp_weight = L5_CM_opp_weight/circuits;


    %% M1 neural network
    for i = 1:circuits
        VL_flexor(i) = (((cerebellum_VL_weight*rand*sum(cerebellar_agonist)*glutamine_mod) +...
            (BG_VL_weight*rand*sum(bg_gpi_fire_agonist)*gaba_mod)) > (1-neuron_fire))&&(rand>thalamus_dbs_i);
        VL_flexor(i) = VL_flexor(i) ||(rand>(1-thalamus_dbs_e));
        VL_extensor(i) = (((cerebellum_VL_weight*rand*sum(cerebellar_antagonist*glutamine_mod)) +...
            (BG_VL_weight*rand*sum(bg_gpi_fire_antagonist) *gaba_mod)) > (1-neuron_fire))&&(rand>thalamus_dbs_i);
        VL_extensor(i) = VL_extensor(i) ||(rand>(1-thalamus_dbs_e));
        
        for k = 1:10
            % CM_diffuse_single_tmp(k) = (((cerebellum_CM_weight*rand*(sum(cerebellar_agonist)+sum(cerebellar_antagonist))*glutamine_mod) + ...
            %     (BG_CM_weight*rand*sum(bg_gpi_fire)*gaba_mod))...
            %     > (1-neuron_fire))&&(rand>thalamus_dbs_i);
            % CM_diffuse_single_tmp(k) = CM_diffuse_single_tmp(k) ||(rand>(1-thalamus_dbs_e));

            CM_diffuse_agonist_tmp(k) = (((cerebellum_CM_weight*rand*sum(cerebellar_agonist)*glutamine_mod) + ...
                (BG_CM_weight*rand*sum(bg_gpi_fire_agonist)*gaba_mod))...
                + ((rand*sum(L5_flexor_prev)*L5_CM_same_weight*glutamine_mod)+(rand*sum(L5_extensor_prev)*L5_CM_opp_weight*glutamine_mod))...
                > (1-neuron_fire))&&(rand>thalamus_dbs_i);
            CM_diffuse_agonist_tmp(k) = CM_diffuse_agonist_tmp(k) ||(rand>(1-thalamus_dbs_e));

            CM_diffuse_antagonist_tmp(k) = (((cerebellum_CM_weight*rand*sum(cerebellar_antagonist)*glutamine_mod) + ...
                (BG_CM_weight*rand*sum(bg_gpi_fire_antagonist)*gaba_mod))...
                + ((rand*sum(L5_extensor_prev)*L5_CM_same_weight*glutamine_mod)+(rand*sum(L5_flexor_prev)*L5_CM_opp_weight*glutamine_mod))...
                > (1-neuron_fire))&&(rand>thalamus_dbs_i);
            CM_diffuse_antagonist_tmp(k) = CM_diffuse_antagonist_tmp(k) ||(rand>(1-thalamus_dbs_e));
        end

        % CM_diffuse_single(i) = mean(CM_diffuse_single_tmp);
        CM_diffuse_agonist(i) = mean(CM_diffuse_agonist_tmp);
        CM_diffuse_antagonist(i) = mean(CM_diffuse_antagonist_tmp);
        
        
        for g = 1:10
            L4_flexor(1,g) = ((VL_L4_weight*rand*VL_flexor(i)*glutamine_mod) > (1-neuron_fire))&&(rand>L4_dbs_i);
            L4_flexor(1,g) = L4_flexor(1,g) ||(rand>(1-L4_dbs_e));
            L4_extensor(1,g) = ((VL_L4_weight*rand*VL_extensor(i)*glutamine_mod) > (1-neuron_fire))&&(rand>L4_dbs_i);
            L4_extensor(1,g) = L4_extensor(1,g) ||(rand>(1-L4_dbs_e));
        end
        
        L23_flexor = zeros(1,4);
        L23_extensor = zeros(1,4);
        for g = 1:4
        
            motor_drive_fire = rand >= (1-motor_drive);
    
            PV_flexor = (L4_PV_weight*rand*sum(L4_flexor,"all")*glutamine_mod) > (1-neuron_fire);
            L23_flexor(g) = (((rand*motor_L23_weight*motor_drive_fire*glutamine_mod)+(rand*L4_L23_weight*sum(L4_flexor)*glutamine_mod)+(rand*PV_flexor*PV_L23_weight*gaba_mod)) > (1-neuron_fire))&&(rand>L23_dbs_i);
            L23_flexor(g) = L23_flexor(g) ||(rand>(1-L23_dbs_e));

            PV_extensor = (L4_PV_weight*rand*sum(L4_extensor,"all")*glutamine_mod) > (1-neuron_fire);
            L23_extensor(g) = (((rand*motor_L23_weight*motor_drive_fire*glutamine_mod)+(rand*L4_L23_weight*sum(L4_extensor)*glutamine_mod)+(rand*PV_extensor*PV_L23_weight*gaba_mod)) > (1-neuron_fire))&&(rand>L23_dbs_i);
            L23_extensor(g) = L23_extensor(g) ||(rand>(1-L23_dbs_e));
        end
        
        L5_flexor(i) = ((rand*L23_L5_weight*sum(L23_flexor,"all")*glutamine_mod) > (1-neuron_fire))&&(rand>L5_dbs_i);
        L5_flexor(i) = L5_flexor(i) ||(rand>(1-L5_dbs_e));
        
        L5_extensor(i) = ((rand*L23_L5_weight*sum(L23_extensor,"all")*glutamine_mod) > (1-neuron_fire))&&(rand>L5_dbs_i);
        L5_extensor(i) = L5_extensor(i) ||(rand>(1-L5_dbs_e));
    end
    CM_nucleus.CM_agonist = mean(CM_diffuse_agonist);
    CM_nucleus.CM_antagonist = mean(CM_diffuse_antagonist);
end

