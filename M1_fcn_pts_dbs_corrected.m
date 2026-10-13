function [L5_flexor,L5_extensor,VL_flexor,VL_extensor,CM_nucleus] = M1_fcn_pts_dbs_corrected(cerebellar_agonist,cerebellar_antagonist,bg_gpi_fire_agonist,bg_gpi_fire_antagonist,motor_drive,circuits,L5_flexor_prev,L5_extensor_prev,weight_params)
%M1_FCN_PTS_DBS_CORRECTED Supplied M1 function with CM DBS parentheses fixed.
% DBS now multiplies the entire CM input sum, as it does for VL.
% Cortical L4/L23/L5 DBS equations otherwise retain the supplied formulation.
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
        
        dbs_thalamus_f = 0;
        dbs_L23_f = 0;
        dbs_L5_f = 0;
        dbs_L4_f = 0;

        dbs_thalamus_a = 0;
        dbs_L23_a = 0;
        dbs_L5_a = 0;
        dbs_L4_a = 0;
        
    
    var_list = { ...
    'neuron_fire','gaba_mod','glutamine_mod', ...
    'cerebellum_VL_weight','BG_VL_weight','VL_L4_weight','L4_PV_weight','motor_L23_weight','L4_L23_weight', ...
    'PV_L23_weight','L23_L5_weight','dbs_thalamus_f','dbs_L23_f','dbs_L5_f','dbs_L4_f',...
    'dbs_thalamus_a','dbs_L23_a','dbs_L5_a','dbs_L4_a',...
     'cerebellum_CM_weight','BG_CM_weight','L5_CM_same_weight','L5_CM_opp_weight'};

for i = 1:numel(var_list)
    v = var_list{i};
    if isfield(weight_params, v)
        eval([v ' = weight_params.(v);']);
    end
end
    

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
        dbs_thalamus_effect = 1+((rand + direction_dbs_logistic_curve(dbs_thalamus_f,0) - 0.5)*dbs_thalamus_a*2.5);
        VL_flexor(i) = (((cerebellum_VL_weight*rand*sum(cerebellar_agonist)*glutamine_mod) +...
            (BG_VL_weight*rand*sum(bg_gpi_fire_agonist)*gaba_mod)) * dbs_thalamus_effect > (1-neuron_fire));
        VL_flexor(i) = VL_flexor(i) || ((rand*dbs_thalamus_a*0.08)>=(direction_dbs_logistic_curve(dbs_thalamus_f,1)));

        dbs_thalamus_effect = 1+((rand + direction_dbs_logistic_curve(dbs_thalamus_f,0) - 0.5)*dbs_thalamus_a*2.5);
        VL_extensor(i) = (((cerebellum_VL_weight*rand*sum(cerebellar_antagonist*glutamine_mod)) +...
            (BG_VL_weight*rand*sum(bg_gpi_fire_antagonist) *gaba_mod)) * dbs_thalamus_effect > (1-neuron_fire));
        VL_extensor(i) = VL_extensor(i) ||((rand*dbs_thalamus_a*0.08)>=(direction_dbs_logistic_curve(dbs_thalamus_f,1)));
        
        for k = 1:10
            % CM_diffuse_single_tmp(k) = (((cerebellum_CM_weight*rand*(sum(cerebellar_agonist)+sum(cerebellar_antagonist))*glutamine_mod) + ...
            %     (BG_CM_weight*rand*sum(bg_gpi_fire)*gaba_mod))...
            %     > (1-neuron_fire))&&(rand>thalamus_dbs_i);
            % CM_diffuse_single_tmp(k) = CM_diffuse_single_tmp(k) ||(rand>(1-thalamus_dbs_e));
            dbs_thalamus_effect = 1+((rand + direction_dbs_logistic_curve(dbs_thalamus_f,0) - 0.5)*dbs_thalamus_a*2.5);
            CM_diffuse_agonist_tmp(k) = (((cerebellum_CM_weight*rand*sum(cerebellar_agonist)*glutamine_mod) + ...
                (BG_CM_weight*rand*sum(bg_gpi_fire_agonist)*gaba_mod)...
                + (rand*sum(L5_flexor_prev)*L5_CM_same_weight*glutamine_mod)...
                + (rand*sum(L5_extensor_prev)*L5_CM_opp_weight*glutamine_mod)) * dbs_thalamus_effect) ...
                > (1-neuron_fire);
            CM_diffuse_agonist_tmp(k) = CM_diffuse_agonist_tmp(k) || ((rand*dbs_thalamus_a*0.08)>=(direction_dbs_logistic_curve(dbs_thalamus_f,1)));

            dbs_thalamus_effect = 1+((rand + direction_dbs_logistic_curve(dbs_thalamus_f,0) - 0.5)*dbs_thalamus_a*2.5);
            CM_diffuse_antagonist_tmp(k) = (((cerebellum_CM_weight*rand*sum(cerebellar_antagonist)*glutamine_mod) + ...
                (BG_CM_weight*rand*sum(bg_gpi_fire_antagonist)*gaba_mod) ...
                + (rand*sum(L5_extensor_prev)*L5_CM_same_weight*glutamine_mod)...
                + (rand*sum(L5_flexor_prev)*L5_CM_opp_weight*glutamine_mod)) * dbs_thalamus_effect)...
                > (1-neuron_fire);
            CM_diffuse_antagonist_tmp(k) = CM_diffuse_antagonist_tmp(k) || ((rand*dbs_thalamus_a*0.08)>=(direction_dbs_logistic_curve(dbs_thalamus_f,1)));
        end

        % CM_diffuse_single(i) = mean(CM_diffuse_single_tmp);
        CM_diffuse_agonist(i) = mean(CM_diffuse_agonist_tmp);
        CM_diffuse_antagonist(i) = mean(CM_diffuse_antagonist_tmp);
        
        
        for g = 1:10
            dbs_L4_effect = (rand + direction_dbs_logistic_curve(dbs_L4_f,0))*dbs_L4_a;
            L4_flexor(1,g) = ((VL_L4_weight*rand*VL_flexor(i)*glutamine_mod + dbs_L4_effect) > (1-neuron_fire));
            L4_flexor(1,g) = L4_flexor(1,g) || ((rand*dbs_L4_a*0.5)>=(direction_dbs_logistic_curve(dbs_L4_f,1)));

            dbs_L4_effect = (rand + direction_dbs_logistic_curve(dbs_L4_f,0))*dbs_L4_a;
            L4_extensor(1,g) = ((VL_L4_weight*rand*VL_extensor(i)*glutamine_mod + dbs_L4_effect) > (1-neuron_fire));
            L4_extensor(1,g) = L4_extensor(1,g) || ((rand*dbs_L4_a*0.5)>=(direction_dbs_logistic_curve(dbs_L4_f,1)));
        end
        
        L23_flexor = zeros(1,4);
        L23_extensor = zeros(1,4);
        for g = 1:4
        
            motor_drive_fire = rand >= (1-motor_drive);
    
            dbs_L23_effect = (rand + direction_dbs_logistic_curve(dbs_L23_f,0))*dbs_L23_a;
            PV_flexor = (L4_PV_weight*rand*sum(L4_flexor,"all")*glutamine_mod) > (1-neuron_fire);
            L23_flexor(g) = (((rand*motor_L23_weight*motor_drive_fire*glutamine_mod)...
                +(rand*L4_L23_weight*sum(L4_flexor)*glutamine_mod)...
                +(rand*PV_flexor*PV_L23_weight*gaba_mod) + dbs_L23_effect) > (1-neuron_fire));
            L23_flexor(g) = L23_flexor(g) || ((rand*dbs_L23_a*0.5)>=(direction_dbs_logistic_curve(dbs_L23_f,1)));

            dbs_L23_effect = (rand + direction_dbs_logistic_curve(dbs_L23_f,0))*dbs_L23_a;
            PV_extensor = (L4_PV_weight*rand*sum(L4_extensor,"all")*glutamine_mod) > (1-neuron_fire);
            L23_extensor(g) = (((rand*motor_L23_weight*motor_drive_fire*glutamine_mod)...
            +(rand*L4_L23_weight*sum(L4_extensor)*glutamine_mod)...
            +(rand*PV_extensor*PV_L23_weight*gaba_mod)+dbs_L23_effect) > (1-neuron_fire));
            L23_extensor(g) = L23_extensor(g) ||((rand*dbs_L23_a*0.5)>=(direction_dbs_logistic_curve(dbs_L23_f,1)));
        end
        
        dbs_L5_effect = (rand + direction_dbs_logistic_curve(dbs_L5_f,0))*dbs_L5_a;

        L5_flexor(i) = ((rand*L23_L5_weight*sum(L23_flexor,"all")*glutamine_mod + dbs_L5_effect) > (1-neuron_fire));
        L5_flexor(i) = L5_flexor(i) || ((rand*dbs_L5_a*0.5)>=(direction_dbs_logistic_curve(dbs_L5_f,1)));
        
        dbs_L5_effect = (rand + direction_dbs_logistic_curve(dbs_L5_f,0))*dbs_L5_a;

        L5_extensor(i) = ((rand*L23_L5_weight*sum(L23_extensor,"all")*glutamine_mod + dbs_L5_effect) > (1-neuron_fire));
        L5_extensor(i) = L5_extensor(i) || ((rand*dbs_L5_a*0.5)>=(direction_dbs_logistic_curve(dbs_L5_f,1)));
    end
    CM_nucleus.CM_agonist = mean(CM_diffuse_agonist);
    CM_nucleus.CM_antagonist = mean(CM_diffuse_antagonist);
end

