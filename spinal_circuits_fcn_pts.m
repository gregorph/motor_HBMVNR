function [spinal_circuit_output] = spinal_circuits_fcn_pts(L5_flexor,L5_extensor,muscleposition,weight_params)
%SPINAL_CIRCUITS_FCN_PTS Summary of this function goes here
%   Detailed explanation goes here
    
        neuron_fire = 0.8;
    
        %Neurotransmitter settings
        mod_inhibition = -0.5; %Glycine
        mod_excitation = 0.5; % Glutamine
    
        %Synaptic weights
        L5_SpI_weight = 2;
        L5_Ia_weight = 1;
        proprioceptor_Ia_weight = 2; %0.4
        Spi_a_weight = 2;
        Ia_a_weight = 1.2;
        a_motorfibre_weight = 1;
    

    var_list = { ...
    'neuron_fire','mod_inhibition','mod_excitation', ...
    'L5_SpI_weight','L5_Ia_weight','proprioceptor_Ia_weight','Spi_a_weight','Ia_a_weight','a_motorfibre_weight', ...
     };

for i = 1:numel(var_list)
    v = var_list{i};
    if isfield(weight_params, v)
        eval([v ' = weight_params.(v);']);
    end
end


    muscleposition_input = 0.04 + 0.32*muscleposition;
    muscleposition_input_rev = (0.04 + 0.32) - muscleposition_input;
    
    flexorpower = zeros(length(L5_flexor), 10);
    extensorpower = zeros(length(L5_extensor), 10);
    for i = 1:length(L5_flexor)
        
        Sp_I_flexor = zeros(1,9);
        Sp_I_extensor = zeros(1,9);
        for g = 1:9
            Sp_I_flexor(1,g) = (rand*mod_excitation*L5_SpI_weight*sum(L5_flexor(i))) > (1-neuron_fire);
            Sp_I_extensor(1,g) = (rand*mod_excitation*L5_SpI_weight*sum(L5_extensor(i))) > (1-neuron_fire);
        end

        muscleproprioceptor_extensor = zeros(1, 5);
        muscleproprioceptor_flexor = zeros(1, 5);

        Ia_flexor = zeros(1, 3);
        Ia_extensor = zeros(1, 3);
        for g = 1:3
            for z = 1:5
                muscleproprioceptor_extensor(1,z) = rand > (1-muscleposition_input_rev);
                muscleproprioceptor_flexor(1,z) = rand > (1-muscleposition_input);
            end
            Ia_flexor(1,g) = (rand*mod_excitation*L5_Ia_weight*sum(L5_flexor(i))) + rand*mod_excitation*proprioceptor_Ia_weight*sum(muscleproprioceptor_flexor,"all") > (1-neuron_fire);
            Ia_extensor(1,g) = (rand*mod_excitation*L5_Ia_weight*sum(L5_extensor(i))) + rand*mod_excitation*proprioceptor_Ia_weight*sum(muscleproprioceptor_extensor,"all") > (1-neuron_fire);
        end
    
        alpha_flexor = zeros(1,3);
        alpha_extensor = zeros(1,3);
        for g = 1:3
            alpha_flexor(1,g) = (rand*mod_excitation*Spi_a_weight*sum(Sp_I_flexor,"all")) + (rand*mod_inhibition*Ia_a_weight*sum(Ia_extensor,"all")) > (1-neuron_fire);
            alpha_extensor(1,g) = (rand*mod_excitation*Spi_a_weight*sum(Sp_I_extensor,"all")) + (rand*mod_inhibition*Ia_a_weight*sum(Ia_flexor,"all")) > (1-neuron_fire);
        end
    
        for g = 1:10
            flexorpower(i,g) = (rand*mod_excitation*a_motorfibre_weight*sum(alpha_flexor,"all")) > (1-neuron_fire);
            extensorpower(i,g) = (rand*mod_excitation*a_motorfibre_weight*sum(alpha_extensor,"all")) > (1-neuron_fire);
        end
    end
    spinal_circuit_output.flexorpower = flexorpower;
    spinal_circuit_output.extensorpower = extensorpower;
end

