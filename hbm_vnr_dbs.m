function [deepcerebelarfire, golgi_cellfire, purkinjefire, granulefire, ...
    interneuronfire, cerebellum_lp] = hbm_vnr_dbs(infolive_firing, ...
    pontine_firing, golgi_cellfire, mod_inhibition, mod_excitation, ...
    neuron_firing, weight_params)
%HBM_VNR_DBS Original cerebellar circuit with DBS applied to the DCN only.
% The first six inputs retain the meaning used in hbm_vnr. Optional fields:
%   weight_params.dbs_DCN_f                 Frequency in Hz (default 130).
%   weight_params.dbs_DCN_a                 Amplitude >= 0 (default 0).
%   weight_params.dbs_DCN_activation_scale  Independent activation (0.08).
% Zero amplitude preserves the original first five outputs and RNG draws.
% The original sixth output was never assigned; it is returned as [] here.
% Requires the project's direction_dbs_logistic_curve(frequency, mode).

if nargin < 6 || isempty(neuron_firing)
    neuron_firing = 0.8;
end
if nargin < 7 || isempty(weight_params)
    weight_params = struct();
end
dbs_DCN_f = 130;
dbs_DCN_a = 0;
dbs_DCN_activation_scale = 0.08;
if isfield(weight_params, 'dbs_DCN_f')
    dbs_DCN_f = weight_params.dbs_DCN_f;
end
if isfield(weight_params, 'dbs_DCN_a')
    dbs_DCN_a = weight_params.dbs_DCN_a;
end
if isfield(weight_params, 'dbs_DCN_activation_scale')
    dbs_DCN_activation_scale = weight_params.dbs_DCN_activation_scale;
end
% Called for every circuit/sample: keep these scalar checks inexpensive.
if ~(isscalar(dbs_DCN_a) && isreal(dbs_DCN_a) && ...
        isfinite(dbs_DCN_a) && dbs_DCN_a >= 0)
    error('HBM:InvalidAmplitude', 'dbs_DCN_a must be finite and nonnegative.');
end
if ~(isscalar(dbs_DCN_f) && isreal(dbs_DCN_f) && ...
        isfinite(dbs_DCN_f) && dbs_DCN_f >= 0)
    error('HBM:InvalidFrequency', 'dbs_DCN_f must be finite and nonnegative.');
end
if ~(isscalar(dbs_DCN_activation_scale) && isreal(dbs_DCN_activation_scale) && ...
        isfinite(dbs_DCN_activation_scale) && dbs_DCN_activation_scale >= 0)
    error('HBM:InvalidScale', 'dbs_DCN_activation_scale must be nonnegative.');
end

% Original input sampling and upstream cerebellar equations, in the same order.
pontinefire = rand > (1-pontine_firing);
infolivefire = rand > (1-infolive_firing);
pontine_granule_mod = rand*(1+mod_excitation)*pontinefire;
golgi_granule_mod = rand*(1+mod_inhibition)*golgi_cellfire;
granulefire = (pontine_granule_mod-golgi_granule_mod) > (1-neuron_firing);

granule_golgi_mod = rand*(1+mod_excitation)*granulefire;
golgi_cellfire = granule_golgi_mod > (1-neuron_firing);
granule_ii_mod = rand*(1+mod_excitation)*granulefire;
interneuronfire = granule_ii_mod > (1-neuron_firing);

infolivemodulation = rand*(1+mod_excitation)*infolivefire;
granulemodulation = rand*(1+mod_excitation)*granulefire;
ii_purkinje_mod = rand*(1+mod_inhibition)*interneuronfire;
purkinjetmp = infolivemodulation + granulemodulation - ii_purkinje_mod;
purkinjefire = purkinjetmp > (1-neuron_firing);

purkinjedcnmodulation = rand*purkinjefire*(1+mod_inhibition);
infolivedcnmodulation = rand*infolivefire*(1+mod_excitation);
pontinedcnmodulation = rand*pontinefire*(1+mod_excitation);
dcn_input = -purkinjedcnmodulation + infolivedcnmodulation + pontinedcnmodulation;

if dbs_DCN_a > 0
    % Same multiplicative input rule as the supplied BG and VL DBS code.
    dbs_DCN_effect = 1 + (rand + ...
        direction_dbs_logistic_curve(dbs_DCN_f, 0) - 0.5)*dbs_DCN_a*2.5;
    deepcerebelarfire = (dcn_input*dbs_DCN_effect) > (1-neuron_firing);
    % Independent activation uses a fresh draw, not the input-modulation draw.
    independent_activation = rand*dbs_DCN_a*dbs_DCN_activation_scale >= ...
        direction_dbs_logistic_curve(dbs_DCN_f, 1);
    deepcerebelarfire = deepcerebelarfire || independent_activation;
else
    deepcerebelarfire = dcn_input > (1-neuron_firing);
end
cerebellum_lp = [];
end
