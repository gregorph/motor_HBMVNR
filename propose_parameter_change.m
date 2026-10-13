function [candidate, change_made, draws] = propose_parameter_change( ...
    current, parameter_index, distribution_details, dbs_indices, max_draws)
%PROPOSE_PARAMETER_CHANGE One normal increment; reject nonpositive DBS values.
% The normal distribution describes a CHANGE, as in the original algorithm.
% For DBS, accept only current_amplitude + change > 0. Negative changes remain
% possible when they reduce an existing amplitude without crossing zero.
% Initial zero amplitudes still mean no DBS. No clipping/absolute value is used.

if nargin < 5
    max_draws = 1000000;
end
mu = distribution_details(parameter_index, 1);
sigma = distribution_details(parameter_index, 2);
validateattributes(mu, {'numeric'}, {'scalar', 'real', 'finite'});
validateattributes(sigma, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
is_dbs = any(parameter_index == dbs_indices);
for draws = 1:max_draws
    change_made = mu + sigma*randn;
    proposed_value = current(parameter_index) + change_made;
    if isfinite(proposed_value) && (~is_dbs || proposed_value > 0)
        candidate = current;
        candidate(parameter_index) = proposed_value;
        return
    end
end
% Do not silently substitute a clipped value or alter the learned distribution.
error('OGL:RejectionLimit', ['No admissible proposal after %d draws. ' ...
    'Check the mean/SD for parameter %d or increase max_proposal_draws.'], ...
    max_draws, parameter_index);
end
