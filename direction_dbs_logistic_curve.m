function [output_weight] = direction_dbs_logistic_curve(frequency_dbs,output_type)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(output_type == 0) %% returns bias direction for DBS
    K =  0.5;          % upper asymptote / plateau
    A = -0.5;          % lower asymptote / plateau
    R = K - A;          % range = 1.9

    % --- Shape knobs (tune these) ---
    nu = 3.0;           % Richards skew/shape: larger => "later/sharper" transition
    B  = 0.1;         % steepness: smaller => flatter near plateaus

    % Enforce y(50) = 0 exactly
    M = 50 - log(2^nu - 1)/B;
else
    K =  0.974;          % upper asymptote / plateau
    A = 0.015;          % lower asymptote / plateau
    R = K - A;          % range = 1.9

    % --- Shape knobs (tune these) ---
    nu = 3.0;           % Richards skew/shape: larger => "later/sharper" transition
    B  = 0.1;         % steepness: smaller => flatter near plateaus

    % Enforce y(50) = 0 exactly
    M = 50 - log(2^nu - 1)/B;
end

    y = @(x) A + R ./ (1 + exp(B .* (x - M))).^(1./nu);
    
    output_weight = y(frequency_dbs);
end