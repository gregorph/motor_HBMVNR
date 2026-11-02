function [value] = put_in_range(value,minimum,maximum)
%PUT_IN_RANGE Summary of this function goes here
%   Detailed explanation goes here
if(value<minimum)
    value = minimum;
end
if(value>maximum)
    value = maximum;
end
end

