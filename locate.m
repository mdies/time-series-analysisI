function [idleft,idright] = locate(value, vector)

% This function locates the position for the closest match of a given value 
% in an array, identifying the indices of the values that would surround it

[d p] = min(abs(vector - value));

if value > vector(p)
    idleft=p;
    idright=p+1;
else
    idright=p;
    idleft=p-1;
end

if idright > length(vector)
    idright=0;
end

end