function weight = adaptive_array_config(matrix_array,frequency,c)
%
%   Decice how many elements to use in order to obtain the best
%   directivity. Very crude function. Weighting-chart desgined from
%   beamforming simulations listening in the broadside direction. More
%   optimal weights exist for a different direction than broadside
%
%
row_elements = matrix_array.row_elements;
column_elements = matrix_array.column_elements;

uni_distance = matrix_array.uni_distance;

lambda = c/frequency;

lambda_rel = uni_distance/lambda;

if lambda_rel > 0.1581
    % mode = 1
    weight = ones(1,row_elements*column_elements);
    return
elseif 0.156>= lambda_rel && lambda_rel > 0.0986
    mode =3;
elseif 0.0986 >= lambda_rel && lambda_rel > 0.085
    mode = 5;
elseif 0.085 >= lambda_rel && lambda_rel > 0.07
    mode = 6;
else
    mode = 7;
end

weight = zeros(1,row_elements*column_elements);
row_lim = ceil((row_elements)/mode);
column_lim = ceil((column_elements)/mode);
for i = 1:row_lim
    for j = 1:column_lim
        element_index = (mode*(i-1))*row_elements + mode*(j-1) +1;
        weight(element_index) = 1;
    end
end


end

