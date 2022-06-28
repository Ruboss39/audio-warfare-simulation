classdef matrix_array
    %MATRIX_ARRAY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        r_prime
        row_elements
        column_elements
        elements
        uni_distance
    end
    
    methods
        function obj = matrix_array(r_a,element_distance,row_elements,column_elements)
            %MATRIX_ARRAY Construct an instance of this class
            %   Detailed explanation goes here
            % 
            % Generates a coordinates for a 2D array. Uni_distance is the
            % uniform distance between the array-elements
            %
            % r_a is the vector pointing to the middle of the array.
            %
            %
            % r_prime is a vector on the form:
            %             x_1   x_2 ... x_n
            % r_prime =   y_1   y_2 ... y_n, where n = rows*columns
            %             z_1   z_2 ... z_n
            %

            r_prime = zeros(3,row_elements*column_elements);

            %matrix array in the xy-plane
            uni_distance = element_distance;      % Give length in mm
            element_index = 1;
            for i = 1:row_elements
                for j = 1:column_elements
                    r_prime(1,element_index) = i*uni_distance + r_a(1);
                    r_prime(2,element_index) = j*uni_distance + r_a(2);
                    element_index = element_index + 1;
                end
            end
            
            %set origin of the matrix to (0,0)
            r_prime(1,:) = r_prime(1,:)-row_elements*uni_distance/2 - uni_distance/2;
            r_prime(2,:) = r_prime(2,:)-column_elements*uni_distance/2 - uni_distance/2;

            obj.r_prime = r_prime;
            obj.row_elements = row_elements;
            obj.column_elements = column_elements;
            obj.elements =  row_elements * column_elements;
            obj.uni_distance = element_distance;
        end
    end
end

