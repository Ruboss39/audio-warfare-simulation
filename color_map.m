classdef color_map
    %COLOR_MAP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        color_data_matrix
    end
    
    methods
        function obj = color_map(data_matrix)
            %COLOR_MAP Construct an instance of this class
            %   Detailed explanation goes here
            obj.color_data_matrix = data_matrix;
        end
    end
end

