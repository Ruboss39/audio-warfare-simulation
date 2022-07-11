classdef island < handle
    %ISLAND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        island_coordinates
        island_index
        island_intensity
        island_neighbours_index
        island_tag
    end
    
    methods
        function obj = island(x,y,row,column,intensity)
            %ISLAND Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
            obj.island_coordinates = [x,y];
            obj.island_index = [row,column];
            obj.island_intensity = intensity;
            end
        end

        function tagged = tag(obj,cell_islands,tag_color)
            tagged = 0;
            if isempty(obj.island_tag ) == 1
                obj.island_tag = tag_color;
                for cell_island = 1:length(cell_islands)
                    temp_index = cell_islands(cell_island).island_index;
                    if isempty(obj.island_neighbours_index) == 0
                        for i = 1:length(obj.island_neighbours_index(:,1))
                            if temp_index(1) == obj.island_neighbours_index(i,1) &&...
                                    temp_index(2) == obj.island_neighbours_index(i,2)
                                if isempty(cell_islands(cell_island).island_tag) == 1
                                    cell_islands(cell_island).tag(cell_islands,tag_color);
                                    tagged = 1;
                                end
                            end
                        end
                    end
                end
            end
        end

        function get_neighbours(obj,intensity_map)
            for row = -1:1
                for column = -1:1
                    if row == 0 && column == 0
                        continue;
                    end
                    neigh_row = row + obj.island_index(1);
                    neigh_column = column + obj.island_index(2);

                    if neigh_row < 1 || neigh_column < 1 
                        continue;
                    end
                    if neigh_row > length(intensity_map(:,1)) || neigh_column > length(intensity_map(1,:))
                        continue;
                    end
                    if intensity_map(neigh_row,neigh_column) > 0
                        obj.island_neighbours_index  = [obj.island_neighbours_index; [neigh_row,neigh_column]];
                    end
                end
            end
            
        end

    end
end

