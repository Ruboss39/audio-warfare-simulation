function locations = find_sources(precision_map,x_vec,y_vec,threshold)
%FIND_SOURCES Summary of this function goes here
%   Detailed explanation goes here

scaling_factor = 1/max(max(precision_map));
map_scaled = scaling_factor*precision_map;

rows = length(precision_map(:,1));
columns = length(precision_map(1,:));

for i = 1:rows
    for j = 1:columns
        if map_scaled(i,j) < threshold
            map_scaled(i,j) = 0;
        end
    end
end

islands_cells_found = 0;
island_cell_count = 0;
coordinates_checked = zeros(rows,columns);

for row = 1:rows
    for column = 1:columns
        if map_scaled(row,column) > 0
            if island_cell_count == 0
                islands_cells_found = island(x_vec(column),y_vec(row),row,column,...
                    map_scaled(row,column));
                island_cell_count = island_cell_count +1;
                islands_cells_found(island_cell_count).get_neighbours(map_scaled);

            else 
            islands_cells_found = [islands_cells_found island(x_vec(column),y_vec(row),...
                row,column,map_scaled(row,column))];
            island_cell_count = island_cell_count +1;
            islands_cells_found(island_cell_count).get_neighbours(map_scaled);
            end
        end
    end
end

tag_color = 1;

for cell_ind = 1:length(islands_cells_found)
    if islands_cells_found(cell_ind).tag(islands_cells_found,tag_color) == 1
        tag_color = tag_color +1;
    end 
end
locations = zeros(tag_color-1,2);


for location_ind = 1:tag_color-1
    temp1 = [];
    temp2 = [];
    temp_elements = 0;
    x_median = 0;
    y_median = 0;
    max_intensity = 0;
    for cell_ind = 1:length(islands_cells_found)
        if (islands_cells_found(cell_ind).island_tag == location_ind)
            x_median = x_median + islands_cells_found(cell_ind).island_coordinates(1)*...
                islands_cells_found(cell_ind).island_intensity;
            temp1 = [temp1 islands_cells_found(cell_ind).island_coordinates(1)];

            y_median = y_median + islands_cells_found(cell_ind).island_coordinates(2)*...
                islands_cells_found(cell_ind).island_intensity;
            temp_elements = temp_elements +1;
            temp2 = [temp2 islands_cells_found(cell_ind).island_coordinates(2)];

            if islands_cells_found(cell_ind).island_intensity > max_intensity
                max_intensity = islands_cells_found(cell_ind).island_intensity;
            end
        end
    end
    %x_median = x_median/(max_intensity*temp_elements);
    %y_median = y_median/(max_intensity*temp_elements);
    x_median = mean(temp1);
    y_median = mean(temp2);
    locations(location_ind,1) = x_median;
    locations(location_ind,2) = y_median;
end

end

