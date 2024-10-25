% Step 1: Define the directory and file name for the CSV and XML files
directory = 'C:\Users\12826\Documents\MATLAB\9757afe8-3bee-4586-80ea-7dc2a1f2f01c';
caseName = 'postablationmap';
csvDirectory = fullfile(directory, caseName, 'Contact_Mapping');
csvFileName = 'Map_PP_omni.csv';
csvFullFilePath = fullfile(csvDirectory, csvFileName);

% Define path for the Contact Mapping Model XML file
contactMappingModelPath = fullfile(directory, caseName, 'Contact_Mapping_Model.xml');

% Step 2: Load the CSV file and extract surfaceX, surfaceY, surfaceZ, and MaxBipolarVoltage
omniData = readtable(csvFullFilePath);
surfaceX = omniData.surfaceX;
surfaceY = omniData.surfaceY;
surfaceZ = omniData.surfaceZ;
Maxpp_Vmax = omniData.pp_Vmax;

% Step 3: Combine the surface coordinates into a single matrix and remove rows with NaN or zeros
points = [surfaceX, surfaceY, surfaceZ];

% Remove rows where surfaceX, surfaceY, or surfaceZ is NaN or zero
validRows = all(~isnan(points) & points ~= 0, 2);  % Find rows with valid values
points = points(validRows, :);  % Keep only valid rows
Maxpp_Vmax = Maxpp_Vmax(validRows);  % Also filter voltage data

% Step 4: Generate KD-tree using valid points
MdlKDT = KDTreeSearcher(points);

% Define the search radius (1.5 mm)
searchRadius = 1.5;

% Prepare a new matrix to store the newly added points
newPoints = [];
newVoltages = [];

% Step 5: Iterate over each point and find neighbors using the KD-tree
for i = 1:size(points, 1)
    % Find the points within the search radius
    idx = rangesearch(MdlKDT, points(i, :), searchRadius);
    
    % Extract the indices of neighboring points
    neighborsIdx = idx{1};
    
    if length(neighborsIdx) > 1  % Check if there are multiple neighbors
        % Calculate the maximum voltage between neighbors
        maxVoltage = max(Maxpp_Vmax(neighborsIdx));
        
        % Find the midpoints between the current point and each of its neighbors
        for j = 2:length(neighborsIdx)
            midpoint = (points(i, :) + points(neighborsIdx(j), :)) / 2;
            newPoints = [newPoints; midpoint];  % Store the new points
            newVoltages = [newVoltages; maxVoltage];  % Assign max voltage to the new point
        end
    end
end

% Step 6: Combine original points with new points for the enhanced resolution
enhancedPoints = [points; newPoints];
enhancedVoltages = [Maxpp_Vmax; newVoltages];

% Step 7: Load and process the XML file
if ~isfile(contactMappingModelPath)
    error('File does not exist: %s', contactMappingModelPath);
end

try
    dxgeo = readstruct(contactMappingModelPath);
    disp('XML file loaded successfully.');
catch ME
    error('Failed to load XML file: %s', ME.message);
end

% Extract data from the XML structure
try
    Vertices_LA = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Vertices.Text));
    Polygons_LA = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Polygons.Text));
    MapColor_LA = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Map_color.Text));
    Normals_LA = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Normals.Text));
    MapStatus_LA = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Map_status.Text));
    ColorHighLow = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Color_high_low.Text));
    SurfaceOfOrigin = str2double(strsplit(dxgeo.DIFBody.Volumes.Volume.Surface_of_origin.Text));
catch ME
    error('Failed to extract data from the XML structure: %s', ME.message);
end

% Expected numbers from the XML structure
numVertices = 62722;  % Number of vertices as per XML
numPolygons = 125436; % Number of polygons as per XML

% Reshape vertices, polygons, and normals into Nx3 matrices
Vertices_LA = reshape(Vertices_LA, 3, []).';  % Convert to Nx3 matrix
Polygons_LA = reshape(Polygons_LA, 3, []).';  % Convert to Nx3 matrix
Normals_LA = reshape(Normals_LA, 3, []).';    % Convert to Nx3 matrix

% Ensure MapColor_LA is Nx1
MapColor_LA = MapColor_LA(:);  % Convert to Nx1 if necessary

% Filter vertices based on MapStatus_LA (optional: only plot active points)
activeIndices = MapStatus_LA == 1;
Vertices_Active = Vertices_LA(activeIndices, :);
MapColor_Active = MapColor_LA(activeIndices);

% Step 8: Custom Colormap Definition
colors = [
    0.5 0 0.5;    % Purple for 10 mV
    0 0 1;        % Blue for 0.5 mV
    0 1 0;        % Green for 0.4 mV (Changed)
    1 1 0;        % Yellow for 0.3 mV (Changed)
    1 0 0;        % Red for 0.1 mV
    0.5 0.5 0.5;  % Gray for 0.05 mV or less
];

% Define the color levels manually (non-linearly spaced)
color_levels = [10, 0.5, 0.4, 0.3, 0.1, 0.05];  % Adjusted levels for new colors

% Normalize the color levels for colormap (0 to 1 range)
color_levels_norm = (color_levels - min(color_levels)) / (max(color_levels) - min(color_levels));

% Interpolate the colors
nColors = 256;  % Number of colors in the colormap
xq = linspace(0, 1, nColors);
custom_cmap = interp1(color_levels_norm, colors, xq, 'linear');

% Step 9: Visualize the original geometry from XML and the enhanced points with scatter3
figure;

% Visualize the geometry with Map Color data
%p = patch('Vertices', Vertices_LA, 'Faces', Polygons_LA, 'FaceVertexCData', MapColor_LA, 'FaceColor', 'interp');
% Plot the high-resolution geometry using patch
shellColour = [0.9020    0.8549    0.7294]
p = patch('Vertices', Vertices_LA, 'Faces', Polygons_LA,'FaceColor', shellColour);
set(p, 'EdgeColor', 'none');
daspect([1 1 1]);
camlight; 
lighting gouraud;
axis tight; 
grid on;
alpha(0.85);
hold on;

% Scatter plot with color based on Maxpp_Vmax values
scatter3(surfaceX(1:end-1), surfaceY(1:end-1), surfaceZ(1:end-1), 50, Maxpp_Vmax, 'filled');  % Use Maxpp_Vmax for color mapping

% Apply the custom colormap and color settings
colormap(custom_cmap);
caxis([ColorHighLow(1) ColorHighLow(2)]);  % Set color axis range based on ColorHighLow

% Add colorbar with custom ticks and labels
colorbar('Ticks', [0.05, 0.1, 0.3, 0.4, 0.5, 10], ...
         'TickLabels', {'0.05 mV', '0.1 mV', '0.3 mV', '0.4 mV', '0.5 mV', '10 mV'});

% Step 9: Visualize the original geometry from XML and the enhanced points with scatter3
figure;

% Plot the high-resolution geometry using patch (with interpolated MapColor data)
p = patch('Vertices', Vertices_LA, 'Faces', Polygons_LA,'FaceColor', shellColour);
set(p, 'EdgeColor', 'none');  % Remove edge lines for a smooth surface
daspect([1 1 1]);  % Ensure aspect ratio is equal for all axes
camlight;  % Add a camera light to the figure
lighting gouraud;  % Apply smooth lighting
axis tight;  % Set axis limits tight to the data
grid on;  % Display grid
alpha(0.85);  % Set transparency to 85% for the surface
hold on;

% Scatter plot with color based on Maxpp_Vmax values (for enhanced points)
scatter3(surfaceX(1:end-1), surfaceY(1:end-1), surfaceZ(1:end-1),50, Maxpp_Vmax, 'filled');  % Size 50, use Maxpp_Vmax for color
saveas(gcf, 'original_3D_Shell_Visualization .fig');  % Save as .fig
% Apply the custom colormap and color settings
colormap(custom_cmap);
caxis([ColorHighLow(1) ColorHighLow(2)]);  % Set color axis range based on ColorHighLow

% Add colorbar with custom ticks and labels
h = colorbar('Ticks', [0.05, 0.1, 0.3, 0.4, 0.5, 10], ...
             'TickLabels', {'0.05 mV', '0.1 mV', '0.3 mV', '0.4 mV', '0.5 mV', '10 mV'});

% Label colorbar
ylabel(h, 'Voltage (mV)');

% Add title to the plot
title('High-Resolution Geometry Visualization with Max Voltages and Scatter Plot');

% Add axis labels
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');




% Step 1: Visualize the shell using patch

figure;

% Plot the high-resolution shell using patch with interpolated voltages
p = patch('Vertices', Vertices_LA, 'Faces', Polygons_LA,'FaceColor', shellColour);
set(p, 'EdgeColor', 'none');  % Remove edge lines for a smooth surface
daspect([1 1 1]);  % Ensure equal scaling for all axes
camlight;  % Add a light source for better 3D visualization
lighting gouraud;  % Apply smooth shading for better visualization
axis tight;  % Fit axes tightly around the data
grid on;  % Show grid
alpha(0.85);  % Set transparency of the surface to 85%

hold on;  % Keep the current plot active so that we can add scatter points

% Step 2: Scatter plot the enhanced points with color based on enhancedVoltages

% Scatter plot of the enhanced points
scatter3(enhancedPoints(:,1), enhancedPoints(:,2), enhancedPoints(:,3), 50, enhancedVoltages, 'filled', 'MarkerEdgeColor', 'k');
% 50 is the marker size, and 'filled' makes the points solid color

% Step 3: Apply the custom colormap and color settings
colormap(custom_cmap);  % Apply the custom colormap defined earlier
caxis([ColorHighLow(1) ColorHighLow(2)]);  % Set color axis range based on voltage range

% Step 4: Add colorbar with custom ticks and labels
h = colorbar('Ticks', [0.05, 0.1, 0.3, 0.4, 0.5, 10], ...
             'TickLabels', {'0.05 mV', '0.1 mV', '0.3 mV', '0.4 mV', '0.5 mV', '10 mV'});

% Label the colorbar
ylabel(h, 'Voltage (mV)');

% Step 5: Add title and axis labels
title('3D Shell Visualization with Enhanced Points and Voltage Mapping');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');

hold off;  % Release the hold on the plot
