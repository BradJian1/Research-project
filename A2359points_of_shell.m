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

% Filter points based on the 'utilized' field if it exists
if ismember('utilized', omniData.Properties.VariableNames)
    utilized = omniData.utilized;
    % Select only rows where utilized == 1
    surfaceX = surfaceX(utilized == 1);
    surfaceY = surfaceY(utilized == 1);
    surfaceZ = surfaceZ(utilized == 1);
    Maxpp_Vmax = Maxpp_Vmax(utilized == 1);
end

% Step 3: Combine the surface coordinates into a single matrix and remove rows with NaN or zeros
points = [surfaceX, surfaceY, surfaceZ];
validRows = all(~isnan(points) & points ~= 0, 2);  % Find rows with valid values
points = points(validRows, :);  % Keep only valid rows
Maxpp_Vmax = Maxpp_Vmax(validRows);  % Also filter voltage data

% Step 4: Generate KD-tree using valid points
MdlKDT = KDTreeSearcher(points);
searchRadius = 1.5;  % Define the search radius (1.5 mm)

% Prepare a new matrix to store the newly added points
newPoints = [];
newVoltages = [];

% Step 5: Iterate over each point and find neighbors using the KD-tree
for i = 1:size(points, 1)
    idx = rangesearch(MdlKDT, points(i, :), searchRadius);
    neighborsIdx = idx{1};
    
    if length(neighborsIdx) > 1
        maxVoltage = max(Maxpp_Vmax(neighborsIdx));
        
        % Find the midpoints between the current point and each of its neighbors
        for j = 2:length(neighborsIdx)
            midpoint = (points(i, :) + points(neighborsIdx(j), :)) / 2;
            newPoints = [newPoints; midpoint];
            newVoltages = [newVoltages; maxVoltage];
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

% Reshape vertices, polygons, and normals into Nx3 matrices
Vertices_LA = reshape(Vertices_LA, 3, []).';  % Convert to Nx3 matrix
Polygons_LA = reshape(Polygons_LA, 3, []).';  % Convert to Nx3 matrix

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
color_levels = [10, 0.5, 0.4, 0.3, 0.1, 0.05];
color_levels_norm = (color_levels - min(color_levels)) / (max(color_levels) - min(color_levels));
nColors = 256;
xq = linspace(0, 1, nColors);
custom_cmap = interp1(color_levels_norm, colors, xq, 'linear');

% Step 9: Visualize the original geometry from XML and the enhanced points
figure;
shellColour = [0.9020, 0.8549, 0.7294];
p = patch('Vertices', Vertices_LA, 'Faces', Polygons_LA, 'FaceColor', shellColour);
set(p, 'EdgeColor', 'none');
daspect([1 1 1]);
camlight; 
lighting gouraud;
axis tight; 
grid on;
alpha(0.85);
hold on;

% Scatter plot with color based on Maxpp_Vmax values
scatter3(enhancedPoints(:, 1), enhancedPoints(:, 2), enhancedPoints(:, 3), 50, enhancedVoltages, 'filled');

% Apply the custom colormap and color settings
colormap(custom_cmap);
caxis([ColorHighLow(1) ColorHighLow(2)]);

% Add colorbar with custom ticks and labels
h = colorbar('Ticks', [0.05, 0.1, 0.3, 0.4, 0.5, 10], ...
             'TickLabels', {'0.05 mV', '0.1 mV', '0.3 mV', '0.4 mV', '0.5 mV', '10 mV'});
ylabel(h, 'Voltage (mV)');

% Add title and axis labels
title('High-Resolution Geometry Visualization with Max Voltages and Scatter Plot');
xlabel('X axis');
ylabel('Y axis');
zlabel('Z axis');

% Save the figure as a .fig file
saveas(gcf, 'original_3D_Shell_Visualization.fig');
% Step 9: Visualize the original geometry from XML and the enhanced points with scatter3
figure;
shellColour = [0.9020, 0.8549, 0.7294];
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
scatter3(surfaceX, surfaceY, surfaceZ, 50, Maxpp_Vmax, 'filled');  % Use Maxpp_Vmax for color mapping
% Apply the custom colormap and color settings
colormap(custom_cmap);
caxis([ColorHighLow(1) ColorHighLow(2)]);