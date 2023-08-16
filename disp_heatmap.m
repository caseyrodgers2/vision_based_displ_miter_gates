%% Displacement Error Heatmap
% Note: This script is meant to be post processing for the
% motion_compensate_homography_klt_video.m script
% @author: Casey Rodgers
% @ date:  12/31/2021

close all
clear all
clc

% Input Result Folder Names
blend_folder = 'blender_results\';     % Blender Results Folder name
field_folder = 'field_results\';       % Field Results Folder name

% Load Blender Data
blend_disp_file = horzcat('blender_5x\', blend_folder, '\disp_data_blender_5x.mat');
load(blend_disp_file);              % Load file
blend_avgX_pix = store_avgX_pix;    % Extract avg disp in pixels
blend_avgY_pix = store_avgY_pix;    % Extract avg disp in pixels

% Load ROI data
load('blender_5x\static_roi_store2'); % Load static_roi_store2
load('blender_5x\nonstatic_roi_store3'); % Load nonstatic_roi_store3
roi_store = [static_roi_store, nonstatic_roi_store];  % Combine roi stores

% Load Field Data
field_disp_file = horzcat('fill_13_south\', field_folder, '\disp_data_fill_13_south.mat');
load(field_disp_file);                   % Load file
field_avgX_pix2 = store_avgX_pix;         % Extract avg disp in pixels
field_avgY_pix2 = store_avgY_pix;         % Extract avg disp in pixels
lastFrame = size(field_avgX_pix2, 1);     % Last frame index for field displacements
field_avgX_pix = field_avgX_pix2(lastFrame, :);  % Extract avg disp in pixels
field_avgY_pix = field_avgY_pix2(lastFrame, :);  % Extract avg disp in pixels

% Note: Assume that roi_centroid are equal for both Blender data and Field
% data

% Load Ground Truth Data (FEM data)
ground_disp_file = horzcat('PixelDisplacement_5times.csv');
load(ground_disp_file);                % Load file
ground_avgX_pix = PixelDisplacement_5times(:, 1)';    % Extract avg disp in pixels
ground_avgY_pix = PixelDisplacement_5times(:, 2)';    % Extract avg disp in pixels

ground_avgX_pix(good_rois == 0) = [];   % Get rid of bad roi's
ground_avgY_pix(good_rois == 0) = [];   % Get rid of bad roi's

% Load image
background_img = imread('fill_13_south\DSC_1666.JPG');  % Load background image


%%% Plot Difference
% Blender minus Field
bl_fi_diff_x = blend_avgX_pix - field_avgX_pix;                             % Find difference for x
bl_fi_diff_y = blend_avgY_pix - field_avgY_pix;                             % Find difference for y
%diff_plot(bl_fi_diff_x, bl_fi_diff_y, roi_centroid, 'Blender minus Field')  % Plot difference

% Blender minus Ground Truth
bl_gt_diff_x = blend_avgX_pix(5:size(blend_avgX_pix, 2)) - ground_avgX_pix;        % Find difference for x
bl_gt_diff_y = blend_avgY_pix(5:size(blend_avgY_pix, 2)) - ground_avgY_pix;        % Find difference for y
%diff_plot(bl_gt_diff_x, bl_gt_diff_y, roi_centroid, 'Blender minus Ground Truth')  % Plot difference

% Field minus Ground Truth
fi_gt_diff_x = field_avgX_pix(5:size(field_avgX_pix, 2)) - ground_avgX_pix;      % Find difference for x
fi_gt_diff_y = field_avgY_pix(5:size(field_avgY_pix, 2)) - ground_avgY_pix;      % Find difference for y
%diff_plot(fi_gt_diff_x, fi_gt_diff_y, roi_centroid, 'Field minus Ground Truth')  % Plot difference


% Plot Displacements
% figure
% h1 = axes;
% scatter3(blend_avgX_pix, roi_centroid(:, 1), roi_centroid(:, 2))
% hold on
% scatter3(field_avgX_pix, roi_centroid(:, 1), roi_centroid(:, 2))
% hold on
% scatter3(ground_avgX_pix, roi_centroid(5:size(roi_centroid, 1), 1), roi_centroid(5:size(roi_centroid, 1), 2))
% hold off
% title('X Displacement')
% xlabel('Displacement (pixels)')
% ylabel('Node X Coordinate (pixels)')
% zlabel('Node Y Coordinate (pixels)')
% legend('Blender', 'Field', 'Ground Truth')
% set(h1, 'Zdir', 'reverse')
% set(h1, 'Ydir', 'reverse')
% 
% figure
% h2 = axes;
% scatter3(blend_avgY_pix, roi_centroid(:, 1), roi_centroid(:, 2))
% hold on
% scatter3(field_avgY_pix, roi_centroid(:, 1), roi_centroid(:, 2))
% hold on
% scatter3(ground_avgY_pix, roi_centroid(5:size(roi_centroid, 1), 1), roi_centroid(5:size(roi_centroid, 1), 2))
% hold off
% title('Y Displacement')
% xlabel('Displacement (pixels)')
% ylabel('Node X Coordinate (pixels)')
% zlabel('Node Y Coordinate (pixels)')
% legend('Blender', 'Field', 'Ground Truth')
% set(h2, 'Zdir', 'reverse')
% set(h2, 'Ydir', 'reverse')

% Update roi_store to be larger for illustration purposes
for i = 1:size(nonstatic_roi_store, 2)
    w_2 = nonstatic_roi_store{i}(3) * 2;
    h_2 = nonstatic_roi_store{i}(4) * 2;
    nonstatic_roi_store{i}(1) = nonstatic_roi_store{i}(1) - (w_2 - nonstatic_roi_store{i}(3)) / 2;
    nonstatic_roi_store{i}(2) = nonstatic_roi_store{i}(2) - (h_2 - nonstatic_roi_store{i}(4)) / 2;
    nonstatic_roi_store{i}(3) = w_2;
    nonstatic_roi_store{i}(4) = h_2;
end



% Create Heatmaps
% Displacements
field_mag = sqrt(field_avgX_pix.^2 + field_avgY_pix.^2);  % Field displacement magnitude
output = heat_map(field_mag(5:size(field_mag, 2)), background_img, nonstatic_roi_store, 0, 'Field Displacement Magnitude');
output = heat_map(field_avgX_pix(5:size(field_mag, 2)), background_img, nonstatic_roi_store, 0, 'Field Displacement X');
output = heat_map(field_avgY_pix(5:size(field_mag, 2)), background_img, nonstatic_roi_store, 0, 'Field Displacement Y');

blend_mag = sqrt(blend_avgX_pix.^2 + blend_avgY_pix.^2);  % Field displacement magnitude
output = heat_map(blend_mag(5:size(field_mag, 2)), background_img, nonstatic_roi_store, 0, 'PBGM Displacement Magnitude');
output = heat_map(blend_avgX_pix(5:size(field_mag, 2)), background_img, nonstatic_roi_store, 0, 'PBGM Displacement X');
output = heat_map(blend_avgY_pix(5:size(field_mag, 2)), background_img, nonstatic_roi_store, 0, 'PBGM Displacement Y');

ground_mag = sqrt(ground_avgX_pix.^2 + ground_avgY_pix.^2);  % Field displacement magnitude
output = heat_map(ground_mag, background_img, nonstatic_roi_store, 0, 'FEM Displacement Magnitude');
output = heat_map(ground_avgX_pix, background_img, nonstatic_roi_store, 0, 'FEM Displacement X');
output = heat_map(ground_avgY_pix, background_img, nonstatic_roi_store, 0, 'FEM Displacement Y');

% Differences
bl_gt_mag = sqrt(bl_gt_diff_x.^2 + bl_gt_diff_y.^2);  % Field displacement magnitude
output = heat_map(bl_gt_mag, background_img, nonstatic_roi_store, 1, 'PBGM minus FEM Difference Magnitude');
output = heat_map(bl_gt_diff_x, background_img, nonstatic_roi_store, 1, 'PBGM minus FEM Difference X');
output = heat_map(bl_gt_diff_y, background_img, nonstatic_roi_store, 1, 'PBGM minus FEM Difference Y');

fi_gt_mag = sqrt(fi_gt_diff_x.^2 + fi_gt_diff_y.^2);  % Field displacement magnitude
output = heat_map(fi_gt_mag, background_img, nonstatic_roi_store, 1, 'Field minus FEM Difference Magnitude');
output = heat_map(fi_gt_diff_x, background_img, nonstatic_roi_store, 1, 'Field minus FEM Difference X');
output = heat_map(fi_gt_diff_y, background_img, nonstatic_roi_store, 1, 'Field minus FEM Difference Y');


%%% Plotting Functions

% Difference Plots
% diff_x = difference in x displacement
% diff_y = difference in y displacement
% roi_centroid = coordinates of roi's centroids
% title_name = String of Title Name. ___ X or Y Displacement
function diff_plot(diff_x, diff_y, roi_centroid, title_name)
    
    % If there is difference in sizes, adjust
    size_diff = size(diff_x, 2) - size(roi_centroid, 1);    % Size difference
    size_diff2 = size(diff_x, 2) - size(roi_centroid, 1);    % Size difference

    % If there is a size difference, then we need to fix the roi_centroid
    % to match the size
    if (size_diff < 0)
        roi_centroid_x = roi_centroid(abs(size_diff)+1:size(roi_centroid, 1), 1);
        roi_centroid_y = roi_centroid(abs(size_diff)+1:size(roi_centroid, 1), 2);

    % Otherwise, we're good
    else
        roi_centroid_x = roi_centroid(:, 1);
        roi_centroid_y = roi_centroid(:, 2);
    end

    % Plot
    figure
    h1 = axes;
    scatter3(diff_x, roi_centroid_x, roi_centroid_y)
    title(horzcat(title_name, ' X Displacement'))
    xlabel('Difference (pixels)')
    ylabel('Node X Coordinate (pixels)')
    zlabel('Node Y Coordinate (pixels)')
    set(h1, 'Zdir', 'reverse')
    set(h1, 'Ydir', 'reverse')
    
    figure
    h2 = axes;
    scatter3(diff_y, roi_centroid_x, roi_centroid_y)
    title(horzcat(title_name, ' Y Displacement'))
    xlabel('Difference (pixels)')
    ylabel('Node X Coordinate (pixels)')
    zlabel('Node Y Coordinate (pixels)')
    set(h2, 'Zdir', 'reverse')
    set(h2, 'Ydir', 'reverse')
end



% Heat Map
% displ = displacements 
% img = image you want heat map on
% roi_store = where displacement heat maps will be plotted
% data_type = 0 - Displacement. 1 - Difference.
% title_str = title string
function output = heat_map(displ, img, roi_store, data_type, title_str)

    % Find max and min values
    displ_wo_nan = displ;       % Copy of displacement matrix
    displ_wo_nan(isnan(displ)) = [];   % Displacement matrix without nan
    max_displ = max(displ_wo_nan);     % Max displacement value
    min_displ = min(displ_wo_nan);     % Min displacement value

    data_type_str = '';

    if data_type == 0
        data_type_str = 'Displ.';

    elseif data_type == 1
        data_type_str = 'Diff.';
    end

    max_str = horzcat('Max ', data_type_str, ': ', num2str(round(max_displ, 2)));  % Max Displ label
    min_str = horzcat('Min ', data_type_str, ': ', num2str(round(min_displ, 2)));  % Min Displ label

    % Color max and min
    max_color = 240 / 360;      % Blue
    min_color = 0;              % Red

    output = img;               % Image
    min_pos = zeros(2, 1);      % Min position
    max_pos = zeros(2, 1);     % Max position

    % Insert one shape per roi
    for i = 1:size(roi_store, 2)

        temp_displ = displ(i);  % Displacement

        % Save min and max displacement positions
        if (temp_displ == min_displ)
            min_pos = roi_store{i}(1:2);

        elseif (temp_displ == max_displ)
            max_pos = roi_store{i}(1:2);

        end


        % Find Color
        if (isnan(temp_displ))
            temp_color = [0, 0, 0];  % If displacement is NaN, then make the square black

        else
            hue = (temp_displ - min_displ) / (max_displ - min_displ) * (min_color - max_color) + max_color; 
            hsv_color = [hue, 1.0, 1.0];              % HSV color
            temp_color = hsv2rgb(hsv_color) * 256;    % Convert from hsv to rgb
        end

        % Insert Shape
        output = insertShape(output, 'FilledRectangle', roi_store{i}, 'Color', temp_color);

    end

    % Make sure labels won't go off image
    offset_label = 50;      % Offset amount for label

    % Check x positions
    if (min_pos(1) + 10*offset_label > size(img, 2))
        min_pos(1) = min_pos(1) - 10 * offset_label;
    else
        min_pos(1) = min_pos(1) + offset_label;
    end

    if (max_pos(1) + 10*offset_label > size(img, 2))
        max_pos(1) = max_pos(1) - 10 * offset_label;
    else
        max_pos(1) = max_pos(1) + offset_label;
    end

    % Check y positions
    if (min_pos(2) + 10*offset_label > size(img, 1))
        min_pos(2) = min_pos(2) - 10 * offset_label;
    else
        min_pos(2) = min_pos(2) + offset_label;
    end

    if (max_pos(2) + 10*offset_label > size(img, 1))
        max_pos(2) = max_pos(2) - 10 * offset_label;
    else
        max_pos(2) = max_pos(2) + offset_label;
    end

    % Insert min and max labels
    output = insertText(output, min_pos, min_str, 'FontSize',150,'BoxColor','white');
    output = insertText(output, max_pos, max_str, 'FontSize',150,'BoxColor','white');

    % Plot Figure
    figure
    imshow(output)
    title(title_str)

    % Make color gradient for scale
    colors = [];                % Inialize color matrix
    for i = max_color:-0.001:min_color
        hsv_color = [i, 1.0, 1.0];            % HSV color
        rgb_color = hsv2rgb(hsv_color);        % RGB color

        colors = vertcat(colors, rgb_color);    % Add to color gradient
    end

    % Color scale
    colormap(colors);
    c = colorbar;

    c.Ticks = linspace(0, 1, 10);
    tick_labels = string(round(linspace(min_displ, max_displ, 10), 2));   % Tick mark labels
    c.TickLabels = tick_labels;
    c.FontSize = 15;

end