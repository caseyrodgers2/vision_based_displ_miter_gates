% Miter Gate: Weld Access Hole Plot Animations
% Note: This script is meant to be post processing for the
% motion_compensate_homography_klt_video.m script
% @author: Casey Rodgers, Shuo Wang
% @ date:  07/07/2021
% version : V1

clear all;
close all;
clc;

%% Load and Input Point Numbers of Interest
% Load the disp_data file
load('disp_data_images.mat');

displType = 1; % Point displacements or ROI displacements? 0 = point. 1 = roi.


% Input point numbers of interest points
point_ind = [443, 448, ... 
    454, 453, 460, 459, ...
    461, 464, 469, ...
    475, 473, 480, 479, ...
    484, 483, 490, 486, ...
    493, 491, 499, ...
    502, 503, 510, 507, ...
    515, ...
    512, ...
    520, ...
    518, ...
    522,  ...
    530, ...
    535, ...
    538, ...
    543, ...
    549, ...
    553, ...
    558];    % Point numbers (aka indices) of interest points


% Gate Values
centerTopPoint = [3736, 428]; % Top point coords of center between two gates (will be used to distinguish points on left gate from right gate)
centerBotPoint = [3435, 3778];  % Bot point coords of center between two gates 

topLeftGateX = 3583;       % X coord for top of left gate in pixels
topLeftGateY = 405;       % Y coord for top of left gate in pixels

topRightGateX = 3833;      % Y coord for top of right gate in pixels
topRightGateY = 471;      % Y coord for top of right gate in pixels


% origin at the bottom left corner of the left leaf
e1_imageplane = [0.643,0,-0.767];%image plane x axis unit vector in global coordinates
e2_imageplane = [-0.324,0.906,-0.272];%image plane y axis unit vector in global coordinates
e3_imageplane = [0.694,0.422,0.583];%image plane z axis unit vector in global coordinates  %camera axis

%if we can get camera location and image center 3D location
% a = [ax,ay,az]; %3D coordinates of the image center
% b = [ax,ay,az]; %3D coordinates of the camera
% v = b-a;
% v0 = 1/norm(v).*v;

% e1_leftleaf = [0.861,0,-0.509];%left leaf x axis unit vector in global coordinates
% e2_leftleaf = [0,1,0];%left leaf y axis unit vector in global coordinates
% e3_leftleaf = [0.509,0,0.861];%left leaf z axis unit vector in global coordinates
% 
% e1_rightleaf = [0.861,0,0.509];%right leaf x axis unit vector in global coordinates
% e2_rightleaf = [0,1,0];%right leaf y axis unit vector in global coordinates
% e3_rightleaf = [-0.509,0,0.861];%right leaf z axis unit vector in global coordinates


% Video Properties
frameRate = 2;          % Frame Rate for plot animations


%% Extract Displacements and Coordinates into Smaller Arrays

% Point Displacements
if (displType == 0)

    % Initialize arrays
    point_coords = zeros(length(point_ind), 2);     % Interested point coords
    %point_displ = zeros(length(point_ind), size(store_disp_mm, 1));      % Interested point displ
    point_displx = zeros(length(point_ind), size(store_dispx2_mm, 1));      % Interested point displ x
    point_disply = zeros(length(point_ind), size(store_dispy2_mm, 1));      % Interested point displ y
    point_scales2 = zeros(length(point_ind), 1);      % Interested point scales

    % Loop and extract values
    for i = 1:length(point_ind)
        index = point_ind(i);                       % Point index
        point_coords(i, :) = points_roi(index, :);    % Extract coords
        %point_displ(i, :) = store_disp_mm(:, index)';     % Extract displ 
        point_displx(i, :) = store_dispx2_mm(:, index)';    % Extract displ x
        point_disply(i, :) = store_dispy2_mm(:, index)';    % Extract displ y
        point_scales2(i) = point_scales(index);             % Extract scales
    end

    % Save arrays
    save('disp_of_interest', 'point_coords', 'point_displ', 'point_displx', 'point_disply', 'point_scales2');

    
% ROI Displacements
else
    point_coords = roi_centroid(num_static_roi+1 : num_static_roi + num_nonstatic_roi, :);

    point_displx = store_avgX_pix .* roi_scales;
    point_displx = point_displx(:, num_static_roi+1 : num_static_roi + num_nonstatic_roi)';

    point_disply = store_avgY_pix .* roi_scales;
    point_disply = point_disply(:, num_static_roi+1 : num_static_roi + num_nonstatic_roi)';
    
    point_scales2 = roi_scales;

end



%% Convert Displacements from Image Coords to Global Coords
displ_globalx = zeros(size(point_displx));   % Empty global x displ
displ_globaly = zeros(size(point_displx));   % Empty global y displ
displ_globalz = zeros(size(point_displx));   % Empty global z displ

% Loop through and calculate
for i = 1:size(point_displx, 1)
    for j = 1:size(point_displx, 2)
        dispx = point_displx(i,j); 
        dispy = point_disply(i,j);

        %project to global
        project_global = dispx.*e1_imageplane + dispy.*e2_imageplane;
        displ_globalx(i,j) = project_global(1);
        displ_globaly(i,j) = project_global(2);
        displ_globalz(i,j) = project_global(3);

%         %project to left leaf
%         A = [e1_leftleaf; e2_leftleaf; e3_leftleaf];
%         project_leftleaf = project_global*inv(A);
% 
%         %project to right leaf
%         B = [e1_rightleaf; e2_rightleaf; e3_rightleaf];
%         project_rightleaf = project_global*inv(B);

    end
end



%% Separate Points by Left and Right Gate
% Find gate center line equation
centerSlope = (centerBotPoint(1) - centerTopPoint(1)) / (centerBotPoint(2) - centerTopPoint(2)); % Slope
centerB = centerTopPoint(1) - centerSlope * centerTopPoint(2);  % Y intercept of line eq
findCenterX = @(yCoord) centerSlope * yCoord + centerB;        % Symbolic expression for finding x coord of center at a certain y coord

% Left Gate
leftNum = 0;                % Number of left points
left_distToTop = [];        % Distance to top
left_displx = [];            % Displ x
left_disply = [];            % Displ y
left_displz = [];            % Displ z
left_scale = [];            % Scale

% Right Gate
rightNum = 0;                % Number of right points
right_distToTop = [];        % Distance to top
right_displx = [];            % Displ x
right_disply = [];            % Displ y
right_displz = [];            % Displ z
right_scale = [];            % Scale

% Loop through points and sort
for i = 1:size(point_displx, 1)
    
    centerXcoord = findCenterX(point_coords(i,2)); % Find x coord of center of gate
    
    % Left gate
    if (point_coords(i,1) < centerXcoord)
        leftNum = leftNum + 1;
        left_distToTop(leftNum) = -sqrt((topLeftGateX - point_coords(i,1)).^2 + (topLeftGateY - point_coords(i,2)).^2);
        left_displx(leftNum,:) = displ_globalx(i,:);
        left_disply(leftNum,:) = displ_globaly(i,:);
        left_displz(leftNum,:) = displ_globalz(i,:);
        left_scale(leftNum,:) = point_scales2(i);
        
    % Right gate
    else
        rightNum = rightNum + 1;
        right_distToTop(rightNum) = -sqrt((topRightGateX - point_coords(i,1)).^2 + (topRightGateY - point_coords(i,2)).^2);
        right_displx(rightNum,:) = displ_globalx(i,:);
        right_disply(rightNum,:) = displ_globaly(i,:);
        right_displz(rightNum,:) = displ_globalz(i,:);
        right_scale(rightNum,:) = point_scales2(i);
    end
        
end

% Convert distance to top of gate from pix to mm to cm
left_distToTop = left_distToTop .* left_scale' ./ 10;
right_distToTop = right_distToTop .* right_scale' ./ 10;


%% Animated Plots

%% X Plot
% Find min and max displ for x axis
minLeftx = min(min(displ_globalx));
maxLeftx = max(max(displ_globalx));

% Video Params
video_name_x = 'animation_x.mp4';
v_x = VideoWriter(video_name_x, 'MPEG-4');
v_x.Quality = 99;
v_x.FrameRate = frameRate;
open(v_x)

% Figure
figure
for i = 1:size(point_displx, 2)
    
    % Plot
    if (i == 1 || i == size(point_displx, 2))
        plot(left_displx(:,i), left_distToTop, '-o', 'Color', '#0002C1')
        hold on
        plot(right_displx(:,i), right_distToTop, '-s', 'Color', '#C20000')
        hold on
    else
        plot(left_displx(:,i), left_distToTop, '-o', 'Color', '#C5C6FF')
        hold on
        plot(right_displx(:,i), right_distToTop, '-s', 'Color', '#FFBFBF')
        hold on
    end

%     plot(left_displx(:,i), left_distToTop, '-o')
%     xlim([minLeftx maxLeftx])
%     hold on
%     plot(right_displx(:,i), right_distToTop, '-s')
%     hold off

    % Axis
    xlim([minLeftx maxLeftx])
    
    % Labels and grid
    title('Weld Access Holes X Movement in Global Coords')
    xlabel('X Displacement (mm)')
    ylabel('Gate Height (cm) (top = 0)')
    legend('Left', 'Right', 'Location', 'eastoutside')
    grid on
    
    % Draw and pause
    drawnow
    pause(1/frameRate)
    
    % Save to video
    frame = getframe(gcf); %get frame
    writeVideo(v_x, frame);
    
end

hold off
close(v_x)  % Close video



%% Z Plot
% Find min and max displ for z axis
minLeftz = min(min(displ_globalz));
maxLeftz = max(max(displ_globalz));

% Video Params
video_name_z = 'animation_z.mp4';
v_z = VideoWriter(video_name_z, 'MPEG-4');
v_z.Quality = 99;
v_z.FrameRate = frameRate;
open(v_z)

% Figure
figure
for i = 1:size(point_displx, 2)
    
    % Plot
    if (i == 1 || i == size(point_displx, 2))
        plot(left_displz(:,i), left_distToTop, '-o', 'Color', '#0002C1')
        hold on
        plot(right_displz(:,i), right_distToTop, '-s', 'Color', '#C20000')
        hold on
    else
        plot(left_displz(:,i), left_distToTop, '-o', 'Color', '#C5C6FF')
        hold on
        plot(right_displz(:,i), right_distToTop, '-s', 'Color', '#FFBFBF')
        hold on
    end
    
%     plot(left_displz(:,i), left_distToTop, '-o')
%     xlim([minLeftz maxLeftz])
%     hold on
%     plot(right_displz(:,i), right_distToTop, '-s')
%     hold off

    % Axis
    xlim([minLeftz maxLeftz])
    
    % Labels and grid
    title('Weld Access Holes Z Movement in Global Coords')
    xlabel('Z Displacement (mm)')
    ylabel('Gate Height (cm) (top = 0)')
    legend('Left', 'Right', 'Location', 'eastoutside')
    grid on
    
    % Draw and pause
    drawnow
    pause(1/frameRate)
    
    % Save to video
    frame = getframe(gcf); %get frame
    writeVideo(v_z, frame);
    
end

hold off
close(v_z)  % Close video



%% Y Plot - Displacement plotted along gate height
% Find min and max displ for y axis
distToTop_all = [left_distToTop, right_distToTop];
biggesty = max(max(abs(displ_globaly)));
minLefty = min(distToTop_all) - biggesty - 100;
maxLefty = max(distToTop_all) + biggesty + 100;

% Video Params
video_name_y = 'animation_y.mp4';
v_y = VideoWriter(video_name_y, 'MPEG-4');
v_y.Quality = 99;
v_y.FrameRate = frameRate;
open(v_y)

% Figure
figure
for i = 1:size(point_displx, 2)
    
    % Plot
    if (i == 1 || i == size(point_displx, 2))
        plot(zeros(length(left_distToTop)), left_distToTop' + left_disply(:,i), '-o', 'Color', '#0002C1')
        hold on
        plot(ones(length(right_distToTop)), right_distToTop' + right_disply(:,i), '-s', 'Color', '#C20000')
        hold on
    else
        plot(zeros(length(left_distToTop)), left_distToTop' + left_disply(:,i), '-o', 'Color', '#C5C6FF')
        hold on
        plot(ones(length(right_distToTop)), right_distToTop' + right_disply(:,i), '-s', 'Color', '#FFBFBF')
        hold on
    end
    
%     plot(zeros(length(left_distToTop)), left_distToTop + left_disply(:,i), '-o')
%     xlim([-1 2])
%     ylim([minLefty maxLefty])
%     hold on
%     plot(ones(length(right_distToTop)), right_distToTop + right_disply(:,i), '-s')
%     hold off

    % Axis
    ylim([minLefty maxLefty])
    
    % Labels and grid
    title('Weld Access Holes Y Movement in Global Coords')
    xlabel('Left Gate at 0. Right Gate at 1.')
    ylabel('Gate Height (cm) (top = 0)')
    %legend('Left', 'Right', 'Location', 'eastoutside')
    grid on
    
    % Draw and pause
    drawnow
    pause(1/frameRate)
    
    % Save to video
    frame = getframe(gcf); %get frame
    writeVideo(v_y, frame);
    
end

hold off
close(v_y)  % Close video



%% Y Plot 2 - Displacement plotted on x axis so it's easier to see
% Find min and max displ for y axis
minLefty2 = min(min(displ_globaly));
maxLefty2 = max(max(displ_globaly));

% Video Params
video_name_y2 = 'animation_y2.mp4';
v_y2 = VideoWriter(video_name_y2, 'MPEG-4');
v_y2.Quality = 99;
v_y2.FrameRate = frameRate;
open(v_y2)

% Figure
figure
for i = 1:size(point_displx, 2)
    
    % Plot
    if (i == 1 || i == size(point_displx, 2))
        plot(left_disply(:,i), left_distToTop, '-o', 'Color', '#0002C1')
        hold on
        plot(right_disply(:,i), right_distToTop, '-s', 'Color', '#C20000')
        hold on
    else
        plot(left_disply(:,i), left_distToTop, '-o', 'Color', '#C5C6FF')
        hold on
        plot(right_disply(:,i), right_distToTop, '-s', 'Color', '#FFBFBF')
        hold on
    end
    
%     plot(left_disply(:,i), left_distToTop, '-o')
%     xlim([minLefty2 maxLefty2])
%     hold on
%     plot(right_disply(:,i), right_distToTop, '-s')
%     hold off

    % Axis
    xlim([minLefty2 maxLefty2])
    
    % Labels and grid
    title('Weld Access Holes Y Movement in Global Coords')
    xlabel('Y Displacement (mm)')
    ylabel('Gate Height (cm) (top = 0)')
    legend('Left', 'Right', 'Location', 'eastoutside')
    grid on
    
    % Draw and pause
    drawnow
    pause(1/frameRate)
    
    % Save to video
    frame = getframe(gcf); %get frame
    writeVideo(v_y2, frame);
    
end

hold off
close(v_y2)  % Close video



%% All X, Y, & Z Plot
% We already found x, y axis, and z axis limits so no need  to find them again.

% Video Params
video_name_all = 'animation_all.mp4';
v_all = VideoWriter(video_name_all, 'MPEG-4');
v_all.Quality = 99;
v_all.FrameRate = frameRate;
open(v_all)

% Figure
figure
h1 = axes;
for i = 1:size(point_displx, 2)
    
    % Plot
    if (i == 1 || i == size(point_displx, 2))
        plot3(left_displz(:,i), left_displx(:,i), left_distToTop' + left_disply(:,i), '-o', 'Color', '#0002C1')
        hold on
        plot3(right_displz(:,i), right_displx(:,i), right_distToTop' + right_disply(:,i), '-s', 'Color', '#C20000')
        hold on
    else
        plot3(left_displz(:,i), left_displx(:,i), left_distToTop' + left_disply(:,i), '-o', 'Color', '#C5C6FF')
        hold on
        plot3(right_displz(:,i), right_displx(:,i), right_distToTop' + right_disply(:,i), '-s', 'Color', '#FFBFBF')
        hold on
    end
    
%     plot3(left_displz(:,i), left_displx(:,i), left_distToTop' + left_disply(:,i), '-o')
%     hold on
%     plot3(right_displz(:,i), right_displx(:,i), right_distToTop' + right_disply(:,i), '-s')
%     hold off
    
    % Axes
    xlim([minLeftz maxLeftz])
    ylim([minLeftx maxLeftx])
    zlim([minLefty maxLefty])
    set(h1, 'Xdir', 'reverse')
    set(h1, 'Ydir', 'reverse')
    
    % Labels and grid
    title('Weld Access Holes 3D Movement in Global Coords')
    ylabel('X Displ (mm)')
    zlabel('Gate Height (cm) (top = 0)')
    xlabel('Z Displ (mm)')
    legend('Left', 'Right', 'Location', 'eastoutside')
    grid on
    
    % Draw and pause
    drawnow
    pause(1/frameRate)
    
    % Save to video
    frame = getframe(gcf); %get frame
    writeVideo(v_all, frame);
    
end

hold off
close(v_all)  % Close video

