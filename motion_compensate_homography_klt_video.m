% Miter Gate Feature based Multi-ROI tracking
% @author: Casey Rodgers, Shaik Althaf V S, Shuo Wang
% @ date:  12/31/2021

clc
clear all
close all

%% Load Images & Set up Result Folder

% Get images
imgFolder = 'blender_5x\';   % Image folder (Folder where the images are)
filePattern = fullfile(imgFolder, '*.jpg');
imgFiles = dir(filePattern);

% Make sure result folder exists / create it
imgFolder_name = imgFolder(1:length(imgFolder)-1);  % Folder name
username = getenv('username');                      % Username
curr_time = datestr(datetime('now'));               % Date and time
result_folder = strcat(imgFolder_name, '_', username, '_', curr_time, '_results\');
result_filePath = strcat(imgFolder, result_folder);
result_filePath = result_filePath(~isspace(result_filePath));   % Get rid of spaces
result_filePath = strrep(result_filePath, '-', '_');            % Get rid of dashes
result_filePath = strrep(result_filePath, ':', '_');            % Get rid of colons

if ~exist(result_filePath, 'dir')
    mkdir(result_filePath)
end



%% User Input Choices

%%% Frames and Video
skipFrames = 0;     % Number of frames you want to skip at the beginning
lastFrame = size(imgFiles,1);  % Last frame for video. Must be number relative to first frame in folder. (ie in 200 image set, lastFrame = 200 for all images)
every_frame_v = 20;    % Want a video of every __ (ie 5) frames. 1 means all frames.
every_fig = 5;      % Want a matlab figure every __ (ie 5) frames of video. 1 means all figures.
onlyLast = 1;       % Only want to do first and last image? 0 = no. 1 = yes
whatFrameR = 5;    % What frame rate do you want? # = number of frames / second

%%% Images 
rotation = 0;     % How much do the images need to be rotated?
undistort = 1;     % Undistort the image? 0 = No. 1 = yes.

%%% ROI's and ROI Points
useLastStaticRoi = 1; % Want to use previously saved static roi? 0 = no. 1 = yes
useLastNonStaticRoi = 0; % Want to use previously saved nonstatic roi? 0 = no. 1 = yes
num_nonstatic_roi = 30 ;  % ENTER number of Nonstatic ROI (Region of Interest)
num_static_roi = 4; % Enter number of static roi
threshold = 0.5; % MaxDistance for calculating homography
homographyinliers = 10; %minimum number of inliers for calculating homography

%%% Use node file for points? (0 = No. 1 = Yes.)
camCoords = [0.279991, -12.1451, 36.0736];      % Camera Global coords
blenderAdj = [27.6568, 0, 0];                   % Adjustement to match Blender movement

% Nonstatic
node_file_nonstatic = 1;                           % For nonstatic points
nodeFile_nonstatic = 'nodePixel.csv';              % Node location file ((# points, 2) size array with pixel locations)
node_nonstatic_roi_halfSize = 15;                   % Half the size of the roi you want around the node
node_coord_file_nonstatic = 'nodeCoordinates.csv'; % Global 3D Nodal coordinates (in INCHES)

% Static
node_file_static = 0;                        % For static points
nodeFile_static = '';                        % Node location file
node_static_roi_halfSize = 15;                % Half the size of the roi you want around the node
node_coord_file_static = '';                 % Global 3D Nodal coordinates (in INCHES)


%%% Choose Manual or Automatic Point Selection
points_select_type = 0; % Choose Manual or Automatic Point Selection for each ROI? 0 = no. 1 = yes.
% If no, then this will be used for static / nonstatic regions:
points_select_static = 1;     % Point Selection type. 0 = Manual. 1 = Automatic
points_select_nonstatic = 1;  % Point Selection type. 0 = Manual. 1 = Automatic.

%%% Choose Feature Selection Type
chooseFeatType = 0;    % Choose feature selection type for each ROI? 0 = no. 1 = yes.
% If no, then this will be used for static / nonstatic regions:
featType_static = 1;     % Feature Selection type. 0 = Strongest. 1 = Uniform
featType_nonstatic = 0;  % Feature Selection type. 0 = Strongest. 1 = Uniform

%%% Choosing Number of Features
chooseNumFeat = 0; % Choose number of features for each ROI? 0 = no. 1 = yes.
% If no, then this will be used for static / nonstatic regions:
Num_features_static = 120 ; % Number of points to use select in per ROI
Num_features_nonStatic = 5 ; % Number of points to use select in per ROI

%%% KLT Tracking
vanish_thres = 250;  % Ignore feature point if motion is greater than __ mm. 
fill_emp = 0;       % Will water block the features? If yes, is it a fill or empty? 0 = fill. 1 = empty. Set to 0 if no water involved.
averaging = 6;      % How many frames on each side do you want to average for the centered moving average? 0 = no averaging.
hist_eq = 1;        % Do you want histogram eqivalization to normalize global lighting/sunlight/shadow effect? 0 = no. 1 = yes.
detector_option = 4 ; % FEature Selection Method - 1. detectMinEigenFeatures, 2. BRSIK, 3. SURF, 4. Harris Corner ; % Recommended - OPTION 4. ( Harrs Corner)
num_for_ref_avg = -1;  % Number of frames at the beginning that should be averaged for the reference coords for tracking
BlockSize_static = 15;  %block size of KLT tracker for static points
BlockSize_nonstatic = 75;   %block size of KLT tracker for nonstatic points (15 for blender)
KLTMaxBidirectionalError = 10;  %Forward-backward error threshold of KLT

%%% Displacement Mean and Standard Deviation
numOutlierRounds = 1;   % Number of times you want to calculate the STD and then throw out values that are higher than 2 STD from mean
numStd = 1;             % Number of standard deviations points must be within to not be thrown out as outliers
stdThr = 0.5;           % Standard Deviation threshold. If the standard deviation for an roi is above this value, then it will get an Nan

%%% Plotting
q_scale  = 400;     % Quiver arrow scale
score_min = 0;    % What minimum score must a point have to be plotted? Range: 0-1. 0 = terrible tracking. 1 = best tracking

%%% Debugging
debug_static = 0; % Whether to check static feature points before looping through photos, 1 check, 0 not check
check_frame_interval = 300; %check every x frames for static/target feature points


%%% Scale factor Estimation Method 
% 1. Use Known dimension at ROI and distance to ROI 
% 2. Use Focal length from Metadat and Distance to ROI
scale_method = 1;

%%% ccompensation_method: (0) None, (1) ImgRegister, (2) Homography BEFORE
% KLT tracking, (3) Homography AFTER KLT tracking
compensation_method = 3;



%% SELECT MULTIPLE ROI and SCALE FACTORS & MOTION COMPENSATION - HOMOGRAPHY

% Figure out starting image (beginning for fill. end for empty.)
start = -1;
if (fill_emp == 0)
    start = skipFrames+1;
else
    start = size(imgFiles,1) - skipFrames;
end

% Initialize the ROI region on Frame 1
filepath = imgFiles(1).folder;
firstFrameName = split(imgFiles(start).name, ["_", "."]);   % Name of first frame with ext
firstFrameNum = str2double(firstFrameName(2));              % Get number of first frame
img0_name = fullfile(filepath, imgFiles(start).name);
img0 = imread(char(strcat(img0_name)));


% Load nodal coordinates (if necessary)
if (node_file_nonstatic == 1)
    nodePixel_nonstatic = load(nodeFile_nonstatic);                 % Load file for pixel locations
    nodeGlobalCoords_nonstatic = load(node_coord_file_nonstatic);    % Load file for 3D global locations for scales
    nodeGlobalCoords_nonstatic = nodeGlobalCoords_nonstatic * 0.0254; % Convert from INCHES to meters
end

if (node_file_static == 1)
    nodePixel_static = load(nodeFile_static);                  % Load file for pixel locations)
    nodeGlobalCoords_static = load(node_coord_file_static);    % Load file for 3D global locations for scales
    nodeGlobalCoords_static = nodeGlobalCoords_static * 0.0254; % Convert from INCHES to meters
end


% Undistort image
if (undistort == 1)
    load('calibration_corrected.mat');
    params = cameraParams8;
    [undistort_img0, newOrigin] = undistortImage(img0, params, 'OutputView', 'full');
    objectFrame = undistort_img0;
else
    objectFrame = img0;
end


% Extract focal length f_pix
info = imfinfo(char(strcat(img0_name)));
% extarct param
%f_mm = info.DigitalCamera.FocalLength;
f_mm = 24;
% ######### FOR NIKON _Z-5 ##############  Change as per Cmaera Model
sensor_width_mm =  35.9;
mm_to_pix = 6016/sensor_width_mm ;
pix_res = 1/mm_to_pix; %mm
f_pix = f_mm * mm_to_pix;



%% Rotate (if necessary)
objectFrame = imrotate(objectFrame,rotation);
A_ref = objectFrame;
total_static_roi = 0;       % Total static roi for plotting later

% Static Regions
if (useLastStaticRoi ~= 1)
    
    [static_roi_store, static_roi_scales, num_static_roi, good_static_rois] = GetROI(node_file_static, node_static_roi_halfSize, A_ref, scale_method, f_pix, nodePixel_static, nodeGlobalCoords_static, camCoords, blenderAdj);

    save(fullfile(imgFolder_name,'static_roi_store2'), 'static_roi_store', 'static_roi_scales', 'StaticObjectFrame', 'num_static_roi', 'good_static_rois');
    
else 
    good_static_rois = [1, 1, 1, 1];
    load(fullfile(imgFolder_name, 'static_roi_store2'))
end


% Non Static Regions
if (useLastNonStaticRoi ~= 1)
    
    [nonstatic_roi_store, nonstatic_roi_scales, num_nonstatic_roi, good_nonstatic_rois] = GetROI(node_file_nonstatic, node_nonstatic_roi_halfSize, A_ref, scale_method, f_pix, nodePixel_nonstatic, nodeGlobalCoords_nonstatic, camCoords, blenderAdj);
        
    save(fullfile(imgFolder_name,'nonstatic_roi_store3'), 'nonstatic_roi_store', 'nonstatic_roi_scales', 'num_nonstatic_roi', 'good_nonstatic_rois');
    
else
    load(fullfile(imgFolder_name, 'nonstatic_roi_store3'))
end
       

% Combine Static and NonStatic ROI
roi_store = [static_roi_store, nonstatic_roi_store];   % ROI store
roi_scales = [static_roi_scales, nonstatic_roi_scales]; % ROI scales
good_rois = [good_static_rois, good_nonstatic_rois];    % Good arrays


% % Adjust scales to be correct
% static_adj = ones(1, num_static_roi);
% nonstatic_adj_left = horzcat(linspace(1, 7.0866/8.83, 9), linspace(9.0678/11.92, 10.34796/14.72, 4), repmat(16/20.34, 1, 5));
% nonstatic_adj_right = horzcat(linspace(1, 7.0866/8.94, 9), linspace(9.0678/11.92, 12.2172984/18.16 , 8), [16/21.6]);
% roi_scale_adj = horzcat(static_adj, zeros(1, num_nonstatic_roi));
% stat_adj_count = 1;
% nonstat_adj_count = 1;
% for i = num_static_roi+1 : num_static_roi+num_nonstatic_roi
%     if mod(i, 2) == 1
%         roi_scale_adj(i) = nonstatic_adj_left(stat_adj_count);
%         stat_adj_count = stat_adj_count + 1;
%     else
%         roi_scale_adj(i) = nonstatic_adj_right(nonstat_adj_count);
%         nonstat_adj_count = nonstat_adj_count + 1;
%     end
% end
% 
% roi_scales = roi_scales .* roi_scale_adj;

% Light Normalization 
if hist_eq == 1   %ROIs change color
    
    resImg_A = lightNormal(A_ref, roi_store, BlockSize_nonstatic, BlockSize_static);
    A_ref = resImg_A;

    figure
    imshow(A_ref);

end



%% Detect interest points in the MULTIPLE object region.

point_scales = [];   % pix2mm for each point
point_num = [];      % Number of Points per ROI
roi_centroid = size(length(roi_store), 2);   % Centroid x, y for each ROI
%figure,

pointImage = A_ref;
points_roi = [];

for  i= 1 : num_nonstatic_roi+num_static_roi
    objectRegion = roi_store{i};
    
    % Find roi centroid
    roi_centroid(i, 1) = roi_store{i}(1) + roi_store{i}(3) ./ 2;
    roi_centroid(i, 2) = roi_store{i}(2) + roi_store{i}(4) ./ 2;

    % Feature selection method
    % Number of features to be used for ROI
    if points_select_type == 1
        if i <= num_static_roi
            prompt = 'Choose point selection method (0 manual, 1 automatic) for Static ROI #' + string(i) + ': ';
        else
            prompt = 'Choose point selection method (0 manual, 1 automatic) for NonStatic ROI #' + string(i-num_static_roi) + ': ';
        end
        point_selection_method = input(prompt);
    else
        if i <= num_static_roi
            point_selection_method = points_select_static;
        else
            point_selection_method = points_select_nonstatic;
        end
    end
    
    % No averaging for klt tracking
    if num_for_ref_avg == -1
        % Manually selecting features to track
        if point_selection_method == 0
            prompt = 'How many points do you want to manually select?: ';
            numPoints = input(prompt);
            
            pointImage2 = insertShape(pointImage,'Rectangle',objectRegion,'Color','red');
            imshow(pointImage2);
            impixelinfo;
            
            points_avg = zeros(numPoints, 2);
            
            for k = 1:numPoints
                point_roi = drawpoint('DrawingArea', objectRegion);
                points_avg(k, :) = point_roi.Position;
            end
            
            point_num = [point_num, numPoints];

%           points = [];
%           [x,y] = getpts;
%           points = vertcat(points, [x,y]);
%           disp('Choose the points you want in the left image only. You may need to zoom and pan to the desired ROI. Then File > Close the cpselect tool when done selecting points.')
%           [selectedMovingPoints,selectedFixedPoints] = cpselect(A_ref,A_ref, 'Wait', true);
%           points = selectedMovingPoints;
%           points_avg = points;

            
        % Using a feature detector to automatically detect features
        else
            points = Feature_detector(A_ref,objectRegion,detector_option);

            % Number of features to be used for ROI
            if chooseNumFeat == 1
                if i <= num_static_roi
                    prompt = 'Choose number of features for Static ROI #' + string(i) + ': ';
                else
                    prompt = 'Choose number of features for NonStatic ROI #' + string(i-num_static_roi) + ': ';
                end
                Num_features = input(prompt);
            else
                if i <= num_static_roi
                    Num_features = Num_features_static;
                else
                    Num_features = Num_features_nonStatic;
                end
            end


            % Feature Selection Type to be used for ROI
            % If you are choosing the feature types manually for each ROI
            if chooseFeatType == 1
                if i <= num_static_roi
                    prompt = 'Choose feature selection type (0 = strongest. 1 = uniform.) for Static ROI #' + string(i) + ': ';
                else
                    prompt = 'Choose feature selection type (0 = strongest. 1 = uniform.) for NonStatic ROI #' + string(i-num_static_roi) + ': ';
                end
                feat_select_type = input(prompt);
            % Choosing the feature types the same way for each type of ROI
            else
                if i <= num_static_roi
                    feat_select_type = featType_static;
                else
                    feat_select_type = featType_nonstatic;
                end
            end
            
            % Select the feature points by either: strongest or uniform
            if feat_select_type == 0
                strongPoints = points.selectStrongest(Num_features);
            else
                strongPoints = selectUniform(points, Num_features, size(A_ref));
            end
            
            points_avg = strongPoints.Location;

        end

    % Yes, averaging for klt tracking
    else
        for j=1:num_for_ref_avg-1

            % Prevent going out of index
            avg_index = 1;
            % Fill / No water
            if fill_emp == 0
                avg_index = start + j;
                if (avg_index > size(imgFiles,1))
                    ref_divide = j;
                    break
                else
                    ref_divide = num_for_ref_avg;
                end

            % Empty
            elseif fill_emp == 1
                avg_index = start - j;
                if (avg_index > size(imgFiles,1))
                    ref_divide = j;
                    break
                else
                    ref_divide = num_for_ref_avg;
                end
            end

            kltImg_name = fullfile(filepath, imgFiles(avg_index).name);
            kltImg = imread(char(strcat(kltImg_name)));
            pointImage = kltImg;
            pointImage = imrotate(pointImage,rotation);

            objectImage = insertShape(pointImage,'Rectangle',objectRegion,'Color','red');
            imshow(objectImage);

            % Good features to track
            if point_selection_method == 1
                points = [];
                % [x,y] = getpts;
                % points = vertcat(points, [x,y]);
                imshow(objectImage);
                disp('Choose similar points in both left and right image. Then File > Close the cpselect tool when done selecting points.')
                [selectedMovingPoints,selectedFixedPoints] = cpselect(objectImage,objectImage, 'Wait', true);
                points = selectedMovingPoints;

                if j == 1
                    points_avg = points;
                else
                    points_avg = points_avg + points;
                end

            else
                %points = detectMinEigenFeatures(rgb2gray(A_ref),'ROI',objectRegion);
                first_points = Feature_detector(A_ref,objectRegion,detector_option);
                [first_features, first_points] = extractFeatures(rgb2gray(A_ref),first_points);

                motion_points = Feature_detector(pointImage,objectRegion,detector_option);
                [motion_features, motion_points] = extractFeatures(rgb2gray(pointImage),motion_points);

                indexPairs = matchFeatures(motion_features, first_features, 'Unique', true); 

                movingPoints = motion_points(indexPairs(:,1), :);
                fixedPoints = first_points(indexPairs(:,2), :);

                % Number of points to use select in per ROI
                if i <= num_static_roi
                    Num_features = Num_features_static;
                else
                    Num_features = Num_features_nonStatic;
                end

                if feat_select_type == 0
                    m_points = movingPoints.selectStrongest(Num_features);
                    f_points = fixedPoints.selectStrongest(Num_features);
                else
                    m_points = selectUniform(movingPoints, Num_features, size(A_ref));
                    f_points = selectUniform(fixedPoints, Num_features, size(A_ref));
                end

                if j == 1
                    points_avg = f_points.Location + m_points.Location;
                else
                    points_avg = points_avg + m_points.Location;
                end

            end % End if selection method

        end % End j for loop

        points_avg = points_avg ./ ref_divide;

    end 


    % Are points for static roi? Add to total number count (for
    % plotting later)
    if i <= num_static_roi
        total_static_roi = total_static_roi + size(points_avg, 1);
    end

    % Add to points_roi vector
    if i == 1
        points_roi = points_avg;
        point_num = [point_num, size(points_avg, 1)];
    else
        points_roi = vertcat( points_roi , points_avg) ;
        point_num = [point_num, size(points_avg, 1)];
    end

    % Add to scale factors
    for j=1:size(points_avg, 1)
        point_scales = [point_scales, roi_scales(i)];
    end

    %Display
    pointImage = insertMarker(pointImage,points_roi,'+','Color','red','size',100);
    imshow(pointImage);

    title('Detected interest points');
    
    disp(' ')

end 

plot_name = fullfile(result_filePath, strcat(imgFolder_name, '_detectedPoints', '.png'));
imwrite(pointImage, plot_name);    % Save image as .png
impixelinfo;



%% Debugging
if debug_static == 1
    disp('Debugging...')
    last = floor((size(imgFiles,1)-start)/check_frame_interval)*check_frame_interval+start;
    %disp(start)
    %disp(check_frame_interval)
    %disp(last)

    % Read reference image and potentially rotate
    disp(' ')
    disp(['fixed image',imgFiles(start).name])
    fixedimage = imread(char(strcat(fullfile(filepath, imgFiles(start).name))));
    fixedimage = imrotate(fixedimage,rotation);

    % Undistort reference image
    if (undistort == 1)
            [fixedimage, newOrigin_a] = undistortImage(fixedimage, params, 'OutputView', 'full');
    end


%         for  j=1:num_static_roi
%             StaticObjectRegion = cell2mat(static_roi_store(:,j));
%             static_points = Feature_detector(fixedimage,StaticObjectRegion,detector_option);
%             [static_features, static_points] = extractFeatures(rgb2gray(fixedimage),static_points);
%             static_pointsReference = static_points;
%             static_featuresReference = static_features;
%             static_points = Feature_detector(movingimage,StaticObjectRegion,detector_option);
%             [static_features, static_points] = extractFeatures(rgb2gray(movingimage),static_points);
%             indexPairs = matchFeatures(static_features, static_featuresReference, 'Unique', true); 
%             
%             if j == 1
%                 movingPoints = static_points.Location(indexPairs(:,1), :);
%                 fixedPoints = static_pointsReference.Location(indexPairs(:,2), :);
%             else
%                 movingPoints = vertcat(movingPoints , static_points.Location(indexPairs(:,1), :));
%                 fixedPoints = vertcat(fixedPoints , static_pointsReference.Location(indexPairs(:,2), :));
%             end
%             %disp(size(movingPoints));
%         end
%         [tform,inlierIdx,status] = estimateGeometricTransform2D(movingPoints,fixedPoints,'rigid','Confidence', 99.9, 'MaxNumTrials', 2000,'MaxDistance',threshold);
%         disp(['status: ',num2str(status)])
%         disp(['number of inliers: ',num2str(sum(inlierIdx))])
%         
%         % If there's not enough matched points, then not much can
%         % be done about it. Continue to next loop
%         if (status == 1)
%             disp('matchedPoints1 and matchedPoints2 inputs do not contain enough points. Moving on to next frame...')
%             continue;
%         end
% 
%         % If there's not enough inliers, then try to fix it
%         if (sum(inlierIdx) < homographyinliers)
%             disp('Not enough inliers found. Attempting to fix...')
%         end
%         turn = 1;
%         while (sum(inlierIdx) < homographyinliers)
%             if (turn == 10)
%                 break
%             else
%                 error_message = strcat('Not enough inliers found. Attempt #', string(turn));
%                 disp(error_message)
%                 Threshold = threshold*(2^turn);
%                 [tform,inlierIdx,status] = estimateGeometricTransform2D(movingPoints, fixedPoints,...
%             'rigid', 'Confidence', 99.9, 'MaxNumTrials', 2000,'MaxDistance',Threshold);
%                 disp(['status: ',num2str(status)])
%                 disp(['number of inliers: ',num2str(sum(inlierIdx))])
%                 turn = turn+1;
%             end
%             if (onlyLast == 1)
%                 break
%             end
%         end
% 
%         if (sum(inlierIdx) < homographyinliers || status ~= 0)
%             disp('Not enough inliers could be found even after increasing threshold. Moving to next frame...')
%             continue;
%         end
%         
%         inlierPtsDistorted = movingPoints(inlierIdx,:);
%         inlierPtsOriginal  = fixedPoints(inlierIdx,:);
%         figure 
%         showMatchedFeatures(rgb2gray(fixedimage),rgb2gray(movingimage),inlierPtsOriginal,inlierPtsDistorted)
%         title(['Feature matching based Homography, Inlier Points of image ',imgFiles(start).name, ' and image ',imgFiles(i).name])

    % Get static and nonstatic points
    fixedPoints = [];
    movingPoints = [];
    for i = start+check_frame_interval: check_frame_interval: last
        disp(' ')
        disp(['moving image',imgFiles(i).name])

        % Get next image
        movingimage = imread(char(strcat(fullfile(filepath, imgFiles(i).name))));
        movingimage = imrotate(movingimage,rotation);
        if (undistort == 1)
            [movingimage, newOrigin_a] = undistortImage(movingimage, params, 'OutputView', 'full');
        end
        
        % Track the points
        tracker_static = vision.PointTracker('MaxBidirectionalError',KLTMaxBidirectionalError,...
                                       'BlockSize',[BlockSize_static,BlockSize_static],'NumPyramidLevels',3,...
                                       'MaxIterations',30);
        tracker_nonstatic = vision.PointTracker('MaxBidirectionalError',KLTMaxBidirectionalError,...
                                       'BlockSize',[BlockSize_nonstatic,BlockSize_nonstatic],'NumPyramidLevels',3,...
                                       'MaxIterations',30);
                                   
                                   
        % First / refernce Frame Feature Points 
        points_in_static=double((points_roi(1:total_static_roi,:)));         
        objectFrame = fixedimage;

        % Initialize the tracker.
        initialize(tracker_static,points_in_static,objectFrame);

        % Track with KLT on Frame B
        [points2_static,validity_static,scores_static] = tracker_static(movingimage);
        points_in_nonstatic=double((points_roi(total_static_roi+1:end,:)));         
        objectFrame = fixedimage;

        % Initialize the tracker.
        initialize(tracker_nonstatic,points_in_nonstatic,objectFrame);

        % Track with KLT on Frame B
        [points2_nonstatic,validity_nonstatic,scores_nonstatic] = tracker_nonstatic(movingimage);
        points_in=double((points_roi));
        points2=[points2_static;points2_nonstatic];
        validity=[validity_static;validity_nonstatic];
        scores=[scores_static;scores_nonstatic];

        % Display with Imfuse overlay and distance
        %  Quiver plot between 2 frames       
        pts1 = double(points_in) .* validity;
        pts2 = double(points2) .* validity;

        pts1 (all(pts1 == 0, 2),:) = [];%delete points that are invalid
        pts2 (all(pts2 == 0, 2),:) = [];%delete points that are invalid

        figure 
        showMatchedFeatures(rgb2gray(fixedimage),rgb2gray(movingimage),pts1,pts2)
        title(['KLT tracking points of image ',imgFiles(start).name, ' and image ',imgFiles(i).name])
        

%         H_mat = (tform.T)';
%         pts2_org = vertcat(double(points2)',ones(1,size(points2,1)));
%         pts2_H  = H_mat * pts2_org;
%         pts2_H_2d = pts2_H(1:2,:)';
%         pts2 = pts2_H_2d;
% 
%         figure 
%         showMatchedFeatures(rgb2gray(fixedimage),rgb2gray(movingimage),pts1,pts2)
%         title('Matched Points, KLT first, then homography')
%         
%         
%         movingimage_warp = imwarp(movingimage,tform,'OutputView',imref2d(size(rgb2gray(fixedimage))));
%         figure;
%         compare_warp = imfuse(fixedimage,movingimage_warp,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
%         imshow(compare_warp);
%         title('Warp fuse')
%         impixelinfo;

        % Prepare for camera compensation
        fixedPoints = pts1(1:sum(validity_static),:);
        movingPoints = pts2(1:sum(validity_static),:);
        [tform,inlierIdx,status] = estimateGeometricTransform2D(movingPoints,fixedPoints,'rigid','Confidence', 99.9, 'MaxNumTrials', 2000,'MaxDistance',threshold);
        
        inlierPtsDistorted = movingPoints(inlierIdx,:);
        inlierPtsOriginal  = fixedPoints(inlierIdx,:);

        % Plot
        figure 
        showMatchedFeatures(rgb2gray(fixedimage),rgb2gray(movingimage),inlierPtsOriginal,inlierPtsDistorted)
        title(['KLT tracking based Homography, Inlier Points of image ',imgFiles(start).name, ' and image ',imgFiles(i).name])
        
        % Camera compensation
        H_mat = (tform.T)';
        pts2_org = vertcat(double(pts2)',ones(1,size(pts2,1)));
        pts2_H  = H_mat * pts2_org;
        pts2_H_2d = pts2_H(1:2,:)';
        pts2 = pts2_H_2d;

        % Plot
        figure 
        showMatchedFeatures(rgb2gray(fixedimage),rgb2gray(movingimage),pts1,pts2)
        title(['KLT tracking based Homography, after warping, image ',imgFiles(start).name, ' and image ',imgFiles(i).name])
        
        movingimage_warp = imwarp(movingimage,tform,'OutputView',imref2d(size(rgb2gray(fixedimage))));
        figure;
        compare_warp = imfuse(fixedimage,movingimage_warp,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        imshow(compare_warp);
        title(['KLT tracking based Homography, warp fuse, image ',imgFiles(start).name, ' and image ',imgFiles(i).name])
        impixelinfo;
    end 
    
    disp(' ')
    disp('press enter to continue')
    pause;
    disp('Continuing...')
else
end



%% Loop through photos
% Initialize matrices
store_disp_x_pix = [];  % Store pixel displacement in x direction for all points
store_disp_y_pix = [];  % Store pixel displacement in y direction for all points
store_scores = [];      % Store tracking scores for all points
store_frameNum= [];     % Store frame number
store_img_index = [];   % Store image indices for images that are going to go in video
store_res_names = [];   % Store image names for quiver plots later
store_pts2 = [];        % Store pts2 for each frame for plotting later
store_pts1 = [];        % Store pts1 for plotting later
    
loop_num = 1;
water_cover = ones(1, size(points_roi, 1));    % Is it covered by water? 1 = no. 0 = yes.

% Loop through images
for k = 1 : lastFrame - skipFrames - 1
%for k = 1 : 6  % to check looping for a few frames 

    if ((mod(k, every_frame_v) >= every_frame_v - averaging) || (mod(k, every_frame_v) <= averaging) ||(mod(k, every_frame_v) == 0) ||  k <= averaging + 1 || k >= lastFrame - skipFrames - 1 - averaging)

        res_name = '';

        % If we only want first and last frames
        if (onlyLast == 1)

            % Last frame for a fill event
            if (fill_emp == 0)
                index = size(imgFiles, 1);

            % Last frame for an empty event
            else
                index = 1;
            end

            res_name = strcat(imgFolder_name, '_First_Last');

        % If we want the in-between frames also
        else
            % Frame for a fill event
            if (fill_emp == 0)
                index = start + k;
                res_name = strcat(imgFolder_name, '_', string(firstFrameNum + k));
                disp(' ')
                disp(strcat('frame number: ', string(firstFrameNum + k)));

            % Frame for an empty event
            else
                index = start - k;
                res_name = strcat(imgFolder_name, '_', string(firstFrameNum - k));
                disp(' ')
                disp(strcat('frame number: ', string(firstFrameNum - k)));
            end
        end

        % File Paths of two images to compare between - change to your folder
        img1 = fullfile(filepath, imgFiles(start).name);
        img2 = fullfile(filepath, imgFiles(index).name);

        % Load images
        A = imread(char(strcat(img1)));
        B = imread(char(strcat(img2)));
        
        % Undistort image
        if (undistort == 1)
            [A, newOrigin_a] = undistortImage(A, params, 'OutputView', 'full');
            [B, newOrigin_b] = undistortImage(B, params, 'OutputView', 'full');
        end

        % Rotate (if necessary)
        A = imrotate(A,rotation);
        B = imrotate(B,rotation);

        % HIST Eqivalization to normalize global lighting / sunlight / shadow
        % effect
        
        if hist_eq == 1   %ROIs change color
            % A = histeq(A);
            % B = histeq(B);

            % Apply light normalization to Image A
            resImg_A = lightNormal(A, roi_store, BlockSize_nonstatic, BlockSize_static);
            A = resImg_A;
            
            % Apply light normalization to Image B
            resImg_B = lightNormal(B, roi_store, BlockSize_nonstatic, BlockSize_static);
            B = resImg_B;

            % figure
            % imshow(A);
            % figure
            % imshow(B);
        end

        A_gray = rgb2gray(A); % "fixed" image
        B_gray = rgb2gray(B); % image to register

        if (mod(k, every_frame_v) == 0 || k == 1 || k == lastFrame - skipFrames - 1)
            store_res_names = [store_res_names, res_name];
            store_img_index = [store_img_index, index];
        end


%         % Light normalization for only the ROI's instead of the whole image
%         if hist_eq == 1
% %             A = histeq(A);
% %             B = histeq(B);
%             for i=1:(size(roi_store,2))
%                 ROI_A = A_gray(roi_store{i}(2):roi_store{i}(2)+roi_store{i}(4), roi_store{i}(1):roi_store{i}(1)+roi_store{i}(3));
%                 histEqROI_A = histeq(ROI_A);
%                 resImg_A = A_gray;
%                 resImg_A(roi_store{i}(2):roi_store{i}(2)+roi_store{i}(4), roi_store{i}(1):roi_store{i}(1)+roi_store{i}(3)) = histEqROI_A;
%                 A_gray = resImg_A;
%             end
%             for i=1:(size(roi_store,2))
%                 ROI_B = B_gray(roi_store{i}(2):roi_store{i}(2)+roi_store{i}(4), roi_store{i}(1):roi_store{i}(1)+roi_store{i}(3));
%                 histEqROI_B = histeq(ROI_B);
%                 resImg_B = B_gray;
%                 resImg_B(roi_store{i}(2):roi_store{i}(2)+roi_store{i}(4), roi_store{i}(1):roi_store{i}(1)+roi_store{i}(3)) = histEqROI_B;
%                 B_gray = resImg_B;
%             end
%             figure
%             imshow(A_gray);
%             figure
%             imshow(B_gray);
%         end

        % Homography BEFORE KLT Tracking
        if compensation_method == 2
            for  i=1:num_static_roi

                StaticObjectRegion = static_roi_store{i};

                static_points = Feature_detector(A_ref,StaticObjectRegion,detector_option);
                [static_features, static_points] = extractFeatures(rgb2gray(A_ref),static_points);

                %motion_points = detectSURFFeatures(rgb2gray(B),'ROI',StaticObjectRegion);
                motion_points = Feature_detector(B,StaticObjectRegion,detector_option);
                [motion_features, motion_points] = extractFeatures(rgb2gray(B),motion_points);

                indexPairs = matchFeatures(motion_features, static_features, 'Unique', true);     

                if i == 1
                        movingPoints = motion_points(indexPairs(:,1), :);
                        fixedPoints = static_points(indexPairs(:,2), :);
                    else
                        movingPoints = vertcat(movingPoints , motion_points(indexPairs(:,1), :));
                        fixedPoints = vertcat(fixedPoints , static_points(indexPairs(:,2), :));
                end
            end

            % Matlab estimate H using 'rigid' transform assumption
            %%%% Change 'rigid' to 'projective' if warp Looks Distorted
            [tforms,inlierIndex,status] = estimateGeometricTransform2D(movingPoints, fixedPoints,...
            'rigid', 'Confidence', 99.9, 'MaxNumTrials', 2000,'MaxDistance',threshold);
            disp(['status: ',num2str(status)])
            disp(['number of inliers: ',num2str(sum(inlierIndex))])

            % If there's not enough matched points, then not much can
            % be done about it. Continue to next loop
            if (status == 1)
                disp('matchedPoints1 and matchedPoints2 inputs do not contain enough points. Moving on to next frame...')
                continue;
            end

            % If there's not enough inliers, then try to fix it
            if (sum(inlierIndex) < homographyinliers)
                disp('Not enough inliers found. Attempting to fix...')
            end
            turn = 1;
            while (sum(inlierIndex) < homographyinliers)
                if (turn == 10)
                    break
                else
                    error_message = strcat('Not enough inliers found. Attempt #', string(turn));
                    disp(error_message)
                    Threshold = threshold*(2^turn);
                    [tforms,inlierIndex,status] = estimateGeometricTransform2D(movingPoints, fixedPoints,...
                'rigid', 'Confidence', 99.9, 'MaxNumTrials', 2000,'MaxDistance',Threshold);
                    disp(['status: ',num2str(status)])
                    disp(['number of inliers: ',num2str(sum(inlierIndex))])
                    turn = turn+1;
                end
                if (onlyLast == 1)
                    break
                end
            end

            if (sum(inlierIndex) < homographyinliers || status ~= 0)
                disp('Not enough inliers could be found even after increasing threshold. Moving to next frame...')
                continue;
            end
            
            % Warp operation : same as Outputniew
            B_wrap = imwarp(B,tforms,'OutputView',imref2d(size(A_gray)));

            if compensation_method ~= 3 % METHOD 3
                % H transform after KLT tracking on tracked Points  Method
                B = B_wrap;
            end


        % Imregister Method
        elseif compensation_method == 1
            [optimizer, metric] = imregconfig('monomodal');
            optimizer.GradientMagnitudeTolerance = 1E-4;
            optimizer.MinimumStepLength = 1E-5;
            optimizer.MaximumStepLength = 6.25E-2;
            optimizer.MaximumIterations = 100;
            optimizer.RelaxationFactor = 5E-1;
            BregionImg = B_gray(region(2):(region(2)+region(4)),region(1):(region(1)+region(3)));
            tforms = imregtform(BregionImg, AregionImg, 'rigid', optimizer, metric);
            
            % Warp operation : same as Outputniew
            B_wrap = imwarp(B,tforms,'OutputView',imref2d(size(A_gray)));
        end
        

        % KLT Tracker
        % Create a tracker object.
        tracker_static = vision.PointTracker('MaxBidirectionalError',KLTMaxBidirectionalError,...
                                       'BlockSize',[BlockSize_static,BlockSize_static],'NumPyramidLevels',3,...
                                       'MaxIterations',30);
        tracker_nonstatic = vision.PointTracker('MaxBidirectionalError',KLTMaxBidirectionalError,...
                                       'BlockSize',[BlockSize_nonstatic,BlockSize_nonstatic],'NumPyramidLevels',3,...
                                       'MaxIterations',30);

        % First / refernce Frame Feature Points
        points_in_static=double((points_roi(1:total_static_roi,:)));        
        objectFrame = A;

        % Initialize the tracker.
        initialize(tracker_static,points_in_static,objectFrame);
        
        % Track with KLT on Frame B
        [points2_static,validity_static,scores_static] = tracker_static(B);
        points_in_nonstatic=double((points_roi(total_static_roi+1:end,:)));         
        objectFrame = A;

        % Initialize the tracker.
        initialize(tracker_nonstatic,points_in_nonstatic,objectFrame);

        % Track with KLT on Frame B
        [points2_nonstatic,validity_nonstatic,scores_nonstatic] = tracker_nonstatic(B);
        points_in=double((points_roi));
        points2=[points2_static;points2_nonstatic];
        validity=[validity_static;validity_nonstatic];
        scores=[scores_static;scores_nonstatic];
          
        %     % Display with Imfuse overlay and distance
        %     % Quiver plot between 2 frames   
        %     pts1 = double(points_in(validity,:));
        %     pts2 = double(points2(validity,:));
        %     dp  = pts2 - pts1;
        % 
        %     out = insertMarker(B,pts2,'+','Color','red','size',20);
        %     out = insertMarker(out,pts1,'+','Color','green','size',20);%% 1st frame
        %     figure;
        %         imshow(out);
        %         hold on 
        %         quiver(pts1(:,1),pts1(:,2),dp(:,1),dp(:,2),2,'LineWidth',1.5 )
        %         axis equal
        %         impixelinfo;

        % Keep only valid points
        pts1 = double(points_in) .* validity;
        pts2 = double(points2) .* validity;

        % Homography AFTER KLT Tracking
        if compensation_method == 3
            fixedPoints = pts1(1:total_static_roi,:);
            fixedPoints (all(fixedPoints == 0, 2),:) = [];
            movingPoints = pts2(1:total_static_roi,:);
            movingPoints (all(movingPoints == 0, 2),:) = [];
            [tforms,inlierIdx,status] = estimateGeometricTransform2D(movingPoints,fixedPoints,'rigid','Confidence', 99.9, 'MaxNumTrials', 2000,'MaxDistance',threshold);
            H_mat = (tforms.T)';
            pts2_org = vertcat(double(points2)',ones(1,size(points2,1)));
            pts2_H  = H_mat * pts2_org;
            pts2_H_2d = pts2_H(1:2,:)';
            pts2 = pts2_H_2d;
        end    

        store_pts2(loop_num, :, :) = pts2;    % Store pts2
        store_pts1 = pts1;                    % Store pts1

        dp  = (pts2 - pts1) .* horzcat(validity, validity);
        
        % calculate abs displacemnt
        disp_pix = abs(sqrt(dp(:,1).^2 + dp(:,2).^2));

        % Check if points are still valid (aka are they covered by water?)
        if (loop_num == 1)
            temp_water_cover = water_cover(1);
        else
            temp_water_cover = water_cover(loop_num-1, :);
        end
        temp_water_cover(disp_pix .* point_scales' > vanish_thres) = 0;
        water_cover(loop_num, :) = temp_water_cover;
        point_scales(temp_water_cover == 0) = 0;

        disp_x_pix = dp(:,1);
        disp_y_pix = dp(:,2);

        % stack disp data to be saved
        try
            store_disp_x_pix(loop_num,:) = disp_x_pix';
            store_disp_y_pix(loop_num,:) = disp_y_pix';
            store_scores(loop_num,:) = scores';
            store_frameNum(loop_num,:) = firstFrameNum + k;
        catch
            disp('Displacement values were not able to be added to matrix with all values. Moving to next frame...')
            continue
        end

        loop_num = loop_num + 1;

    else
        continue;
    end
    
end


% Averaging across time and for ROI
if (onlyLast == 1)
    averaging = 0;
end
[store_avgX_pix, store_avgY_pix, store_stdX_pix, store_stdY_pix, store_frameNum_roi] = averageROIDispl(store_disp_x_pix, store_disp_y_pix, store_scores, averaging, water_cover, num_nonstatic_roi, num_static_roi, point_num, numOutlierRounds, numStd, stdThr, store_frameNum);


% Plot result video
createResultVideo(result_filePath, imgFolder_name, whatFrameR, store_avgX_pix, store_avgY_pix, roi_centroid, q_scale, rgb2gray(A_ref), store_img_index, store_pts2, store_pts1, every_fig, store_res_names, filepath, imgFiles, onlyLast)


% keep saving in every loop, if interruptions in between some data gets
% saved
% save disp___filename 
try
    disp_data_name = fullfile(result_filePath, strcat('disp_data_', imgFolder_name));
    save(disp_data_name, 'store_disp_x_pix', 'store_disp_y_pix', 'store_scores', 'store_frameNum', 'store_img_index', 'store_res_names', 'store_pts2', 'store_pts1', 'total_static_roi', 'score_min', 'result_filePath', 'imgFolder_name', 'points_roi', 'point_scales', 'roi_centroid', 'store_avgX_pix', 'store_avgY_pix', 'store_stdX_pix', 'store_stdY_pix', 'num_static_roi', 'num_nonstatic_roi', 'roi_scales', 'good_rois', 'store_frameNum_roi')
catch
    disp('End of program has been reached and no displacements were successfully found. Stopping program...')
    return
end



%% Plot Displacement Vectors for each of the tracked Points wrt frame number

% If we're not just doing the first and last frame, then make plots
if (onlyLast ~= 1)

    %%% Per ROI plots
    
    % Static 
    % Static ROI Mean Displ for X and Y (pix)
    PlottingROIMethod(store_frameNum_roi, store_avgX_pix, 1, num_static_roi, result_filePath, imgFolder_name)
    PlottingROIMethod(store_frameNum_roi, store_avgY_pix, 1, num_static_roi, result_filePath, imgFolder_name)
 
    % Static ROI Std Displ for X and Y (pix)
    PlottingROIMethod(store_frameNum_roi, store_stdX_pix, 1, num_static_roi, result_filePath, imgFolder_name)
    PlottingROIMethod(store_frameNum_roi, store_stdY_pix, 1, num_static_roi, result_filePath, imgFolder_name)

    % Non Static
    % NonStatic ROI Mean Displ for X and Y (pix)
    PlottingROIMethod(store_frameNum_roi, store_avgX_pix, 0, num_static_roi, result_filePath, imgFolder_name)
    PlottingROIMethod(store_frameNum_roi, store_avgY_pix, 0, num_static_roi, result_filePath, imgFolder_name)

    % NonStatic ROI Std Displ for X and Y (pix)
    PlottingROIMethod(store_frameNum_roi, store_stdX_pix, 0, num_static_roi, result_filePath, imgFolder_name)
    PlottingROIMethod(store_frameNum_roi, store_stdY_pix, 0, num_static_roi, result_filePath, imgFolder_name)

    % NonStatic ROI Mean Displ for X and Y (mm)
    store_avgX_mm = store_avgX_pix .* roi_scales;   % Displ converted from pix to mm
    store_avgY_mm = store_avgY_pix .* roi_scales;   % Displ converted from pix to mm
    PlottingROIMethod(store_frameNum_roi, store_avgX_mm, 0, num_static_roi, result_filePath, imgFolder_name)
    PlottingROIMethod(store_frameNum_roi, store_avgY_mm, 0, num_static_roi, result_filePath, imgFolder_name)

    % NonStatic ROI Std Displ for X and Y (mm)
    store_stdX_mm = store_stdX_pix .* roi_scales;   % Displ converted from pix to mm
    store_stdY_mm = store_stdY_pix .* roi_scales;   % Displ converted from pix to mm
    PlottingROIMethod(store_frameNum_roi, store_stdX_mm, 0, num_static_roi, result_filePath, imgFolder_name)
    PlottingROIMethod(store_frameNum_roi, store_stdY_mm, 0, num_static_roi, result_filePath, imgFolder_name)

    
    %%% Per point plots
    
    % Displacement X and Y plots
    %PlottingPointMethod(store_frameNum, store_dispx2_mm, 0, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    %PlottingPointMethod(store_frameNum, store_dispy2_mm, 0, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    
    % Displacement Plots (Static)
    %PlottingPointMethod(store_frameNum, store_disp_mm, 1, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    %PlottingPointMethod(store_frameNum, store_disp_pix, 1, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    
    % Displacement Plots (Non Static)
    %PlottingPointMethod(store_frameNum, store_disp_mm, 0, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    %PlottingPointMethod(store_frameNum, store_disp_pix, 0, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)

    % Angle Plots
    %PlottingPointMethod(store_frameNum, store_ang_rad, 1, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    %PlottingPointMethod(store_frameNum, store_ang_rad, 0, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)

    % Score Plots
    %PlottingPointMethod(store_frameNum, store_scores2, 1, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)
    %PlottingPointMethod(store_frameNum, store_scores2, 0, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)


end



%% FUNCTIONS

% Feature selection method Function
function detected_points = Feature_detector(img,input_roi,detector_option)

    if detector_option == 1 
        % Good features to track
        detected_points = detectMinEigenFeatures(rgb2gray(img),'ROI',input_roi);

    elseif detector_option == 2
        % detectBRISKFeatures
        detected_points = detectBRISKFeatures(rgb2gray(img),'ROI',input_roi);

    elseif detector_option == 3
        % SURF
        detected_points = detectSURFFeatures(rgb2gray(img),'ROI',input_roi);

    elseif detector_option == 4
        %Hrrsis corner
        detected_points = detectHarrisFeatures(rgb2gray(img),'ROI',input_roi);

    else
        % DEFAULT -  %Hrrsis corner
        %detected_points = detectHarrisFeatures(rgb2gray(img),'ROI',input_roi);
        detected_points = detectHarrisFeatures(rgb2gray(img),'ROI',input_roi,'MinQuality',0.001);
    end

end



% Get ROI Function
function [roi_store2, roi_scales2, num_roi2, good_rois] = GetROI(node_file, node_file_roi_halfSize, A_ref, scale_method, f_pix, nodePixel, nodeCoords, camCoords, blenderAdj)

    roi_store2 ={};   % Nonstatic roi
    roi_scales2 = [];    % pix2mm for each roi
    node_coords2 = [];   % usable node coords for scales
    
    % Use file nodes for ROI's
    if (node_file == 1)

        tot_roi = size(nodePixel, 1);   % Total number of nodes
        num_roi2 = 0;                   % Number of usable nodes
        good_rois = zeros(1, size(nodePixel, 1));  % Array with good rois. 0 = bad. 1 = good.

        for i= 1 : tot_roi

            min_x = nodePixel(i, 1) - node_file_roi_halfSize;   % Min x for roi
            min_y = nodePixel(i, 2) - node_file_roi_halfSize;   % Min y for roi

            width = 2 * node_file_roi_halfSize;     % Width of roi
            height = 2 * node_file_roi_halfSize;     % Height of roi

            % If roi goes off the image, then skip to the next one
            if (min_x < 1 || min_y < 1 || min_x + width > size(A_ref, 2) || min_y + height > size(A_ref, 1))
                continue
            end

            % Make sure the roi won't go off the image
%             % Left cut off
%             if (min_x < 1)
%                 width = width + (1 - min_x);
%             end
% 
%             % Top cut off
%             if (min_y < 1)
%                 height = height + (1 - min_y);
%             end
% 
%             % Right cut off
%             if (min_x + width > size(A_ref, 2))
%                 min_x = min_x - (min_x + width - size(A_ref, 2));
%             end
% 
%             % Bottom cut off
%             if (min_y + height > size(A_ref, 1))
%                 min_y = min_y - (min_y + height - size(A_ref, 1));
%             end

            num_roi2 = num_roi2 + 1;    % Increase number of roi
            good_rois(i) = 1;

            objectRegion = [min_x, min_y, width, height];
            roi_store2{num_roi2} = objectRegion;
            node_coords2(num_roi2, :) = nodeCoords(i, :);

        end

        % Scales
        prev = 0;
        roi_scales2 = ScalingMethod(3, f_pix, prev, node_coords2, camCoords, blenderAdj);
        roi_scales2 = roi_scales2';         % Transpose to be consisten with static ones


    % Choose ROI's by hand
    elseif (node_file == 0)

        % ROI select for KLT Tracking Region selection
        % use mouse
        
        figure; 
        imshow(A_ref);
        title('Select ROI for KLT Tracking')
        disp('Zoom in as required and Press ENTER')
        pause;%press enter to continue

        for i= 1 : num_roi2
            % Get region of interest
            polygon_roi = drawrectangle;
            objectRegion = reshape(polygon_roi.Position', 1, []);
            objectImage = insertShape(objectFrame,'Rectangle',objectRegion,'Color','red');
            objectFrame  = objectImage;
            imshow(objectFrame);  
            roi_store2{i} = objectRegion;
            
            prev = 0;
            if (i ~= 1)
                prev = roi_scales2(i-1);
            end
            roi_scales2(i) = ScalingMethod(scale_method, f_pix, prev, nodeCoords, camCoords, blenderAdj);
        end

    end
    
end



% Scaling Method Function
% Scale method 1 = Enter the known distance between 2 points in the image
% Scale method 2 = Enter the known distance from the camera to the roi
% Scale method 3 = Automatically find scale based on 3D global node
%   coordinates
function scale = ScalingMethod(scale_method, f_pix, prev, nodeCoords, camCoords, blenderAdj)

    % Scale method = 3
    if scale_method == 3
        x_dist = nodeCoords(:, 1) + blenderAdj(1) - camCoords(1);   % Global x distance from camera to point
        y_dist = nodeCoords(:, 2) + blenderAdj(2) - camCoords(2);   % Global y distance from camera to point
        z_dist = nodeCoords(:, 3) + blenderAdj(3) - camCoords(3);   % Global z distance from camera to point

        scale = sqrt(x_dist.^2 + y_dist.^2 + z_dist.^2 ) / f_pix;   % Calculate scale (m/pix)
        scale = scale * 1000;                                       % Convert scale to mm/pix

    else
        % Enter Known distance in (m)
        % Figure out distance prompt based on scale method
        dist_prompt = '';
        if scale_method == 1
            dist_prompt = 'Enter the known distance in meters (or type -1 if it is the same as previous): ';
        elseif scale_method == 2
            dist_prompt = 'Enter the known distance from the camera to the roi in meters: ';
        end
    
        dist_m = input(dist_prompt);
        dist_mm  =  dist_m * 1000;
        % If it's a number, then we are inputing a new scale factor
        % If it's "-1", then we are just using the previous pix2mm number
        if (dist_m ~= -1)  
            if scale_method == 1
                % Select 2 points of known distance ( in the same plane of motion if possible) for ROI
                l = drawline('StripeColor','red');
                % Estimate scaling factor
                pos_diff = (l.Position(2,:) - l.Position(1,:)); 
                dist_pix = abs(sqrt(pos_diff(1).^2 + pos_diff(2).^2));
                pix2mm = dist_mm / dist_pix;
            end        
        else
            pix2mm = prev;
        end
    
        if scale_method == 1
            scale = pix2mm;
        else
            % metadata - focal length method
            scale = dist_m*1000/f_pix; %  Z_roi_mm 
        end     

    end

end



% Light Normalization Function
function resImg_A = lightNormal(A_ref, roi_store, BlockSize_nonstatic, BlockSize_static)

    resImg_A = A_ref;

    % Find blocksize (aka amount we need to let the light normalization
    % extend
    if (BlockSize_nonstatic > BlockSize_static)
        blockSize = BlockSize_nonstatic + 2;
    else
        blockSize = BlockSize_static + 2;
    end

    % Loop through each ROI and apply light normalization
    for i=1:(size(roi_store,2))

        % Make light normalization area slightly larger than roi area
        start1 = round(roi_store{i}(2) - blockSize);
        if (start1 < 1)
            start1 = 1;
        end

        end1 = round(roi_store{i}(2) + roi_store{i}(4) + blockSize);
        if (end1 > size(A_ref, 1))
            end1 = size(A_ref, 1);
        end

        start2 = round(roi_store{i}(1) - blockSize);
        if (start2 < 1)
            start2 = 1;
        end

        end2 = round(roi_store{i}(1) + roi_store{i}(3) + blockSize);
        if (end2 > size(A_ref, 2))
            end2 = size(A_ref, 2);
        end

%         disp(start1)
%         disp(end1)
%         disp(start2)
%         disp(end2)
        
        % Apply Light Normalization
        ROI_A = A_ref(start1:end1, start2:end2,:);
        histEqROI_A = histeq(ROI_A);
        resImg_A(start1:end1, start2:end2,:) = histEqROI_A;
        
%         figure
%         imshow(A_ref);

    end
end



% Plotting Function based on Points
% x = x axis values (usually = frame_vals)
% y = y axis values (usually = displacement, angle, or scores)
% static = Plot only static points? 0 = no. 1 = yes. Else = Plot all points
% total_static_roi = Total number of static points.
% store_scores2 = Scores for all points
% score_min = Minimum score for a point to plot
% result_filePath = Where to save files.
function PlottingPointMethod(x, y, static, total_static_roi, store_scores2, score_min, result_filePath, imgFolder_name)

    % Figure out label names
    x_name = inputname(1);
    y_name = inputname(2);
    [x_axis_name, x_orig_name] = AxisName(x_name);
    [y_axis_name, y_orig_name] = AxisName(y_name);
    title_name = y_axis_name;
    plotTitle = y_orig_name;
    
    % Figure out what points to plot
    beginPoints = -1;
    lastPoints = -1;
    if static == 0
        beginPoints = total_static_roi + 1;
        lastPoints = size(y,2);
        title_name = strcat(title_name, ' for Non-Static Points');
        plotTitle = strcat(plotTitle, '_nonStaticPts');
    elseif static == 1
        beginPoints = 1;
        lastPoints = total_static_roi;
        title_name = strcat(title_name, ' for Static Points');
        plotTitle = strcat(plotTitle, '_staticPts');
    else
        beginPoints = 1;
        lastPoints = size(y,2);
        title_name = strcat(title_name, ' for All Points');
        plotTitle = strcat(plotTitle, '_allPts');
    end
    
    f_plot = figure;
        for i = beginPoints:lastPoints
            % Extract valid displacements and frame numbers
            y_row = y(:,i);
            y_vals = y_row(store_scores2(:,i) > score_min);
            x_vals = x(store_scores2(:,i) > score_min);

            txt = ['P',num2str(i)];
            plot(x_vals, y_vals,'DisplayName',txt,'LineWidth',2)
            hold on
        end
        hold off
        legend show
        legend('Location','eastoutside')
        title(title_name)
        grid on
        xlabel(x_axis_name)
        ylabel(y_axis_name)
            set(gcf,'color','w');

    plot_name = fullfile(result_filePath, strcat(imgFolder_name, '_', plotTitle, '.png'));
    saveas(f_plot, plot_name);    % Save plot as .png
    
    fig_name = fullfile(result_filePath, strcat(imgFolder_name, '_', plotTitle, '.fig'));
    savefig(f_plot, fig_name)     % Save plot as .fig
    
end



% Average ROI displacments (in Pixels)
% store_disp_x_pix = pixel displacements in x direction.
% store_disp_y_pix = pixel displacements in y direction.
% store_scores = scores for tracking
% averaging = how many pixels on each side you want to average
% water_cover = is the displacement covered by water? 
% num_nonstatic_roi = number of nonstatic roi
% num_static_roi = number of static roi
% point_num = number of points in each roi
% numOutlierRounds = number of rounds to get of outliers
% numStd = number of std away from mean that points need to be to not get
%       thrown out
% stdThr = standard deviation threshold
% store_frameNum = frame numbers (there should be one per time frame)
function [store_avgX_pix, store_avgY_pix, store_stdX_pix, store_stdY_pix, store_frameNum_roi] = averageROIDispl(store_disp_x_pix, store_disp_y_pix, store_scores, averaging, water_cover, num_nonstatic_roi, num_static_roi, point_num, numOutlierRounds, numStd, stdThr, store_frameNum)
    
    store_avgX_pix = [];    % Average x pixel displacement
    store_avgY_pix = [];    % Average y pixel displacement
    store_stdX_pix = [];    % Average x pixel std
    store_stdY_pix = [];    % Average y pixel std
    store_frameNum_roi = []; % Frame Numbers after averaging over time

    count = 1;  % Counter

    while count <= size(store_disp_x_pix, 1)

        % Initialize arrays
        row_avgX_pix = [];    % Average x pixel displacement
        row_avgY_pix = [];    % Average y pixel displacement
        row_stdX_pix = [];    % Average x pixel std
        row_stdY_pix = [];    % Average y pixel std

        start_ind = count - averaging;  % Start index
        end_ind = count + averaging;    % End index

        store_frameNum_roi = [store_frameNum_roi, store_frameNum(count)];

        % If it's the beginning or end, then our start and stop points will
        % need to be tweaked
        if (start_ind <= 0)
            start_ind = 1;
        elseif (end_ind > size(store_disp_x_pix, 1))
            end_ind = size(store_disp_x_pix, 1);
        end

        % Get matrix of displacements we want average of
        temp_dispx = store_disp_x_pix(start_ind:end_ind, :);
        temp_dispy = store_disp_y_pix(start_ind:end_ind, :);
        temp_score = store_scores(start_ind:end_ind, :);
        
        % Find average by dividing by number of nonzero elements
        % Sum matrices to get a column vector (each row is one point)
        disp_x_pix2 = sum(temp_dispx, 1)';
        disp_y_pix2 = sum(temp_dispy, 1)';
        score_2 = sum(temp_score, 1)';
        
        % Count number of zeros in displ array
        numZeros_x = sum(temp_dispx == 0)';
        numZeros_y = sum(temp_dispy == 0)';
        numZeros_score = sum(temp_score == 0)';
        
        % Divide each element by number of non-zero elements in array
        disp_x_pix2 = disp_x_pix2 ./ (size(temp_dispx, 1) - numZeros_x);
        disp_y_pix2 = disp_y_pix2 ./ (size(temp_dispy, 1) - numZeros_y);
        score_2 = score_2 ./ (size(temp_score, 1) - numZeros_score);
        
        % Account for water cover
        disp_x_pix2(water_cover(count, :) == 0) = 0;
        disp_y_pix2(water_cover(count, :) == 0) = 0;

        %size(disp_x_pix2)
        %size(disp_y_pix2)
        
        % Find average and std for each roi
        for j = 1 : num_nonstatic_roi+num_static_roi

            if j == 1
                startPoint = 1;                         % Start index of points of current roi
            else
                startPoint = sum(point_num(1:j-1));     % Start index of points of current roi
            end
            endPoint = sum(point_num(1:j));             % End index of points of current roi
        
            % Extract array 
            temp_roiDisplX = disp_x_pix2(startPoint : endPoint); % X displacements array for roi
            temp_roiDisplY = disp_y_pix2(startPoint : endPoint); % Y displacements array for roi
            
            % Get rid of non valid points
            temp_roiDisplX(temp_roiDisplX == 0) = [];
            temp_roiDisplY(temp_roiDisplY == 0) = [];
            
            % Find Initial Averages
            avgX_pix = mean(temp_roiDisplX);
            avgY_pix = mean(temp_roiDisplY);

            %size(avgX_pix)
            %size(avgY_pix)

            % Find Initial Standard Deviations
            stdX_pix = std(temp_roiDisplX);
            stdY_pix = std(temp_roiDisplY);
            
            % Go through and get rid of outliers
            for i = 1:numOutlierRounds

                % Get rid of displacements 2 std's above or below mean
                temp_roiDisplX(temp_roiDisplX > avgX_pix + numStd * stdX_pix) = [];
                temp_roiDisplX(temp_roiDisplX < avgX_pix - numStd * stdX_pix) = [];
                
                temp_roiDisplY(temp_roiDisplY > avgY_pix + numStd * stdY_pix) = [];
                temp_roiDisplY(temp_roiDisplY < avgY_pix - numStd * stdY_pix) = [];
                
                % Find Updated Averages
                avgX_pix = mean(temp_roiDisplX);
                avgY_pix = mean(temp_roiDisplY);

                % Find Updated Standard Deviations
                stdX_pix = std(temp_roiDisplX);
                stdY_pix = std(temp_roiDisplY);
                
            end

            % Assign NaN to roi's with std's above the threshold
            avgX_pix(stdX_pix > stdThr) = NaN;
            avgY_pix(stdY_pix > stdThr) = NaN;
            stdX_pix(stdX_pix > stdThr) = NaN;
            stdY_pix(stdY_pix > stdThr) = NaN;
            
            % Replace NaN avg and std with NaN
            avgX_pix(isnan(avgX_pix)) = NaN;
            avgY_pix(isnan(avgY_pix)) = NaN;
            stdX_pix(isnan(stdX_pix)) = NaN;
            stdY_pix(isnan(stdY_pix)) = NaN;
            
            % Add averages to avg array vector for one time frame
            row_avgX_pix = [row_avgX_pix, avgX_pix];
            row_avgY_pix = [row_avgY_pix, avgY_pix];

            % Add standard deviations to std array
            row_stdX_pix = [row_stdX_pix, stdX_pix];
            row_stdY_pix = [row_stdY_pix, stdY_pix];
            
        end

        % Add averages to avg array
        store_avgX_pix = [store_avgX_pix; row_avgX_pix];
        store_avgY_pix = [store_avgY_pix; row_avgY_pix];

        % Add standard deviations to std array
        store_stdX_pix = [store_stdX_pix; row_stdX_pix];
        store_stdY_pix = [store_stdY_pix; row_stdY_pix];

        % Update count
        count = count + averaging * 2;

        if (averaging == 0)
            break;
        end

    end
end



% Create Result Video 
% result_filePath = file path where results will go
% imgFolder_name = folder name where all the input images are
% whatFrameR = desired video framerate
% store_avgX_pix = average pixel displacement in x direction
% store_avgY_pix = average pixel displacement in y direction
% roi_centroid = roi centroids coordinates (# roi's x 2)
% q_scale = how much to scale the quiver plot arrows
% A_gray = reference image in grayscale
% img_index = all image indices for imgFile for video
% pts2 = points locations for images (# images x # points x 2)
% pts1 = points locations for reference image (# points x 2)
% every_fig = how often a figure is wanting to be saved
% res_name = string frame number array (# images x 1)
% filepath = filepath to images
% imgFiles = imgFile array with all images in folder
% onlyLast = Calculate only last frame? Or all?
function createResultVideo(result_filePath, imgFolder_name, whatFrameR, store_avgX_pix, store_avgY_pix, roi_centroid, q_scale, A_gray, img_index, pts2, pts1, every_fig, res_name, filepath, imgFiles, onlyLast)

    % Create VideoWriter object to make video
    if (onlyLast ~= 1)
        video_name = fullfile(result_filePath, strcat(imgFolder_name, '_comp_v2.mp4'));
        v = VideoWriter(video_name, 'MPEG-4');
        v.Quality = 99;
        v.FrameRate = whatFrameR;
        open(v)
    end

    % Go through all images and plot result
    for i = 1:length(img_index)

        % Quiver Plot text 
        disp_pix2 = sqrt(store_avgX_pix(i, :).^2 + store_avgY_pix(i, :).^2);
        num_pts = length(disp_pix2);
        text_str = cell(num_pts,1);
        arrow_num = cell(num_pts,1);
        for j = 1:num_pts
            text_str{j}= [num2str(disp_pix2(j),'%0.0f'),'mm'];
            arrow_num{j}= [num2str((j),'%0d')];
        end
        
    %     % ROI statistics text boxes
    %     roi_stats_str = cell(num_static_roi+num_nonstatic_roi, 1);
    %     for i=1:num_static_roi+num_nonstatic_roi
    %         roi_stats_str{i} = ['Average X: ', num2str(store_avgX_pix(loop_num,i), '%0.2f'), ...
    %             '. STD X: ', num2str(store_stdX_pix(loop_num,i), '%0.2f'), ...
    %             '. Average Y: ', num2str(store_avgY_pix(loop_num,i), '%0.2f'), ...
    %             '. STD Y: ', num2str(store_stdY_pix(loop_num,i), '%0.2f')];
    %     end

        % Load image 2
        img2 = fullfile(filepath, imgFiles(img_index(i)).name);
        B = imread(char(strcat(img2)));
        B_gray = rgb2gray(B);
        
        C = imfuse(A_gray, B_gray, 'falsecolor', 'Scaling', 'joint', 'ColorChannels', [1 2 0]);
    
        temp_pts2(:, 1) = pts2(i, :, 1);
        temp_pts2(:, 2) = pts2(i, :, 2);
        out = insertMarker(C,temp_pts2,'+','Color','red','size',20);
        out = insertMarker(out,pts1,'+','Color','green','size',20);%% 1st frame
        %out = insertText(out, pts1 + 25 , text_str,'FontSize',50,'BoxColor','red'); % Value
        %out = insertText(out, pts1 - 200 , arrow_num,'FontSize',50,'BoxColor','w'); % Point number
    
        % Plot scales for each point
    %         if loop_num == 1
    %             pix_res_str = cell(num_pts,1);
    %             for i=1:num_pts
    %                 pix_res_str{i}= [num2str(point_scales(i),'%0.0f'),'mm-1pix'];
    %             end
    %             out = insertText(out, pts1 - 500 , pix_res_str,'FontSize',50,'BoxColor','blue');
    %         end
    
        % ROI Statistic labels
        %out = insertText(out, roi_centroid + 25, roi_stats_str,'FontSize',50,'BoxColor','w');
        
        f1= figure;
        imshow(out, 'Border', 'tight');
        hold on 
        
        % For each point
        %q=quiver(pts1(:,1),pts1(:,2),q_scale*disp_x_mm2,q_scale*disp_y_mm2,0,'LineWidth',1.5 );
        
        % For each roi, Plot
        q=quiver(roi_centroid(:,1),roi_centroid(:,2),q_scale.*store_avgX_pix(i,:)',q_scale.*store_avgY_pix(i,:)',0,'LineWidth',1.5 );
        
        q.AutoScale = 'off';
        %axis equal
        %impixelinfo;
        pause(0.01);
    
        % Save current fig as .fig file
        if ((mod(i, every_fig) == 0) ||  i == 1 || i == size(img_index, 1))
            fig_outName = fullfile(result_filePath,strcat(res_name(i), '.fig'));
            savefig(f1, fig_outName)
        end

        % Save result photo and add to video 
        % Save to result folder and add to video
        % Save current fig as .png file
        outName = fullfile(result_filePath,strcat(res_name(i), '.png'));
        saveas(f1, outName);    % f1 = w/ quiver.  f2 = w/o quiver

        if (onlyLast ~= 1)
            output_img = imread(outName);
            writeVideo(v, output_img)
        else
            break
        end

        close all
    
    end

    % Close video
    if (onlyLast ~= 1)
        close(v)
    end

end



% Plotting Function based on ROI
% x = x axis values (usually = frame_vals)
% y = y axis values (usually = displacement, angle, or scores)
% static = Plot only static points? 0 = no. 1 = yes. Else = Plot all points
% num_static_roi = Total number of static ROI's.
% result_filePath = Where to save files.
function PlottingROIMethod(x, y, static, num_static_roi, result_filePath, imgFolder_name)

    % Figure out label names
    x_name = inputname(1);
    y_name = inputname(2);
    [x_axis_name, x_orig_name] = AxisName(x_name);
    [y_axis_name, y_orig_name] = AxisName(y_name);
    title_name = y_axis_name;
    plotTitle = y_orig_name;
    
    % Figure out what points to plot
    beginPoints = -1;
    lastPoints = -1;
    if static == 0
        beginPoints = num_static_roi + 1;
        lastPoints = size(y,2);
        title_name = strcat(title_name, ' for Non-Static Points');
        plotTitle = strcat(plotTitle, '_nonStaticPts');
    elseif static == 1
        beginPoints = 1;
        lastPoints = num_static_roi;
        title_name = strcat(title_name, ' for Static Points');
        plotTitle = strcat(plotTitle, '_staticPts');
    else
        beginPoints = 1;
        lastPoints = size(y,2);
        title_name = strcat(title_name, ' for All Points');
        plotTitle = strcat(plotTitle, '_allPts');
    end
    
    f_plot = figure;
        for i = beginPoints:lastPoints
            % Extract valid displacements and frame numbers
            y_row = y(:,i);
            y_vals = y_row;
            x_vals = x;

            txt = ['ROI',num2str(i)];
            plot(x_vals, y_vals,'DisplayName',txt,'LineWidth',2)
            hold on
        end
        hold off
        legend show
        legend('Location','eastoutside')
        title(title_name)
        grid on
        xlabel(x_axis_name)
        ylabel(y_axis_name)
            set(gcf,'color','w');

    plot_name = fullfile(result_filePath, strcat(imgFolder_name, '_', plotTitle, '.png'));
    saveas(f_plot, plot_name);    % Save plot as .png
    
    fig_name = fullfile(result_filePath, strcat(imgFolder_name, '_', plotTitle, '.fig'));
    savefig(f_plot, fig_name)     % Save plot as .fig
    
end



% Create Axis String Names
function [axis_name, orig_name] = AxisName(var_name)
    orig_name = var_name(7:length(var_name));      % Get rid of 'store_' at the beginning
    x_name_arr = strsplit(orig_name, '_');     % Split string by '_' to extract units (if they exist)
    x_name_str = x_name_arr(1);             % Axis name
    x_units = x_name_arr(length(x_name_arr));   % Possibly axis units
    
    axis_name = x_name_str;           % Set Axis name to axis name
    
    % Check to see if there's units and add them if yes
    if (strcmp(x_name_str, x_units) == 0)
        x_units_str = strcat(' (', x_units, ')');
        axis_name = strcat(x_name_str, x_units_str);   
    end
end
    