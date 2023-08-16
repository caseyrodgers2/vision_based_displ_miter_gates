%% Create Video from a folder of images
% @author: Casey Rodgers
% @ date:  06/09/2021
% version : V1

clear all
close all
clc

% Get images
imgFolder = 'imgs\';   % Image folder (Folder where the images are)
filePattern = fullfile(imgFolder, '*.jpg');
imgFiles = dir(filePattern);
imgFolder_name = imgFolder(1:length(imgFolder)-1);  % Folder name

% Create VideoWriter object to make video
video_name = fullfile(imgFolder, strcat(imgFolder_name, '_video_file.mp4'));
v = VideoWriter(video_name, 'MPEG-4');
v.Quality = 80;
v.FrameRate = 15;
open(v)

% Loop through images
for k = 1 : size(imgFiles,1)
%for k = 1 : 2
    % Load image
    img1 = fullfile(imgFolder_name, imgFiles(k).name);
    img_arr = imread(char(img1));
    img_arr2 = imresize(img_arr,[1080, NaN]);
    img_arr2_hist = histeq(img_arr2);
    text_str = num2str(k);
    img_arr3 = insertText(img_arr2_hist, [size(img_arr2_hist,2)- 150,25] ,text_str,'FontSize',25,'BoxColor','yellow');   
%     figure,
%         imshow(img_arr3)
    % Add image to video
    writeVideo(v, img_arr3)
    disp(k)
%     if (mod(k, 5) == 0)
%         disp(k)
%     end
    
%     pause(0.5);    
end

% Close VideoWriter
close(v)
