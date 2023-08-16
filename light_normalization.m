%% Perform light normalization on an entire image
% @author: Casey Rodgers
% @ date:  04/28/2022

% Image folder and image filenames
root_foldername = "imgs\";
img_filenames = ["TotalFull", "TotalFullDarker", "TotalFullShadow"];
img_type = ".jpg";  % Image type

% Loop through images
for i = 1:length(img_filenames)
    % Read in image
    img_filepath = strcat(root_foldername, img_filenames(i), img_type);  % Put together image filepath
    orig_img = imread(img_filepath);    % Read in image
    figure
    imshow(orig_img);  % Show original image
    
    % Light Normalization 
    new_img = histeq(orig_img);
    figure
    imshow(new_img);

    % Save new image
    imwrite(new_img, strcat(root_foldername, img_filenames(i), "_lightNormal", img_type));

end