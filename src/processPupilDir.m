function [arrayeyes, error] = processPupilDir(path)
%PROCESSPUPILDIR Directory analysis and calling the image-processing function
% 
%   Input:
%       path                       absolute path of the directory to be analyzed
%
%   Output:
%       arrayeyes                  array of "eye"-typed structures
%       error                      vector with the calculation error of the pupil center
%
% 
% 
%   Author: jlnkls
%
%   26/01/2016


%% File extension definition
image.Extension = 'png';
data.Extension = 'mat';

%% Extraction of image filenames
image.Files = dir( fullfile(path,['*.' image.Extension]) );
[~, image.Files] = cellfun(@fileparts, {image.Files.name}, 'UniformOutput', false);

image.Number = length(image.Files);                                                      % Número de archivos a analizar

%% Extraction of filename with "ground-truth" centers
data.File = dir( fullfile(path,['*.' data.Extension]) );
[~, data.File] = cellfun(@fileparts, {data.File.name}, 'UniformOutput', false);

data.File = sprintf(data.File{:});
data.File = strcat(data.File, '.', data.Extension);

%% Loading filename with "ground-truth" centers
load(strcat(path, data.File));

%% Information initialization of images to be processed and error
arrayeyes = struct('center', cell(1,image.Number), ...
    'error', cell(1,image.Number), 'contour', cell(1,image.Number), ...
    'filename', cell(1,image.Number));

error = NaN(1,image.Number);

%% Processing images one by one
for index = 1:(image.Number)
    
    image.Name = image.Files(index);
    image.Name = sprintf(image.Name{:});
    image.Name = strcat(image.Name, '.', image.Extension);
    
    image.Current = imread(strcat(path, image.Name));

    arrayeyes(1,index) = processImage(image.Current,index,[arrayeyes.center]);
    arrayeyes(1,index).filename = image.Name;
   
    error(1,index) = sqrt((sum(arrayeyes(index).center - centers(index,:))) .^ 2);
        
end



end



