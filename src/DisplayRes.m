function [] =  DisplayRes(path)
%DISPLAYRES Image display function, with true centers and
%           overlapping centers detected
% 
%   Input:
%       path                       absolute path of the directory to be analyzed
%
%   Output:
%       ~
%
% 
% 
%   Author: jlnkls
%
%   26/01/2016


%% Image analysis and processing
[arrayeyes, error] = processPupilDir(path);

%% File extension definition
data.Extension = 'mat';

%% Extraction of filename with "ground-truth" centers
image.Number = length(arrayeyes);

data.File = dir( fullfile(path,['*.' data.Extension]) );
[~, data.File] = cellfun(@fileparts, {data.File.name}, 'UniformOutput', false);

data.File = sprintf(data.File{:});
data.File = strcat(data.File, '.', data.Extension);

%% Loading filename with "ground-truth" centers
load(strcat(path, data.File));

%% Processing images one by one
for index = 1:(image.Number)

    image.Current = imread(strcat(path, arrayeyes(index).filename));
    
    figure(1);
    imshow(image.Current);
    hold on;
    
    image.IDSplit = strsplit(arrayeyes(index).filename, '-');
    image.ID = strcat(image.IDSplit(2), '-' , image.IDSplit(3));

    text(0.02, 1-0.1, image.ID, ...
        'BackgroundColor',[96/255 159/255 96/255],'Color','White', ...
        'FontSize',14,'Units','normalized');
    text(1-0.02, 1-0.1, strcat('Error: ', ...
         num2str(error(index),'%.1f')), ... 
        'BackgroundColor', [163/255 99/255 99/255], ...
        'HorizontalAlignment', 'right', 'Color','White', ...
        'FontSize',14,'Units','normalized');
    
    plot(arrayeyes(index).center(1),arrayeyes(index).center(2), ...
        'Marker','.','Color','cyan');
    plot(centers(index,1), centers(index,2), 'Marker', '+', 'Color', 'r');
      
    
    pause;
    
    clf;
    
    
    
end


close all;


end