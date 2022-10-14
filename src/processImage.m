function [ eye ] = processImage(originalImage, nImage, centers)
%processImage Obtaining the center of a pupil by processing
%             image
%
%   Input:
%       originalImage             original eye image to be processed
%       nImage                    number of the image to be processed
%       centers                   centers from past images
%   
%   Output:
%       eye                       structure that contains information on the center and 
%                                 pupil outline
%
%   Function calls:
%      regionGrowing
%      ILPF 
% 
%   Author: jlnkls
%
%   26/01/2016


%% Image processing

% Transformation of the 'centers' vector into an Nx2 matrix
centers = vec2mat(centers,2);

% Conversion of original image to grayscale and numeric type 'double'
image.gray = double(rgb2gray(originalImage));

% Median filtering
image.medianFilter = medfilt2(image.gray);

% Image normalization (range: [0-1])
image.medianFilter = image.medianFilter/max(image.medianFilter(:));

% Parameter definition of the disc-shaped structuring element with a radius of 5 pixels
strElement.radius = 5;
strElement.approximation = 0;
strElement.disk = strel('disk', strElement.radius, ...
    strElement.approximation);

% Erosion with the previously defined structuring element
image.erode = imerode(image.medianFilter,strElement.disk);

% Definition of low pass filter parameters
lowPassFilter.size = [10, 10];
lowPassFilter.cutOff = 5;
lowPassFilter.mask = ILPF( lowPassFilter.size, lowPassFilter.cutOff);

% LPF
image.lowPassFilter_1 = imfilter(image.erode, lowPassFilter.mask);

% Definition of low pass filter parameters
lowPassFilter.size = [5, 5];
lowPassFilter.cutOff = 15;
lowPassFilter.mask = ILPF( lowPassFilter.size, lowPassFilter.cutOff);

% LPF
image.lowPassFilter_2 = imfilter(image.lowPassFilter_1, ...
    lowPassFilter.mask);

% Morphological closure with structuring element of radius 5
image.close = imclose(image.lowPassFilter_2,strElement.disk);

% LPF
image.lowPassFilter_3 = imfilter(image.close,lowPassFilter.mask);

% Image normalization (range: [0-1])
image.lowPassFilter_3 = (image.lowPassFilter_3) / ...
    (max(image.lowPassFilter_3(:)));

%% Image thresholding and segmentation

% Obtaining dimensional values of the image
[image.rows, image.columns, ~] = size(originalImage);

% Cropping of image margins, avoiding the appearance of 
% spurious values (elimination of the first and last 10 rows and columns)
image.lowPassFilter_3(1:end,1:10) = 1;
image.lowPassFilter_3(1:end,((image.columns-10):image.columns)) = 1;
image.lowPassFilter_3(1:10,1:end) = 1;
image.lowPassFilter_3(((image.rows-10):image.rows),1:end) = 1;

% Thresholding with level set to 40% of the maximum value ("blank")
image.thresh = image.lowPassFilter_3;
image.thresh(image.lowPassFilter_3<=0.4) = 1;
image.thresh(image.lowPassFilter_3>0.4) = 0;

%% Eye-tracking pupil segmentation

% Circular Hough Transform (CHT)
[houghTransform.centers, houghTransform.radii, houghTransform.metric] = ...
    imfindcircles(image.thresh,[10, 15]);

% If the result of the CHT does not find any circle,
% circular regions will be searched by alternative methods
if (isempty(houghTransform.centers))
    
    % Obtaining statistics of the thresholded regions
    statistics = regionprops(image.thresh, 'all');
     
    % Calculation of the circularity of the connected regions found
    Regions.perimeter = [statistics.Perimeter];
    Regions.area = [statistics.Area];
    Regions.circularity = (Regions.perimeter .^ 2) ./ (4 * pi * Regions.area);  

    % If the previous pre-processing does not allow to find the pupil, 
    % a processing will be carried out by means of 'top-hat'
    if isempty(statistics)
    
        % Definition of a structuring element 'disc' of radius 400        
        strElement.disk = strel('disk', 400);
        image.topHat_1 = image.lowPassFilter_3;
        
        % Top-Hat (histogram smoothing)
        image.topHat_1 = imtophat(image.topHat_1, strElement.disk);
        [~, topHat.indice] = min(image.topHat_1(:));
        
        % Growth of region, from seed located in the 
        % darkest region of the image after Top-Hat
        [topHat.rows, topHat.columns] = ...
            ind2sub(size(image.topHat_1),topHat.indice);
        [~, image.topHat_2] = ...
            regionGrowing(image.topHat_1, [topHat.rows, topHat.columns]);
        
        % Obtaining statistics of the processed image        
        statistics = regionprops(image.topHat_2, 'all');           

        % Calculation of the circularity of the connected regions found        
        Regions.perimeter = [statistics.Perimeter];
        Regions.area = [statistics.Area];
        Regions.circularity = (Regions.perimeter  .^ 2) ./ ...
            (4 * pi * Regions.area);
    
        % If the pupil is still not found, 
        % the image will be partitioned by Otsus method
        if isempty(statistics)          

            % Threshold of the image by partitions of 10x10 pixels, 
            % with 5x5 pixels of overlap
            Otsu.function = @(block_struct) ...
                im2bw(block_struct.data,graythresh(block_struct.data));
            image.otsu = blockproc(image.topHat_2, [10 10], ...
                Otsu.function, 'BorderSize', [5 5]);

            % Image dilation (union of possible disconnected pixels after block processing)
            image.otsu = imdilate(image.otsu,strel('square',2));
            % Image erosion (shortening of enlarged elements after dilation)
            image.otsu = imerode(image.otsu,strel('square',1));
            % Image inversion
            image.otsu = logical(1-image.otsu);

            % Obtaining statistics of the thresholded regions            
            statistics = regionprops(image.otsu, 'all');     

            % Calculation of the circularity of the connected regions found           
            Regions.perimeter = [statistics.Perimeter];
            Regions.area = [statistics.Area];
            Regions.circularity = (Regions.perimeter  .^ 2) ./ ...
                (4 * pi * Regions.area);

        end  
    
    end
     
    % In the event that the statistics return values, 
    % and that at least 5 images have been analyzed, 
    % it will be checked that the information is correct, 
    % comparing it with the previous position of the pupil
    if ((nImage>5) && (~isempty(statistics)))
        
        % In 'nEyeCalculated', the difference between 
        % the position of the current image and 
        % the position of the last image in which 
        % a center has been obtained will be obtained     
        nEyeCalculated = 1;
            
        for i=((nImage-1):-1:(nImage-5))        
                   
            if (~isnan(centers(i,1)))
                    nEyeCalculated = nImage-i;
                    break
            end
        end
            
            
        % In the event that several regions are obtained, the one closest to 
        % the point obtained in the previous image will be selected
        if ((length(statistics)>1) && (any(~isnan(centers(nImage-5:nImage-1,1)))))
                    
                       
            % Iterating over all the regions, 
            % and calculation of the distance from the center of the previous image
            clear distance;
            
            for i = 1:length(statistics)
                Regions.Centroid = [statistics.Centroid];
                Regions.distance(i) = ...
                    sqrt((statistics(i).Centroid(1) - ...
                    centers(nImage-nEyeCalculated,1))^2 + ...
                    (statistics(i).Centroid(2) - ...
                    centers(nImage-nEyeCalculated,1))^2);
            end            
            
            % The closest center will be chosen as the center of the pupil
            [~, Regions.index] = nanmin(Regions.distance);
            
        else
        
            % If there is no information about the previous image, 
            % the region with circularity closest to unity will be selected
            [~, Regions.index] = min(abs(Regions.circularity-1));    
        
        end
           
        % If the distance between the calculated center and 
        % the previous center is greater than 30 pixels, 
        % the center will be discarded.
        % If a center has not been obtained in the previous image, 
        % the distance with the last image obtained will be calculated   
        Regions.maxDistance = 30;
        Regions.distanceNImage = ...
            sqrt((statistics(Regions.index).Centroid(1) - ...
            centers(nImage-nEyeCalculated,1))^2 + ...
            (statistics(Regions.index).Centroid(2) - ...
            centers(nImage-nEyeCalculated,2))^2);
        
        % Discarding the center
        if (Regions.maxDistance<Regions.distanceNImage) 
            
            Regions.index=[];     

        end
        
    else
        
        % Choice of the center of the region with circularity closest to the unit
        [~, Regions.index] = min(abs(Regions.circularity-1));   
      
    end
       
%% Assigning output parameters
    
    % No pupil detected --> no center
    if isempty(Regions.index)
    
        eye.center(1,1:2) = [NaN, NaN];     
        eye.error = true;
        eye.contour = [];
    
    % There is a pupil detected --> a center exists
    else
    
        % Center of the pupil (centroid)
        eye.center(1,1:2) = statistics(Regions.index).Centroid;
        eye.error = false;

        % Pupil outline (region outline)
        statistics.contour = zeros(image.rows, image.columns);
        statistics.contour(statistics.PixelIdxList) = 1;
        statistics.contour = edge(statistics.contour);
        [statistics.contour_R ,statistics.contour_C]  = ...
            find(statistics.contour == 1);
        eye.contour = [statistics.contour_R, statistics.contour_C];

    end

else
       
    % Hough transform
    eye.center(1,1:2) = houghTransform.centers;
    eye.error = false;
    
    % Hough Transform Contour Calculation
    houghTransform.th = 0:pi/50:2*pi;
    houghTransform.xunit = (houghTransform.radii * ...
        cos(houghTransform.th)) + houghTransform.centers(1,1);
    houghTransform.yunit = (houghTransform.radii * ...
        sin(houghTransform.th)) + houghTransform.centers(1,2);
    
    eye.contour = [houghTransform.xunit; houghTransform.yunit]';


    
end
   
% Name of the processed file (assigned afterwards)
eye.filename = [];

end

