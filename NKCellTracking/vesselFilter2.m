function [BW,maskedImage] = segmentImage(X)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 30-Mar-2023
%----------------------------------------------------

X= imgaussfilt(X);
% Adjust data to span data range.
X = imadjust(X);

% Threshold image - adaptive threshold
BW = imbinarize(X, 'adaptive', 'Sensitivity', 0.75, 'ForegroundPolarity', 'bright');

% Open mask with square
width = 5;
se = strel('disk', width);
BW = imopen(BW, se);

% Dilate mask with square
width = 4;
se = strel('disk', width);
BW = imdilate(BW, se);

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

