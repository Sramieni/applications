function [outImage, backgroundLevel, binImageEdge] = ptGetProcessedImage (varargin)
% ptGetProcessedImage uses a number of image filtering functions to improve the
% image quality of the phase-contrast cell images so that an optimal segmentation 
% can take place. In particular the functions imVarianceImage, imMinimumThreshold 
% and im2bw are used.
%
% SYNOPSIS       [outImage, backgroundLevel, binImageEdge] = ptGetProcessedImage (varargin)
%
% INPUT          inputImage:     the image that has to be processed
%                greyMax:        the maximum greylevel value of the image (eg 4095 for a 12-bit image)
%                kernelSizeBg:   this should be an odd value and is used to calculate the variance image
%                                used for background subtraction (21 is a good value)
%                kernelSizeEdge: this should be an odd value and is used to
%                                calculate the variance image used for the edge image
%
% OUTPUT         outImage       : the processed image which can be used for segmentation
%                backgroundLevel: the average level of the subtracted background
%                binImageEdge   : image containing edges found in the variance image
%
% DEPENDENCIES   ptGetProcessedImage uses { imVarianceImage
%                                           imMinimumThreshold
%                                           im2bw }
%                                  
%                ptGetProcessedImage is used by { ptTrackCells }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Apr 04          Initial release
% Andre Kerstens        May 04          Function now uses automatic thresholding to find the background
% Andre Kerstens        Jun 04          Added edge image to output of function
% Andre Kerstens        Jul 04          Fixed bug where imMinimumThreshold for the edge images was 
%                                       provided the var image for bg's
% Andre Kerstens        Aug 04          Removed normalization at the end of function
% Andre Kerstens        Aug 04          Fixed bug with input image index calculation

% Test the number of input variables
if nargin < 4
   error('The input image, the maximum grey level (eg 4095 for a 12 bit image) and the kernel sizes for background subtraction and edge image have to be provided. See help ptGetProcessedImage.');
end

% Get the input variables
inputImage   = varargin{1};
greyMax      = varargin{2};
kernelSizeBg = varargin{3};
kernelSizeEdge = varargin{4};
%counter = varargin{5};

% Calculate variance images used for background subtraction and for edge
% detection
varImageBg = imVarianceImage (inputImage, kernelSizeBg);
varImageEdge = imVarianceImage (inputImage, kernelSizeEdge);

% Calculate the optimum thresholds to segment this image (bg subtraction
% and edge thresholds)
[thresholdBg, JBg] = imMinimumThreshold (varImageBg, greyMax);
[thresholdEdge, JEdge] = imMinimumThreshold (varImageEdge, greyMax);

% Get the binary image where all of the background is 1
binImageBg = ~im2bw (varImageBg, thresholdBg);

% Get the binary edge image
binImageEdgeTemp = im2bw (varImageEdge, thresholdEdge);

% Make the edge pixels 0 again using a zero mask
mask = zeros (size (binImageEdgeTemp));
mask(kernelSizeEdge+1:end-kernelSizeEdge, kernelSizeEdge+1:end-kernelSizeEdge) = ...
    binImageEdgeTemp(kernelSizeEdge+1:end-kernelSizeEdge, kernelSizeEdge+1:end-kernelSizeEdge);
binImageEdge = mask;

% Show only the background pixels in the original image (by multiplying with 
% the binary image). And then fetch the index of these specific pixels
newImage = inputImage .* binImageBg;
backOnlyImage = newImage (find (binImageBg));

% Get the x and y coordinates of these pixels
[x,y] = find (binImageBg);

% Now we start to solve the equation I(x,y) = f(x,y) where f(x,y) is ax^2+bxy+cy^2+d
% For this we need to create a matrix A = [x^2, xy, y^2, 1] to solve the coeff. 
% a, b, c and d (1 can be generated by ones(size(x),1)). The equation can then be solved
% by A \ backOnlyImage
one = ones(size(x,1),1);
A = [x.^2, x.*y, y.^2, one];
coeff = A \ backOnlyImage;

% Now calculate the background image using the estimated coefficients (we need
% the size of the whole input image this time so that we can subtract later).
defImage = ones (size(inputImage,1), size(inputImage,2));
[xi,yi] = find (defImage);
backgroundImage = coeff(1).*(xi.^2) + coeff(2).*(xi.*yi) + coeff(3).*(yi.^2) + coeff(4);
bgImage = reshape(backgroundImage,size(inputImage));

%filename = ['bgImage_' num2str(counter) '.mat'];
%fprintf (1, 'filename = %s\n', filename);
%save (filename, 'bgImage');

% The average background level is approx. equal to the constant in the equation
backgroundLevel = coeff(4);

% Now we only have to subtract the background from the image and we're done
outImage = inputImage - bgImage;

% Let's normalize the image back to [0..1] again
%imageMinimum = min (min (outImage));
%imageMaximum = max (max (outImage));
%outImage = (outImage - imageMinimum) / (imageMaximum - imageMinimum);

% The background level will be approx. 0 in the subtracted image
backgroundLevel = 0;
