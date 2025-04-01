% Step 1: Select the OME-TIFF file
[file, path] = uigetfile('*.ome.tif', 'Select an OME-TIFF file.');
if isequal(file, 0)
    error('No file selected. Exiting script.');
end

pathToOMEFile = fullfile(path, file);

% Step 2: Load the OME-TIFF file using Bio-Formats
MD = bfImport(pathToOMEFile);

% Step 3: User selects series (time-lapse frame)
numSeries = numel(MD);
fprintf('Number of time-lapse frames (series) available: %d\n', numSeries);

for i = 1:numSeries
    fprintf('Series %d corresponds to Time-lapse Frame %d\n', i, i);
end

seriesIndex = input(sprintf('Enter the series number (1-%d): ', numSeries));
curMD = MD(seriesIndex);

% Step 4: Extract Image Properties
sizeX = curMD.imSize_(1);
sizeY = curMD.imSize_(2);
sizeZ = curMD.zSize_;
sizeC = numel(curMD.channels_);

% Set sizeT to 1 manually
sizeT = 1;
fprintf('Forcing sizeT = 1 (time points are separate series).\n');

% Step 5: User selects actin channel
fprintf('This series contains %d channels.\n', sizeC);
channelIndex = input('Enter the channel number corresponding to actin: ');

% Step 6: Load the Full Z-Stack at Once (Fixing `t` Issue)
fprintf('Loading full Z-stack for series %d, channel %d...\n', seriesIndex, channelIndex);
tIndex = 1;  % Since time is not used, set it to 1
curStack = curMD.channels_(channelIndex).loadStack(tIndex); % FIX: Provide `tIndex`

% Check if stack is empty
if isempty(curStack)
    error('Failed to load Z-stack for the selected series and channel.');
end

% Step 7: Perform Maximum Intensity Projection (MIP)
mipImage = max(curStack, [], 3);
mipImage = mat2gray(mipImage);

% Enhance actin structures
sensitivity = 0.4;
T = adaptthresh(mipImage, sensitivity, 'ForegroundPolarity', 'bright');
enhancedActin = imbinarize(mipImage, T);

% Display results
figure('Name', 'Actin Enhancement', 'NumberTitle', 'off');
subplot(1,2,1); imshow(mipImage, []); title('Original MIP');
subplot(1,2,2); imshow(enhancedActin, []); title('Enhanced Actin');

% Step 8: Save the Processed Image
outputName = input('Enter output filename (without extension): ', 's');
outputPath = fullfile(path, [outputName, '_MIP.tif']);
imwrite(mipImage, outputPath);
disp(['MIP image saved at: ', outputPath]);
