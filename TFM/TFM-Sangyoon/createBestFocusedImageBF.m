function MO = createBestFocusedImageBF(pathToOMEFile, seriesIndex, zRange, iChan)
% Function to process OME-TIFF files, select position (series), Z-stack,
% and create MovieObject MO with best-focused images.

% Prompt user to select file if path not provided
if nargin < 1
    [file, path] = uigetfile('*.ome.tif', 'Select an OME-TIFF file.');
    pathToOMEFile = fullfile(path, file);
end

% Load the OME-TIFF file using Bio-Formats
MD = bfImport(pathToOMEFile);

% Validate seriesIndex
if nargin < 2 || isempty(seriesIndex)
    numSeries = numel(MD);
    fprintf('File contains %d series (positions).\n', numSeries);
    seriesIndex = input('Enter the series (position) index to process (default: all): ');
    if isempty(seriesIndex)
        seriesIndex = 1:numSeries;
    end
end

% Validate Z-range
if nargin < 3 || isempty(zRange)
    fprintf('Enter Z-range (start and end Z-slices).\n');
    zStart = input('Start Z-slice: ');
    zEnd = input('End Z-slice: ');
    zRange = zStart:zEnd;
end

% Loop through selected series
MO_list = cell(numel(seriesIndex), 1);
for sIdx = 1:numel(seriesIndex)
    curMD = MD(seriesIndex(sIdx));
    nFrames = curMD.nFrames_;
    nChan = numel(curMD.channels_);
    curPath = fullfile(curMD.outputDirectory_, sprintf('Series_%d', seriesIndex(sIdx)));
    mkdir(curPath);

    % Validate Z-range
    if max(zRange) > curMD.zSize_
        error('Specified zRange exceeds available Z-stack slices for series %d.', seriesIndex(sIdx));
    end

    % Select channel if not provided
    if nargin < 4 || isempty(iChan)
        fprintf('This series contains %d channels.\n', nChan);
        iChan = input('Enter bead channel number (default: final channel): ');
        if isempty(iChan)
            iChan = nChan;
        end
    end

    % Validate channel number
    if iChan > nChan || iChan < 1
        error('Invalid channel number specified for series %d. Available channels: 1-%d.', seriesIndex(sIdx), nChan);
    end

    % Initialize channel directories
    channelDirs = cell(nChan, 1);
    for c = 1:nChan
        channelDirs{c} = fullfile(curPath, sprintf('Channel_%d', c));
        mkdir(channelDirs{c});
    end

    % Process frames
    for t = 1:nFrames
        fprintf('Processing series %d, frame %d of %d.\n', seriesIndex(sIdx), t, nFrames);

        % Process each channel
        for c = 1:nChan
            % Load stack and find best-focused image
            curBeadStack = curMD.channels_(iChan).loadStack(t);
            thresVariance = 0.8; % Default threshold for variance
            applySobel = true;   % Use Sobel filter for focus detection
            averagingRange = findBestFocusFromStack(curBeadStack, thresVariance, applySobel);

            curStack = curMD.channels_(c).loadStack(t);
            selectedStack = curStack(:, :, averagingRange);
            focusedImage = mean(selectedStack, 3); % Calculate best-focused image

            % Save the best-focused image in the channel folder
            imwrite(uint16(focusedImage), ...
                fullfile(channelDirs{c}, sprintf('img_t%03d.tif', t)), ...
                'Compression', 'none');
        end
    end

    % Create MovieData for the processed series
    channels = arrayfun(@(c) Channel(channelDirs{c}), 1:nChan);
    MO = MovieData(channels, curPath);
    MO.setPath(curPath);
    MO.setFilename('movieData.mat');
    MO.sanityCheck();
    MO.save();
    MO_list{sIdx} = MO;
end

% Combine MovieData into MovieList if multiple series are processed
if numel(seriesIndex) > 1
    MO = MovieList([MO_list{:}], fileparts(curPath));
    MO.setPath(fileparts(curPath));
    MO.setFilename('movieList.mat');
    MO.sanityCheck();
    MO.save();
else
    MO = MO_list{1};
end

fprintf('Processing complete. Best-focused images saved in structured folders.\n');
end
