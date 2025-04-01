function [] = erodeForceField(MD, band)
% erodeForceField overlays eroded force fields and computes strain energy.
%
% INPUT:
%   MD   : MovieData object or path
%   band : erosion width in microns (default = 4)
%
% OUTPUT:
%   Saves overlay images and strainEnergy.mat in /heatmap

    if nargin < 2 || isempty(band)
        band = 1;
    end

    % Load MovieData if given as path
    if ischar(MD)
        load(MD, 'MD');
    end

    % Get TFM package
    pkgIdx = MD.getPackageIndex('TFMPackage');
    assert(pkgIdx > 0, 'TFMPackage not found.');
    tfmPkg = MD.getPackage(pkgIdx);
    ffProc = tfmPkg.getProcess(4); % Force Field
    displProc = tfmPkg.getProcess(3); % Displacement

    % Load shifted, filtered force field
    sFF = load(ffProc.outFilePaths_{1});
    assert(isfield(sFF, 'forceFieldShifted'), 'Missing forceFieldShifted. Run filterMovieForceField first.');
    forceField = sFF.forceFieldShifted;

    % Load displacement field directly from .mat file
    sDisp = load(displProc.outFilePaths_{1});
    displFieldRaw = sDisp.displField;  % Use displFieldShifted if needed

    % Load ROI mask
    roiMaskPath = fullfile(MD.outputDirectory_, 'roiMask.tif');
    assert(isfile(roiMaskPath), 'ROI mask not found. Run changeROIMask.m');
    roiMask = imread(roiMaskPath) > 0;

    % Resize mask if needed
    imgSize = [MD.imSize_(2), MD.imSize_(1)];
    if any(size(roiMask) ~= imgSize)
        roiMask = imresize(roiMask, imgSize, 'nearest');
    end

    % Apply erosion
    pixelSize = MD.pixelSize_; %conversts the microns to pixels 
    erosionRadius = max(1, round(band / pixelSize));
    erodedMask = imerode(roiMask, strel('disk', erosionRadius));
    fprintf('Erosion applied with radius = %d pixels.\n', erosionRadius);

    % Output folder
    outputDir = fullfile(MD.outputDirectory_, 'heatmap');
    if ~exist(outputDir, 'dir'); mkdir(outputDir); end
    nFrames = MD.nFrames_;
    strainEnergy = zeros(nFrames, 1);
    cellChan = min(2, numel(MD.channels_));

    % Loop over frames
    for i = 1:nFrames
        if i > numel(forceField)
            warning('Frame %d skipped: forceField missing.', i);
            continue;
        end
        ff = forceField(i);
        if isempty(ff.pos) || isempty(ff.vec)
            fprintf('Frame %d skipped: no force vectors.\n', i);
            continue;
        end

        % Get displacement for frame i
        try
            if iscell(displFieldRaw)
                df = displFieldRaw{i};
            elseif isstruct(displFieldRaw) && numel(displFieldRaw) >= i
                df = displFieldRaw(i);
            else
                error('Unknown displacement field format.');
            end
        catch
            warning('Displacement data missing for frame %d. Skipping.', i);
            continue;
        end

        % Extract vector data
        x = ff.pos(:,1); y = ff.pos(:,2);
        fx = ff.vec(:,1); fy = ff.vec(:,2);
        ux = df.vec(:,1); uy = df.vec(:,2);

        % Convert to pixels
        xPix = round(x / pixelSize);
        yPix = round(y / pixelSize);
        valid = xPix > 0 & xPix <= imgSize(2) & yPix > 0 & yPix <= imgSize(1);

        % Mask filtering
        maskIdx = sub2ind(size(erodedMask), yPix(valid), xPix(valid));
        inMask = erodedMask(maskIdx);

        x = x(valid); y = y(valid);
        fx = fx(valid); fy = fy(valid);
        ux = ux(valid); uy = uy(valid);

        x = x(inMask); y = y(inMask);
        fx = fx(inMask); fy = fy(inMask);
        ux = ux(inMask); uy = uy(inMask);

        % Strain energy
        strainEnergy(i) = 0.5 * sum(fx .* ux + fy .* uy) * 1e-15;
        fprintf('Frame %d: %d vectors in ROI, SE = %.2f fJ\n', i, numel(x), strainEnergy(i));

        % Load background image
        img = MD.channels_(cellChan).loadImage(i);
        imgNorm = mat2gray(double(img));

        % Plot overlay
        fig = figure('Visible', 'off');
        imshow(imgNorm); hold on;

        if ~isempty(fx)
            mags = sqrt(fx.^2 + fy.^2);
            cmap = jet(256);
            normMags = round(rescale(mags, 1, 256));
            colors = cmap(normMags, :);
            scale = 5;

            for k = 1:length(x)
                quiver(x(k)/pixelSize, y(k)/pixelSize, fx(k)*scale, fy(k)*scale, ...
                       'Color', colors(k,:), 'LineWidth', 1);
            end
        end

        title(sprintf('Frame %d - Eroded Force Vectors', i));
        saveas(fig, fullfile(outputDir, sprintf('overlay_frame%03d.tif', i)));
        close(fig);
    end

    % Save strain energy result
    save(fullfile(outputDir, 'strainEnergy.mat'), 'strainEnergy');
    fprintf('? Eroded overlays and strain energy saved in: %s\n', outputDir);
end
