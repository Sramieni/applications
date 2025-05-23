function projectDynROI(MD,varargin)
% WIPS: 
%   -   respect computeMIPProcess Specs
%   -   using dynROI class instead of ad hoc object
%   -   comment and clarify options.
%   -   shorter function 
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('MD'); 
ip.addOptional('insetDynROI',[]);     % Every pixel mapped inside the ROI is rendered in a MIP mask
ip.addOptional('dynROI',[]);          % This dynROI define the boundary of the current view (default, insetDynROI) only used if 'crop' is set to manifold
ip.addOptional('renderedChannel',1:length(MD.channels_));
ip.addOptional('crop','manifold');  % 'manifold': crop around p.dynROI, 'full': show the full volume.
ip.addOptional('FoF',[]);           % frame of refenrence for the projection
ip.addOptional('transType','affineOnePass');
ip.addOptional('processFrame',1:MD.nFrames_);
ip.addOptional('channelRender','grayRed');
ip.addOptional('intMinPrctil',[1 99.9]);
ip.addOptional('intMaxPrctil',[100 100]);
ip.addOptional('name',[]);
ip.addOptional('processSingleProj',[]);
ip.addOptional('suppressROIBorder',false);
ip.addOptional('format','png');
ip.addOptional('rawTIFF',false);
ip.addOptional('processMaskVolume',[])
ip.addOptional('processRenderer',[]);
ip.addOptional('fringeWidth',[]);
ip.addOptional('insetFringeWidth',20);
ip.addOptional('maxMIPSize',max([400,MD.imSize_,ceil(MD.zSize_*MD.pixelSizeZ_/MD.pixelSize_)]));
ip.parse(MD,varargin{:});
p=ip.Results;

%% Fill gaps in dynROI
insetDynROI=p.insetDynROI;
dynROI=p.dynROI;
fillTrackGaps(insetDynROI);
% fillTrackGaps(dynROI);

if(isempty(dynROI))
    if(~isempty(p.FoF))
        dynROI=p.FoF.applyBase(insetDynROI,'');
    else
        dynROI=insetDynROI;
    end
end

rawTIFF=p.rawTIFF;

%% Define size of mapping areas for insetDynROI and dynROI.
insetFringeWidth=p.insetFringeWidth;
if(isempty(p.fringeWidth))
    fringeWidth=insetFringeWidth;
else
    fringeWidth=p.fringeWidth;
end
processFrame=p.processFrame;

%% Set normalization value
minIntensityNorm=zeros(1,numel(MD.channels_));
maxIntensityNorm=zeros(1,numel(MD.channels_));
for chIdx=1:length(MD.channels_)
    vol=MD.getChannel(chIdx).loadStack(1); 
    minIntensityNorm(chIdx)=[ prctile(vol(:),p.intMinPrctil(chIdx))];
    maxIntensityNorm(chIdx)=[ prctile(vol(:),p.intMaxPrctil(chIdx))];
end

%% Define the static Rectangular cuboid that contains the pixel to be projected in the frame of reference.
%% Accordinly,  the coordinate of this cube are specified such as the origin of the frame of reference is the zero.

%% In the manifold crop case, the boundaries are given by the transform coordinate along the manifold polygon
%% in the full case, one have to estimate the maximum rectangle cuboid contained that can descibe the extremum coordinate of the original volume
if(strcmp(p.crop,'manifold')&&(~isempty(dynROI)))
    
    minX=MD.getDimensions('X')+1;
    minY=MD.getDimensions('Y')+1;
    minZ=(MD.pixelSizeZ_/MD.pixelSize_)*MD.getDimensions('Z')+1;
    maxX=0; 
    maxY=0;
    maxZ=0;
    for iP=1:length(dynROI)
        minX=floor(min(min(dynROI(iP).x),minX));
        minY=floor(min(min(dynROI(iP).y),minY));
        minZ=floor(min(min(dynROI(iP).z),minZ));
        maxX=ceil(max(max(dynROI(iP).x),maxX));
        maxY=ceil(max(max(dynROI(iP).y),maxY));
        maxZ=ceil(max(max(dynROI(iP).z),maxZ));
    end
    maxXBorder=(maxX+fringeWidth);
    maxYBorder=(maxY+fringeWidth);
    maxZBorder=(maxZ+fringeWidth);
    minXBorder=(minX-fringeWidth);
    minYBorder=(minY-fringeWidth);
    minZBorder=(minZ-fringeWidth);
else
    if(~isempty(p.FoF))
        maxXBorder=MD.getDimensions('X')-p.FoF.origin(1,1);
        maxYBorder=MD.getDimensions('Y')-p.FoF.origin(1,2);
        maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_)-p.FoF.origin(1,3);
        minXBorder=1-p.FoF.origin(1,1);
        minYBorder=1-p.FoF.origin(1,2);
        minZBorder=1-p.FoF.origin(1,3);
    else
        maxXBorder=MD.getDimensions('X');
        maxYBorder=MD.getDimensions('Y');
        maxZBorder=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
        minXBorder=1;
        minYBorder=1;
        minZBorder=1;
    end
end

tform=affine3d();

if(~isempty(p.FoF))
    B=p.FoF.getBase(1);
    tform.T(1:3,1:3)=B;
end
format=p.format;
if(rawTIFF)
    format='tif';
end

projectDynROIProcess=p.processSingleProj;
set(projectDynROIProcess,'ref',p.FoF);
set(projectDynROIProcess,'nFrames',length(p.processFrame));
projectDynROIProcess.setBoundingBox( ...
   [minXBorder maxXBorder],...
   [minYBorder maxYBorder],...
   [minZBorder maxZBorder] );

if(~isempty(p.processSingleProj))
    %% simulate run to comply with movieViewer requirement
    p.processSingleProj.setFunName((@(x) x));
    p.processSingleProj.run();
    %% Save the current function run for futur rerun
    p.processSingleProj.setFunName((@(p) projectDynROI(MD,varargin{:},'processSingleProj',p)));
end


format='png';
% Standardized output for processRenderer
outFilePathsRenderer = cell(1, 5);

if(~isempty(p.processRenderer))
    set(p.processRenderer,'ref',p.FoF);
    set(p.processRenderer,'nFrames',length(p.processFrame));
    p.processRenderer.setBoundingBox( ...
       [minXBorder maxXBorder],...
       [minYBorder maxYBorder],...
       [minZBorder maxZBorder] );
    %% simulate run to comply with movieViewer requirement
    p.processRenderer.setFunName((@(x) x));
    p.processRenderer.run();
    %% Save the current function run for futur rerun
    p.processRenderer.setFunName((@(p) projectDynROI(MD,varargin{:},'processRenderer',p)));
end


% Optional process for saving sparse volume
if(~isempty(p.processMaskVolume))
    outFilePathsMaskVolume=cell(2,numel(MD.channels_));
    frameNb=length(p.processFrame);
    for i = 1:numel(MD.channels_);    
        outFilePathsMaskVolume{1,i} = [outputDirSingleProj filesep 'volMask'  filesep 'ch' filesep 'frame_nb%04d.mat'];
        outFilePathsMaskVolume{2,1} = [outputDirSingleProj filesep 'volMask' filesep 'limits.mat'];
    end;
    for mIdx=1:numel(p.processFrame)
        mkdirRobust(fileparts(outFilePathsMaskVolume{1,i}));
    end
    save(outFilePathsMaskVolume{2,1},'minXBorder', 'maxXBorder','minYBorder','maxYBorder','minZBorder','maxZBorder','frameNb');
    p.processMaskVolume.setOutFilePaths(outFilePathsMaskVolume);
end

% warp, crop, fuse and save each time point
% vol = cell(processFrame,1);
for fIdx = processFrame
    fprintf('.') 
    % produce a ROI mask using the 1D polygon (segment defined by the extremities of the insetDynROI).
    % todo: N Channel (now 2).
    mask=ones(MD.imSize_(2),MD.imSize_(1),MD.zSize_);
    maskedVol=[];
    if(~isempty(insetDynROI))
         % Collect relative frameIdx
        pIndices=nan(1,length(insetDynROI));
        for polIdx=1:length(insetDynROI)
            F=insetDynROI(polIdx).f;
            pIdx=find(F==fIdx,1);
            if isempty(pIdx)
                if(fIdx>max(F))   pIdx=length(F);  else   pIdx=1; end;
            end
            pIndices(polIdx)=pIdx;
        end

        %% Building mask in the 1D case
        nextPoint=length(insetDynROI);
        PCurrent=[insetDynROI(1).x(pIndices(1)) insetDynROI(1).y(pIndices(1)) insetDynROI(1).z(pIndices(1))];
        KCurrent=[insetDynROI(nextPoint).x(pIndices(nextPoint)) insetDynROI(nextPoint).y(pIndices(nextPoint)) insetDynROI(nextPoint).z(pIndices(nextPoint))];

        vol = MD.getChannel(1).loadStack(fIdx);

        % Building mask for both channel on the whole volume
        % NOTE: in order to apply fringe isotropically, we need the mask to
        % be isotropized briefly.
        mask=zeros(size(vol,1),size(vol,2),ceil(size(vol,3)*MD.pixelSizeZ_/MD.pixelSize_));
        sampling=100;
        xSeg=round(linspace(PCurrent(1),KCurrent(1),sampling));
        ySeg=round(linspace(PCurrent(2),KCurrent(2),sampling));
        zSeg=round(linspace(PCurrent(3),KCurrent(3),sampling));
        indx=sub2ind(size(mask),ySeg,xSeg,zSeg);
        
        mask(indx)=1;
        
        %mask=imdilate(mask,IMSphere);  %ones(cubeHalfWidth,cubeHalfWidth,round(cubeHalfWidth*MD.pixelSize_/MD.pixelSizeZ_)));
        %% If no transform are needed, now to save on bwdist.

        distMap=mask;
        distMap=bwdist(distMap);
        mask(distMap<(insetFringeWidth+1))=1;
        [y x z]=...
            ndgrid( linspace(1,size(mask,1),size(vol,1)),...
                    linspace(1,size(mask,2),size(vol,2)),...
                    linspace(1,size(mask,3),size(vol,3)));
        mask=interp3(mask,x,y,z);
    end
    mips=cell(3,length(MD.channels_));
    if(~isempty(insetDynROI))
        if(isempty(p.FoF))&&(strcmp(p.crop,'manifold'))
            aminXBorder=max(1,minXBorder);
            aminYBorder=max(1,minYBorder);
            aminZBorder=max(1,minZBorder);
            amaxXBorder=min(size(mask,2),maxXBorder);
            amaxYBorder=min(size(mask,1),maxYBorder);
            amaxZBorder=min(size(mask,3)*MD.pixelSizeZ_/MD.pixelSize_,maxZBorder);
            
            XCropMask=false(1,size(mask,2));
            XCropMask(aminXBorder:amaxXBorder)=true;
            
            YCropMask=false(1,size(mask,1));
            YCropMask(aminYBorder:amaxYBorder)=true;
            
            ZCropMask=false(1,size(mask,3));
            ZCropMask(ceil((aminZBorder:amaxZBorder)*MD.pixelSize_/MD.pixelSizeZ_))=true;
            
            ZCropVol=false(1,size(vol,3));
            ZCropVol(ceil((aminZBorder:amaxZBorder)*MD.pixelSize_/MD.pixelSizeZ_))=true;
            
            mask(:,:,~ZCropMask)=[];
            mask(~YCropMask,:,:)=[];
            mask(:,~XCropMask,:)=[];
        end
    end
    
    for chIdx=1:length(MD.channels_)
        vol=MD.getChannel(chIdx).loadStack(fIdx);
        if(~isempty(insetDynROI))
            if(isempty(p.FoF))&&(strcmp(p.crop,'manifold'))
                aminXBorder=max(1,minXBorder);
                aminYBorder=max(1,minYBorder);
                aminZBorder=max(1,minZBorder);
                amaxXBorder=min(size(vol,2),maxXBorder);
                amaxYBorder=min(size(vol,1),maxYBorder);
                amaxZBorder=min(size(vol,3)*MD.pixelSizeZ_/MD.pixelSize_,maxZBorder);
                
                XCropVol=false(1,size(vol,2));
                XCropVol(aminXBorder:amaxXBorder)=true;
                
                YCropVol=false(1,size(vol,1));
                YCropVol(aminYBorder:amaxYBorder)=true;
                
                ZCropVol=false(1,size(vol,3));
                ZCropVol(ceil((aminZBorder:amaxZBorder)*MD.pixelSize_/MD.pixelSizeZ_))=true;

                vol(:,:,~ZCropVol)=[];
                vol(~YCropVol,:,:)=[];
                vol(:,~XCropVol,:)=[];
                
            end
            
            maskedVol=vol;
            maskedVol(~mask)=0;
        end
        
        % If needed the map must rotated before cropped (efficiency)
        % Rotation depends on the FrameOfRef associated to the tracks the compose the dynanimc polygon
        % Cropping area according to the polygon OVER TIME plus added vizualiation margin
        % Rotation will use imwarp
        % Can we use imwar for cropping too ?
        
        %% if a FoF is specified, warp and crop data according to the
        tform=affine3d();
        warpedVol=vol;
        warpedMaskedVol=[];
        warpedMask=mask;
        if(~isempty(insetDynROI))
            warpedMaskedVol=maskedVol;
        else
            warpedMaskedVol=zeros(size(warpedVol));
        end
        if(~isempty(p.FoF))
            B=p.FoF.getBase(fIdx);
            tform.T(4,[1 2 3])=(-p.FoF.getOrigAtFrame(fIdx)+p.FoF.origin(1,:))*B;
            tform.T(1:3,1:3)=B;
            %
            tformTransOnly=affine3d();
            tformTransOnly.T(4,[1 2 3])=(-p.FoF.getOrigAtFrame(fIdx));
            
            %
            tformRelTransOnly=affine3d();
            tformRelTransOnly.T(4,[1 2 3])=(-p.FoF.origin(1,:)+p.FoF.getOrigAtFrame(fIdx));
            
            tformRotOnly=affine3d();
            B=p.FoF.getBase(fIdx);
            tformRotOnly.T(1:3,1:3)=B;
            
            tformRotOnlyInit=affine3d();
            B=p.FoF.getBase(1);
            tformRotOnlyInit.T(1:3,1:3)=B;
            
            orig=p.FoF.getOrigAtFrame(fIdx);
            
            inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
                [1 MD.getDimensions('X')],[1 MD.getDimensions('Y')],[1 MD.getDimensions('Z')*MD.pixelSizeZ_/MD.pixelSize_]);
            
            
            switch p.transType
                case 'affineOnePass'
                    maxXBorderFull=MD.getDimensions('X');
                    maxYBorderFull=MD.getDimensions('Y');
                    maxZBorderFull=MD.getDimensions('Z')*(MD.pixelSizeZ_/MD.pixelSize_);
                    minXBorderFull=1;
                    minYBorderFull=1;
                    minZBorderFull=1;
                    
                    minXBorderCurr=minXBorderFull - orig(1); maxXBorderCurr=maxXBorderFull -  orig(1);
                    minYBorderCurr=minYBorderFull - orig(2); maxYBorderCurr=maxYBorderFull  - orig(2);
                    minZBorderCurr=minZBorderFull - orig(3); maxZBorderCurr=maxZBorderFull -  orig(3);
                    
                    inputRef=imref3d([ MD.getDimensions('Y') MD.getDimensions('X') MD.getDimensions('Z')], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    minXBorderCurr=minXBorder ;
                    maxXBorderCurr=maxXBorder ;
                    minYBorderCurr=minYBorder ;
                    maxYBorderCurr=maxYBorder ;
                    minZBorderCurr=minZBorder ;
                    maxZBorderCurr=maxZBorder ;
                    
                    rotOutputRef=imref3d([    ceil(maxYBorderCurr-minYBorderCurr) ...
                        ceil(maxXBorderCurr-minXBorderCurr) ...
                        ceil(maxZBorderCurr-minZBorderCurr) ], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    
                    warpedVol=imwarp(vol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                    
                    if(~isempty(insetDynROI))
                        if(p.suppressROIBorder)
                            warpedMask=imwarp(mask,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                            warpedMaskedVol=warpedVol(warpedMask==1);
                        else
                            warpedMaskedVol=imwarp(maskedVol,inputRef,tformRotOnly,'OutputView',rotOutputRef);
                        end
                    else
                        warpedMaskedVol=zeros(size(warpedVol));
                    end
                    
                case 'translation'
                    disp(num2str(fIdx))
                    minXBorderCurr=minXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    maxXBorderCurr=maxXBorder ;%+ orig(1) - p.FoF.origin(1,1);
                    minYBorderCurr=minYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    maxYBorderCurr=maxYBorder ;%+ orig(2) - p.FoF.origin(1,2);
                    minZBorderCurr=minZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    maxZBorderCurr=maxZBorder ;%+ orig(3) - p.FoF.origin(1,3);
                    
                    %             [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformTransOnly,[minXBorder maxXBorder], [minYBorder maxYBorder], [minZBorder maxZBorder]);
                    %             minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
                    %             minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
                    %             minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);
                    %
                    outputRef=imref3d([ ceil(maxYBorderCurr-minYBorderCurr) ...
                        ceil(maxXBorderCurr-minXBorderCurr) ...
                        ceil(maxZBorderCurr-minZBorderCurr) ], ...
                        [minXBorderCurr maxXBorderCurr], [minYBorderCurr maxYBorderCurr], [minZBorderCurr maxZBorderCurr]);
                    
                    warpedVol=imwarp(vol,inputRef,tformTransOnly,'OutputView',outputRef);
                    warpedMaskedVol=imwarp(maskedVol,inputRef,tformTransOnly,'OutputView',outputRef);
                    
                case 'transCrop'
                    %% to be updated
                    [xLimitsOut,yLimitsOut,zLimitsOut] = outputLimits(tformRelTransOnly,[minXBorder maxXBorder], [minYBorder maxYBorder], [minZBorder maxZBorder]);
                    minXBorderCurr=xLimitsOut(1); maxXBorderCurr=xLimitsOut(2);
                    minYBorderCurr=yLimitsOut(1); maxYBorderCurr=yLimitsOut(2);
                    minZBorderCurr=zLimitsOut(1); maxZBorderCurr=zLimitsOut(2);
                    
                    maskcrop=maskedVol;
                    nullMaskXY=(squeeze(any(maskcrop,3)));
                    YNull=~(squeeze(any(any(mask,3),2)));
                    XNull=~(squeeze(any(any(mask,3),1)));
                    ZNull=~(squeeze(any(any(mask,1),2)));
                    
                    YNull= zeros(1,size(maskcrop,1));
                    YNull(1:minYBorderCurr)=1;
                    YNull(maxYBorderCurr:end)=1;
                    YNull=logical(YNull);
                    
                    XNull= zeros(1,size(maskcrop,2));
                    XNull(1:minXBorderCurr)=1;
                    XNull(maxXBorderCurr:end)=1;
                    XNull=logical(XNull);
                    
                    ZNull= zeros(1,size(maskcrop,3));
                    ZNull(1:ceil(minZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_))=1;
                    ZNull(ceil(maxZBorderCurr*MD.pixelSize_/MD.pixelSizeZ_):end)=1;
                    ZNull=logical(ZNull);
                    
                    maskcrop(:,:,ZNull)=[];
                    maskcrop(YNull,:,:)=[];
                    maskcrop(:,XNull,:)=[];
                    
                    warpedMaskedVol=maskcrop;
                    
                    warpedVol=vol;
                    warpedVol(:,:,ZNull)=[];
                    warpedVol(YNull,:,:)=[];
                    warpedVol(:,XNull,:)=[];
                    
                otherwise
                    error('unknown trans type');
            end
        end

        
        %% Create MIPS for each channel, fuse mask and full volume
        ZRatio=1;
        switch p.transType
            case 'none'
            case 'transCrop'
                ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
            otherwise
                ZRatio=1;
        end;
        if(isempty(p.FoF))
            ZRatio=MD.pixelSizeZ_/MD.pixelSize_;
        end
        [fullmaxXY,fullmaxZY,fullmaxZX,~]=computeMIPs(warpedVol,ZRatio, ...
            minIntensityNorm(chIdx),maxIntensityNorm(chIdx),'raw',rawTIFF);
        [maxXY,maxZY,maxZX,~]=computeMIPs(warpedMaskedVol,ZRatio, ...
                minIntensityNorm(chIdx),maxIntensityNorm(chIdx),'raw',rawTIFF);
        [maskXY,maskZY,maskZX,~]=computeMIPs(warpedMask,ZRatio, ...
                minIntensityNorm(chIdx),maxIntensityNorm(chIdx),'raw',true);

        % Fuse ROI and context
        if(p.suppressROIBorder)
            % Create MIP of ROI and context
            maxXY(imerode(maskXY,ones(3))<1)=fullmaxXY(imerode(maskXY,ones(3))<1);
            maxZY(imerode(maskZY,ones(3))<1)=fullmaxZY(imerode(maskZY,ones(3))<1);
            maxZX(imerode(maskZX,ones(3))<1)=fullmaxZX(imerode(maskZX,ones(3))<1);         
        else
            % Create MIP of ROI and context

            maxXY(maxXY==0)=fullmaxXY(maxXY==0);
            maxZY(maxZY==0)=fullmaxZY(maxZY==0);
            maxZX(maxZX==0)=fullmaxZX(maxZX==0);                        
        end

        
        %% Resize and fuse channel MIPS
        maxMIPSize=p.maxMIPSize;
        [sX,sY,sZ]=size(warpedMaskedVol);
         sZ=sZ*ZRatio;
        resizeScale=maxMIPSize/max([sX,sY,sZ]);
        
        XYMax=imresize(maxXY,resizeScale,'nearest');
        ZYMax=imresize(maxZY,resizeScale,'nearest');
        ZXMax=imresize(maxZX,resizeScale,'nearest');
        
        if(rawTIFF)
            mips{1,chIdx} = uint8((2^8-1)*mat2gray(XYMax,double([minIntensityNorm(chIdx),maxIntensityNorm(chIdx)])));
            mips{2,chIdx} = uint8((2^8-1)*mat2gray(ZYMax,double([minIntensityNorm(chIdx),maxIntensityNorm(chIdx)])));
            mips{3,chIdx} = uint8((2^8-1)*mat2gray(ZXMax,double([minIntensityNorm(chIdx),maxIntensityNorm(chIdx)])));
        else
            mips{1,chIdx} = XYMax;
            mips{2,chIdx} = ZYMax;
            mips{3,chIdx} = ZXMax;
        end
            
       
        if(~isempty(p.processSingleProj))
            projectDynROIProcess.saveFrame(chIdx,fIdx,maxXY,maxZY,maxZX);
        end
        


        %% save sparse Mask volume
        if(~isempty(p.processMaskVolume))
            sparseMask=ndSparse(double(warpedMaskedVol));
            saveMask(sprintfPath(p.processMaskVolume.outFilePaths_{1,chIdx},fIdx),sparseMask)
        end
    end
    
    if(~isempty(p.processRenderer))
        %% fuse volume if 2 channels
        if(length(p.renderedChannel)==2)
            XYProj=renderChannel(mips{1,1},mips{1,2},p.channelRender);
            ZYProj=renderChannel(mips{2,1},mips{2,2},p.channelRender);
            ZXProj=renderChannel(mips{3,1},mips{3,2},p.channelRender);
        else
            XYProj=repmat(mips{1,p.renderedChannel(1)},1,1,3);
            ZYProj=repmat(mips{2,p.renderedChannel(1)},1,1,3);
            ZXProj=repmat(mips{3,p.renderedChannel(1)},1,1,3);
        end
        
        %% write images
        if(~isempty(p.processRenderer))
            p.processRenderer.saveFrame(1,fIdx,XYProj,ZYProj,ZXProj);
        end
    end
end

if(~isempty(p.processRenderer)) 
    ProjAnimation(p.processRenderer,'ortho').saveVideo([p.processRenderer.getOutputDir()  '.avi']);
end

function RGBVol=renderChannel(ch1,ch2,type,varargin)
    switch(type)
        case 'grayRed'
          RGBVol=grayRedRender(ch1,ch2);
        case 'greenRed'
            RGBVol=greenRedRender(ch1,ch2);
        otherwise
            error('unknown channel renderer');
    end

function RGBVol=greenRedRender(greenCh,redCh)
    RGBVol=repmat(greenCh,1,1,3);
    %RGBThree(:,:,1)=max(rMaxXY,rmaxXYKin);
    RGBVol(:,:,1)=redCh;
    RGBVol(:,:,2)=greenCh;
    RGBVol(:,:,3)=0;


function RGBVol=grayRedRender(grayCh,redCh)
    RGBVol=repmat(grayCh,1,1,3);
    RGBVol(:,:,1)=max(grayCh,redCh);
    
    

function saveMask(filepath,sparseMask)
    save(filepath,'sparseMask');