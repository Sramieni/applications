classdef FlowAnalysisProcess < DataProcessingProcess
    % Concrete class for a process analyzing flow
    %
    % Sebastien Besson, 5/2011
    
    properties (SetAccess = protected)  
        flowLimits_
        speedMapLimits_ 
        protSpeedMapLimits_
    end
    
    methods
        function obj = FlowAnalysisProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = FlowAnalysisProcess.getName;
                super_args{3} = @analyzeMovieFlow;
                if isempty(funParams)
                    funParams=FlowAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@DataProcessingProcess(super_args{:});
            
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'speedMap','Md','Ms','E','S','img3C_map','img3C_SNR','protSpeedMap'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            
            % Data loading
            output = ip.Results.output;
            if ischar(output), output = {output}; end
            
            % Read file name
            outFileNames = arrayfun(@(x) x.name,...
                dir([obj.outFilePaths_{1,iChan} filesep '*.mat']),'Unif',false);
            for j=1:numel(output)
                varargout{j} = cell(size(iFrame));
            end
            
            % Load output
            for i=1:numel(iFrame)
                kineticMapFile= [obj.outFilePaths_{1,iChan}...
                    filesep outFileNames{iFrame(i)}(1:end-4) '.mat'];
                s = load(kineticMapFile,output{:});
                for j=1:numel(output)
                    varargout{j}{i} = s.(output{j});
                end
            end
            if numel(iFrame)==1,
                for j=1:numel(output)
                    varargout{j} = varargout{j}{1};
                end
            end
            
        end
        
        function setSpeedMapLimits(obj,speedMapLimits)
            obj.speedMapLimits_=speedMapLimits;
        end
        
        function setProtSpeedMapLimits(obj,speedMapLimits)
            obj.protSpeedMapLimits_=speedMapLimits;
        end
        
        function setFlowLimits(obj,flowLimits)
            obj.flowLimits_=flowLimits;
        end        
        
        function output = getDrawableOutput(obj)
            colors = hsv(numel(obj.owner_.channels_));
            output(1).name='Speed map';
            output(1).var='speedMap';
            output(1).formatData=[];
            output(1).type='image';
            output(1).defaultDisplayMethod=@(x)ImageDisplay('Colormap','jet',...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.speedMapLimits_{x});
            output(2).name='Interpolated vectors';
            output(2).var='Md';
            output(2).formatData=@(x)[x(:,[2 1]) x(:,[4 3])-x(:,[2 1])];
            output(2).type='overlay';
            output(2).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',colors(x,:),...
                'Colormap',jet,'CLim',obj.flowLimits_{x});
            output(3).name='Noise vectors';
            output(3).var='Ms';
            output(3).formatData=@(x)[x(:,[2 1]) x(:,[4 3])-x(:,[2 1])];
            output(3).type='overlay';
            output(3).defaultDisplayMethod=@(x) VectorFieldDisplay('Color',colors(x,:));
            output(4).name='Error circles';
            output(4).var='E';
            output(4).formatData=@(x) [x(:,1:3) x(:,3)];
            output(4).type='overlay';
            output(4).defaultDisplayMethod=@(x) RectangleDisplay('Color',colors(x,:),...
                'Curvature',[1 1]);
            output(5).name='SNR circles';
            output(5).var='S';
            output(5).formatData=@(x) [x(:,1:3) x(:,3)];
            output(5).type='overlay';
            output(5).defaultDisplayMethod=@(x) RectangleDisplay('Color',colors(x,:),...
                'Curvature',[1 1]);
            output(6).name='Error map';
            output(6).var='img3C_map';
            output(6).formatData=[];
            output(6).type='image';
            output(6).defaultDisplayMethod=@ImageDisplay;
            output(7).name='SNR map';
            output(7).var='img3C_SNR';
            output(7).formatData=[];
            output(7).type='image';
            output(7).defaultDisplayMethod=@ImageDisplay;
            
            output(8).name='Protrusion Speed map';
            output(8).var='protSpeedMap';
            output(8).formatData=[];
            output(8).type='image';
            output(8).defaultDisplayMethod=@(x)ImageDisplay('Colormap',BlueBlackRedColorMap,...
                'Colorbar','on','Units',obj.getUnits,'CLim',obj.protSpeedMapLimits_{x}); %speedMapLimits_{x}); %
        end  
        
    end
    
    methods (Static)
        function name =getName()
            name = 'Flow Analysis';
        end
        function h = GUI()
            h= @flowAnalysisProcessGUI;
        end
        
        function procNames =getFlowProcesses()
            procNames = {'SpeckleTrackingProcess';...
                'FlowTrackingProcess'};
        end

        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'flowAnalysis'];
            funParams.FlowProcess = 'SpeckleTrackingProcess';
            funParams.timeWindow = 3;
            funParams.corrLength = 33;
            funParams.gridSize = 11;
            funParams.noise = 1;
            funParams.error = 1;
            funParams.driftCorrection = false;
        end 
        function units = getUnits(varargin)
            units = 'Speed (nm/min)';
        end
    end
end