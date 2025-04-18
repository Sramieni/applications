%Script to extract compTracks from receptorInfoLabeled located in the
%indicated source folders. Note that receptorInfoLabeled has multiple entries
%for different label ratios used, which will each get its own directory.
%
%Khuloud Jaqaman, May 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161109/';

%Define strings for directory hierarchy as needed
 rDDir = {'rD4'};%,'rD60','rD80','rD120','rD140','rD160'};
 aPDir = {'aP0p0'};
outDirNum =1:5;
lRDir = {'lR0p2','lR0p4'};

%Define number of label ratio
numLabelRatio = length(lRDir);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
      
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        compTracksVec = cell(numLabelRatio,1);
        
        %Original output is not organized by label ratio since
        %receptorInfoLabeled for each label ratio is saved as a struct.
        %Iterate through the outputs, to pull out each receptorInfoLabeled
        %then the compTracks.  
        for outDirIndx = 1 : length(outDirNum)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx))];
            
            %Load receptorInfoLabeled

            tempRecepInfo = load([currDir,filesep,...
                'receptorInfoLabeled',int2str(outDirNum(outDirIndx)),'.mat']);
            
            %Pull out compTracks for each labelRatio defined above
            [compTracksVec{:}] = tempRecepInfo.receptorInfoLabeled(1:numLabelRatio).compTracks;

            %For each label ratio, the inner most directory, create the
            %directory and save compTracks.            
            for lRDirIndx=1:numLabelRatio
                
                fprintf('\n   Out = %d, lR = %s ',outDirIndx,lRDir{lRDirIndx});
                
                currOutDir = [currDir,filesep,lRDir{lRDirIndx}];
                
                %Create the direcotry
                mkdir(currOutDir)
                
                %Write compTracks
                compTracks = compTracksVec{lRDirIndx};
                save([currOutDir,'/compTracks'],'compTracks','-v7.3');
                
                clear compTracks
                
                fprintf('... done.');
                
            end %for each labelRatio
            
            clear compTracks tempRecepInfo
            
        end %for each outDir
        
        clear compTracksVec
                            
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear



    
