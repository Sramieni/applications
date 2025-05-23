%%load data

clear; clc;
stiffness=["0.6","1.3","2.6","6","12.7"];
numStiff=length(stiffness);
structCell=cell(numStiff,1);
cmap=jet(5);
%for d=[1:numDrugs]
    for s=[1:numStiff]
        if usejava('desktop')
            %[fileSFolders, pathSFolders] = uigetfile('*.mat',['Select ' char(drugs(d)) ' ' char(stiffness(s))]);
       
            [fileSFolders, pathSFolders]=uigetfile('title',['Select ' char(stiffness(s))]);
        try 
            dataFile=load(append(pathSFolders,fileSFolders));
            data=dataFile.data;
            structCell{s}=data;
            %structCell{s,d}=load(append(stiffnessFolder,'/',f(d).name));
        catch
            disp(['Error :: failed to load file '  fileSFolders])
        end
        end
    end

%% Figures
    
len=zeros(1,numStiff);
for i=1:numStiff
    mm=[];
    for j=1:length(structCell{i})
        directory=path2dir(structCell{i}(j).file);
        speed=load([directory '/WindowingPackage/window_sampling/Speed map - channel 1.mat']);
        speedOut = selectWindows(speed,structCell{i}(j).selected);
        if ~isempty(speedOut)
            len(i)=len(i)+size(structCell{1}(1).selected,2);
            [pp(j,:),ff(j,:),dd(j,:),m]=windowFFT(speedOut);
            mm=[mm,m];
        end
    end
    d(i,:)=sum(dd,1);
    fp(i,:)=mean(pp,'omitnan');
    mom{i}=mm;
end

figure();
hold on;
for p=1:numStiff
    tt=plot(ff(1,:),normalize(fp(p,:)),'Color',cmap(p,:),'LineWidth',2,'Marker','o');
    uistack(tt,'top')
end
ylabel("Response")
xlabel("Frequency (Hz)")
lab={['0.6'  'kPa'],['1.3'  'kPa'],['2.6' 'kPa'],['6'  'kPa'],['12.7' 'kPa']};
legend(lab,'Location','eastoutside')

d=round(d);
maxLength=max(sum(d,2));
density=nan(5,maxLength);
for i = 1:5
    k=1;
    for j = 1:length(d(i,:))
        density(i,k:k+d(i,j)-1)=ff(i,j);
        k=k+d(i,j);
    end
end

n=size(density,1);
densityCell=cell(1,n);
names={'0.6','1.3','2.6','6','12.7'};


for i=1:n
    nDen=density(i,:);
    nDen=imresize(nDen(~isnan(nDen)),[1,len(i)*16]);
    densityCell{i}=nDen;
end
figure();
boxPlotCellArray(mom,names,1,1,1);
%boxPlotCellArray(densityCell,names,1,1,1);
xlabel('Stiffness (kPa)');
ylabel('Frequency');

figure();
errorBarPlotCellArray(mom,names,1);

xlabel('Stiffness (kPa)');
ylabel('Frequency');


%% functions
function [directory] = path2dir(path)
    directory=strsplit(path,filesep);
    directory=directory(1:end-1);
    directory=strjoin(directory,filesep);
    directory=char(directory);

end

function [speedOut] = selectWindows(speed,windows)
        
    speedOut = zeros(length(windows),31);
    for i = 1:length(windows)
        speedOut(i,:)=squeeze(speed.samples.avg(windows{i}(1),windows{i}(2),:));
    end

end

function [edgeOut,speedOut,iOut] = threshCell(edge,speed,t)

    protrudingWindows = any(edge>t,2);
    timeProt=sum(edge>t,2)>7;
    %protrudingWindows=protrudingWindows(timeProt);
    sum(protrudingWindows);
    %maxOutCols = quantile(edge,0.9,2);
    %meanOutCols = mean(edge,2);
    maxOutCols=quantileIndex(edge,.7);
%     numberWindows=sum(maxOutCols(protrudingWindows))
    edgeOut=edge(maxOutCols(protrudingWindows),:);
    
    speedOut=speed(maxOutCols(protrudingWindows),:);
    iOut=maxOutCols(protrudingWindows);
end

function indexes = quantileIndex(arr,p)
    perc=prctile(arr,p);
    indexes=any(arr<=perc,2);
end

function [] = mapFlowAndProtrusion(edge,speed)
    figure, imagesc(edge),
    axis xy
    GreenBlackRedColorMap; caxis([-900 900]),colorbar
    title('Protrusion')
    % flow speed
  
    figure, imagesc(speed),
    axis xy
    colormap jet; caxis([0 2000]),colorbar
    title('Actin speed')
end

function [mm] = hhtSpeed(speed,fig,c)
    [u,n]=size(speed);
    if u==0
        mm=[0,0,0,0,0,0];
        return
    end
    cmap=jet(5);
    figure(fig);
    hold on
    time_frame=31;
    fs=1/6;
    specRes=20;
    for i=1:u
        speed(find(isnan(speed)))=0;
        [instAmp_time, instFreq_time,hht]=hilbertHuangTransform(speed(i,:), fs, specRes);
        if size(instFreq_time,1) < 6
           for l = size(instFreq_time,1):6
               instAmp_time(l,:)=zeros(1,time_frame-1); instFreq_time(l,:)=zeros(1,time_frame-1);
           end
        end
        for p=1:6
            
            m(i,p)=mean(instFreq_time(p,:),'omitnan'); 
            st(i,p)=std(instFreq_time(p,:),'omitnan');
        end
    end
    [l,~]=size(m);
    for x=1:l
                errorbar(m(x,:),st(x,:),'Color',cmap(c,:),'LineWidth',0.1);
                hold on;
    end
    
    xlabel('IMF number','FontSize',20);
    ylabel('Frequency (Hz)','FontSize',20);
    axis([0.5 6.5 0-0.005 fs/2+0.005]);
    set(gca, 'XTick', [1:6], 'XTickLabel', [1:6]);
    set(gca, 'FontSize', 20); 
    mm=mean(m,'omitnan');
end


function [p1,f,d,mp] = windowFFT(speed)
            fs=1/6;
           
            [n,m]=size(speed);
            if n==0
                d=zeros(1,16);
                p1=NaN(1,16);
                f=fs/31*(0:floor((31-1)/2));
                return
            end
            for i=1:n
                l=m;
                ff=fft(speed(i,:)-mean(speed(i,:),'omitnan'));
                p=abs(ff);
                %p=normalize(p);
                pp(i,:)=p(1:floor((l+1)/2));
                f=fs/l*(0:floor((l-1)/2));
                
            end
            mp=sum((f.*pp)./sum(pp,2),2);
            mp=mp(~isnan(mp))';
            d=sum(pp,1,'omitnan');
            p1=mean(pp,'omitnan');
            %d=downsample(d,100);
end
