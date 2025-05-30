function [deviationLevels,indWrongVectors, indMissing, MSE,DTM,Accu] =  compareFields(field1,field2,plotResult,sc)
%function [deviationLevels,indMissing,MSE] =  compareFields(field1,field2)
%compares field1 against field2 and reports how much each vector deviates from
%field2 (deviationLevels), indexes of wrong and missing vectors, and mean
%squared errors.
%   input:
%       field1:         a field to be tested. it has Pos and Vec
%                       fields. e.g.: displField(1)
%       field2:         a field to be used as ground truth. it has x_mat_u,
%                       y_mat_u, ux and uy,  as fields.
%       plotResult:     true if you want to see the results (default:
%                       false)
%       sc:             scale for the quiver 1 for default
%   output:
%       deviationLevels     a Nx1 vector where each component corresponds
%                           to each position in field1, compared to
%                           matching vector in field2. It is calculated as
%                           vector difference divided by original vector
%                           magnitude.
%       indWrongVectors     a Nx1 logical array where true indicates
%                           vectors in field1 that have more than 10% of
%                           deviation.
%       indMissing          a Nx1 logical array where true indicates
%                           vectors in field1 that have NaNs
%       MSE                 Mean-squared-error, i.e., sum of vector
%                           deviations divided by sum of vector norms of
%                           original vectors (in field2) at corresponding
%                           positions
%       DTM                 Deviation of tractio magnitude (from Sabass et
%                           al 2008 BiophysJ)
%       Accu                Accuracy of field1 compared to field2
%                           quantified as sum(1-deviationLevel)/N.
% Sangyoon Han, March 2021

if nargin<3
    plotResult =  false;
end
%% get the field1-corresponding vectors from field2
nPoints = size(field1.pos,1);
bead_x  = field1.pos(:,1);
bead_y  = field1.pos(:,2);
ux2 = zeros(size(bead_x));
uy2 = zeros(size(bead_y));
x_mat_u = field2.x_mat_u;
y_mat_u = field2.y_mat_u;
ux = field2.ux;
uy = field2.uy;

for k=1:nPoints
    [~,indcol_closest_x] = min(abs(x_mat_u(1,:)-bead_x(k)),[],2);
    [~,indrow_closest_y] = min(abs(y_mat_u(:,1)-bead_y(k)),[],1);
    row_bottom = max(1,indrow_closest_y-2);
    row_top = min(size(x_mat_u,2),indrow_closest_y+2);
    col_bottom = max(1,indcol_closest_x-2);
    col_top = min(size(y_mat_u,1),indcol_closest_x+2);
    loc_xmat = x_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ymat = y_mat_u(row_bottom:row_top,col_bottom:col_top);
    loc_ux = ux(row_bottom:row_top,col_bottom:col_top);
    loc_uy = uy(row_bottom:row_top,col_bottom:col_top);
    ux2(k) = interp2(loc_xmat,loc_ymat,loc_ux,bead_x(k),bead_y(k));
    if isnan(ux2(k))
        ux2(k) = ux(indrow_closest_y,indcol_closest_x);
    end
    uy2(k) = interp2(loc_xmat,loc_ymat,loc_uy,bead_x(k),bead_y(k));
    if isnan(uy2(k))
        uy2(k) = uy(indrow_closest_y,indcol_closest_x);
    end
end
% hold on;
% quiver(bead_x,bead_y,ux2,uy2,0,'Color','m')
%% deviationLevels
ux1  = field1.vec(:,1);
uy1  = field1.vec(:,2);
% We replace NaNs in ux1 and uy1 with zeros to include them for deviation
% calculation.
ux1d = ux1; ux1d(isnan(ux1d))=0;
uy1d = uy1; uy1d(isnan(uy1d))=0;

deviation =  ((ux1d-ux2).^2 + (uy1d-uy2).^2).^0.5;
magOrg  = mean((ux2.^2 + uy2.^2).^0.5)+0.5*std((ux2.^2 + uy2.^2).^0.5); 
deviationLevels = deviation / magOrg; %
% % I will use the average magnitude to normalize the deviation
% %To prevent over quantification for small vectors, I'll use magOrg for
% %small vectors but use their own mag for large enough vectors
% magD2  = (ux2.^2 + uy2.^2).^0.5;
% noiseFloorFactor = 0.8;
% floorMag = noiseFloorFactor*magOrg;
% magInteg = floorMag.* (magD2 < floorMag) + magD2.* (magD2 >= floorMag);
% deviationLevelsForAcc = deviation ./ magInteg; %/magOrg; %

if nargin <4 
    sc=1;
end
if  plotResult
    figure
    quiverColormap(bead_x,bead_y,sc*ux1,sc*uy1,'customIntensity',deviationLevels);
%     cBar=colorbar; colormap jet
%     cBar.Limits=[min(deviationLevels) max(deviationLevels)];
%     quiverColormap(bead_x,bead_y,ux1-ux2,uy1-uy2);
% quiverColormap(bead_x,bead_y,ux2,uy2);
end
%% indWrongVectors
thresDev = 0.5;
indWrongVectors = deviationLevels>thresDev;
%% indMissing
indMissing = isnan(ux1);
%% MSE
% MSE = nansum(deviation)/nansum(mag);
MSE = nansum(deviationLevels)/numel(ux1d); % the higher MSE is, the more deviation is.
% MSE = nansum(deviation)/numel(ux1d); % the higher MSE is, the more deviation is. 
% Now the unit is the same as the ux1.
%% DTM
deviationMag =  abs((ux2.^2 + uy2.^2).^0.5 - (ux1d.^2 + uy1d.^2).^0.5)./(ux2.^2 + uy2.^2).^0.5;
numValVecs = sum(~isnan(deviationMag));
DTM  = sum(deviationMag)/numValVecs;
%% Accuracy
% thresGoodDev=0.5;
Accu = nansum(1-deviationLevels)/numel(ux1d); 
% Accu = nansum(deviationLevels<thresGoodDev)/numel(ux1d); 
% debug
% figure, quiver(field1.pos(:,1),field1.pos(:,2),field1.vec(:,1),field1.vec(:,2),0)
% hold on
% quiver(field1.pos(:,1),field1.pos(:,2),ux2,uy2,0)
end

