
[XYZ] = matrix2coords(matrixx);

XYZ(XYZ(:,3)== 0, :) = [];

[maxR, maxID] = max(XYZ(:,3));
centerX = XYZ(maxID, 1);
centerY = XYZ(maxID, 2);

% Find distance from center to all points
for k = 1 : length(XYZ)
  distances(k) = sqrt((centerX- XYZ(k,1))^2 + (centerY- XYZ(k,2))^2);
end
 [maxDistance, indexOfMax] = max(distances);

 maxX = XYZ(indexOfMax,1);
 maxY = XYZ(indexOfMax,2);
 % Find rotation angle from center to farthest away point
theta = atan((centerY-centerY)/(0-centerX)) - atan((centerY-maxY)/(centerX-maxX));

%rotate the matrix around the center
%move all points so center is origin
tXY(:,1) = XYZ(:,1)-centerX;
tXY(:,2) = XYZ(:,2)-centerY;
% rotate by theta around (0,0)
xRot     = tXY(:,1)*cos(theta) - tXY(:,2)*sin(theta);
yRot     = tXY(:,1)*sin(theta) + tXY(:,2)*cos(theta);


 % estimating center and dispersion along x axis
xd = [xRot, XYZ(:,3)];
xdsort = sortrows(xd);
x_collapse = -100;
xz_collapse = [0,0];
for i = 1:length(xdsort)     
    if round(xdsort(i,1)*1000)/1000 == round(x_collapse(length(x_collapse),1)*1000)/1000;
        xz_collapse(length(xz_collapse),2)= xz_collapse(length(xz_collapse),2)+ xdsort(i,2);
    else
        xz_collapse((length(xz_collapse)+1),:) = xdsort(i,:);
        x_collapse((length(x_collapse)+1),1) = xdsort(i,1);
    end
end
    x_collapsed = x_collapse(2:end);
    xz_collapsed = xz_collapse(3:end,:);
    
centriod_x(cell_n) = wmean(xz_collapsed(:,1),xz_collapsed(:,2));
dispersion_x(cell_n) = std(xz_collapsed(:,1),xz_collapsed(:,2));
 
 % estimating center and dispersion along y axis
yd = [yRot, XYZ(:,3)];
ydsort = sortrows(yd);

y_collapse = -100;
yz_collapse = [0,0];
for i = 1:length(ydsort)     
    if round(ydsort(i,1)*1000)/1000 == round(y_collapse(length(y_collapse),1)*1000)/1000;
        yz_collapse(length(yz_collapse),2)= yz_collapse(length(yz_collapse),2)+ ydsort(i,2);
    else
        yz_collapse((length(yz_collapse)+1),:) = ydsort(i,:);
        y_collapse((length(y_collapse)+1),1) = ydsort(i,1);
    end
end
    y_collapsed = y_collapse(2:end);
    yz_collapsed = yz_collapse(3:end,:);
    
centriod_y(cell_n) = wmean(yz_collapsed(:,1),yz_collapsed(:,2));
dispersion_y(cell_n) = std(yz_collapsed(:,1),yz_collapsed(:,2));

obs = matrixx(:); %reduce to one dimension
    [xgrid, ygrid]= meshgrid(1:size(matrixx,2),1:size(matrixx,1));
    x(:,1) = xgrid(:);
    x(:,2) = ygrid(:);
    coeff = nlinfit(x,obs,@rf2fit,[mean_spont(cell_n) max(obs)-mean_spont(cell_n) ...
        centroid_x(cell_n) dispersion_x(cell_n) centroid_y(cell_n) dispersion_y(cell_n)]);
    est = rf2fit(coeff,x);
    figure; subplot(2,1,2); imagesc(reshape(est,[size(matrixx,1) size(matrixx,2)]));
    axis equal;  axis([1 n_col 1 n_rows]); colorbar; title('rf fit');
    rf_amp = coeff(2);
    spont = coeff(1);
    x0 = coeff(3)*spacing;
    wx0 = abs(coeff(4))*spacing;
    y0 = coeff(5)*spacing;
    wy0 = abs(coeff(6))*spacing;
    note_text = sprintf('RF fit: Spont %0.2g, Amp %0.2g, X0 %0.2gdeg, stdX %0.2gdeg, Y0 %0.2gdeg, stdY %0.2gdeg', ...
        spont, rf_amp, x0, wx0, y0, wy0);
%     figure(RF_figure);
    axes('Position',[0 0 1 1],'Visible','off'); 
    text(.025,0.025,note_text,'FontSize',10);

 
