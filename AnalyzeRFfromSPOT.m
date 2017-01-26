% script finds ON and OFF receptive fields using spot stimulus

clear all
cells=0;
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block};

[afname, apname] = uigetfile('*.mat','analysis data');
oldapname = apname;
load(fullfile(apname,afname));   %%% need to copy pname, or it can get written over in load
apname = oldapname;
AnalysisFileName = ['RF_spot_analysis_', Block_Name];
stimulus = TDT2mat(['C:\data\TDT\',Tank_Name], Block_Name, 'TYPE', 2); %may need to change directory if data moved

duration = 0.5; %500 ms
waitInt = 0.05; %not used
blankframe=1; %1 = yes
spacing = 5; %5 degrees


%xTrg can only hold up to 255, made extra fTrg of two to indicate stimcond+
%255
StimOnset = stimulus.epocs.xTrg.onset;
fTrg = stimulus.epocs.fTrg.data;
high_cond = 0;
if max(fTrg) > 1
    gridX = 27; % number of X positions on grid
    gridY = 15; % number of Y positions on grid
    n_col = 27;
    n_rows = 15;
    for i = 1:length(fTrg)
        if fTrg(i)== 1
            high_cond = [high_cond; 1];
        else
            high_cond(length(high_cond)) = 2;
        end
    end   
    high_cond = high_cond(2:length(high_cond));
    for i = 1:length(high_cond)
        StimCond(i)= stimulus.epocs.xTrg.data(i) + ((high_cond(i) - 1)*255);
    end
else
    StimCond = stimulus.epocs.xTrg.data;
    gridX = 17; % number of X positions on grid
    gridY = 13; % number of Y positions on grid
    n_col = 17;
    n_rows = 13;
end

bluemap = [0, 0, 0
    0, 0, 0.1
    0, 0, 0.3
    0, 0, 0.5
    0, 0, 0.8
    0, 0, 1.0];

redmap = [0, 0, 0
    0.1, 0, 0
    0.3, 0, 0
    0.5, 0, 0
    0.8, 0, 0
    1.0, 0, 0];


nCond = gridX*gridY;
if blankframe
    nCond = nCond+1;
end
matrixx = zeros(n_rows,n_col);
ONmatrixx = zeros(n_rows,n_col);
OFFmatrixx = zeros(n_rows,n_col);
n_cells = length(spikeT);
excel_variables = cell(12, n_cells);
for m = 1: n_cells
    rf_cell(m).condition(1).trial{1}={};
    rf_cell(m).condition(:).all_trials = [];
    thresholded(m).condition(:)= [];
    block_spikes = find(spikeT{1,m}>= (block*(10^5) - 10^5) & spikeT{1,m} < (block*(10^5)));
    block_spikes2 = spikeT{m}(block_spikes);
    spikeTimes{m} = block_spikes2-((block-1)*10^5);
end
%%
%fill rf.cells with spiketimes by condition, ON OFF

for i = 1:n_cells
    cellspike = spikeTimes{i};
    trialpercond = zeros(nCond,1);
    for j = 1: length(StimOnset)
        startT = StimOnset(j);
        stopT = StimOnset(j)+(2*duration);
        trial_no = trialpercond(StimCond(j))+1;
        epochspikes = cellspike(cellspike > startT & cellspike < stopT);
        zeroedTs = epochspikes-startT; %zeroed spike times
        on_spikes = cellspike(cellspike > (startT + 0.05) & cellspike < (startT + 0.25)); %200ms after first 50ms
        off_spikes = cellspike(cellspike > (startT + 0.55) & cellspike < (startT + 0.75));
        
        rf_cell(i).condition(StimCond(j)).trial{trial_no}={zeroedTs};
        rf_cell(i).condition(StimCond(j)).all_trials = [rf_cell(i).condition(StimCond(j)).all_trials, (zeroedTs)];
        rf_cell(i).condition(StimCond(j)).ONcount{trial_no} = length(on_spikes)*5;%times 5 so it's spikes per second
        rf_cell(i).condition(StimCond(j)).OFFcount{trial_no} = length(off_spikes)*5;
        rf_cell(i).condition(StimCond(j)).EPOCHcount{trial_no} = length(epochspikes); %for whole trial: 1 second
        trialpercond(StimCond(j))= trial_no;
    end
    
    for k = 1: nCond %number of conditions
        rf_cell(i).condition(k).ONmean = mean(cell2mat(rf_cell(1,i).condition(1,k).ONcount));
        rf_cell(i).condition(k).OFFmean = mean(cell2mat(rf_cell(1,i).condition(1,k).OFFcount));
        ONmeans(i,k) = mean(cell2mat(rf_cell(1,i).condition(1,k).ONcount)); 
        OFFmeans(i,k) = mean(cell2mat(rf_cell(1,i).condition(1,k).OFFcount));
        ON_OFFmeans(i,k) = (ONmeans(i,k)+ OFFmeans(i,k))/2;
    end
    ON_OFFmeans(isnan(ON_OFFmeans))= 0;
   %correlation coefficient
   meanRon = mean(ONmeans(i,1:(n_rows*n_col)));
   meanRoff = mean(OFFmeans(i,1:(n_rows*n_col)));
   for k = 1:(n_rows*n_col)
       numerators(k)= (ONmeans(i,k)- meanRon)*(OFFmeans(i,k)-meanRoff);
       onDenomenator(k) = (ONmeans(i,k)-meanRon)^2;
       offDenomenator(k) = (OFFmeans(i,k)-meanRoff)^2;
   end
   correlation_coef(i) = (sum(numerators))/sqrt((sum(onDenomenator))*sum(offDenomenator));
   
%%       
%calculate clustering ratio
if max(ON_OFFmeans) > 1; 
matri = reshape((ON_OFFmeans(i,1:(n_rows*n_col))), n_rows, n_col);
[matrixx2] = matrix2coords(matri);
matrixx2(matrixx2(:,3)== 0, :) = [];

    XYZs = sortrows(matrixx2, -3); %sorted by highest response
    Zs = XYZs(:,3);
    x0 = XYZs(1,1);
    y0 = XYZs(1,2);

if length(Zs)> 3
    
Idx = find(Zs >= Zs(4));
for id = 1:length(Idx)
        dist2ctr(id,1) = XYZs(Idx(id),1);
        dist2ctr(id,2) = XYZs(Idx(id),2);
        dist2ctr(id,3) = sqrt((x0 - XYZs(Idx(id),1))^2 + (y0 - XYZs(Idx(id),2))^2);
end
dist2ctr = sortrows(dist2ctr, 3);
combs = [1,2;1,3;1,4;2,3;2,4;3,4];
for c = 1:6
    rf_clust(c) = sqrt((dist2ctr(combs(c,1),1)- dist2ctr(combs(c,2),1))^2+ (dist2ctr(combs(c,1),2)- dist2ctr(combs(c,2),2))^2);
end
rf_clust_ratio(i) = mean(rf_clust) *5; %in degrees
else
    rf_clust_ratio(i) = 0;
end
else
    rf_clust_ratio(i) = 0;
end

   
   
 %%  
   %thresholding
  spont_firing = cell2mat(rf_cell(1,i).condition(1,((n_rows*n_col)+1)).EPOCHcount); %1 second
  mean_spont(i) = mean(spont_firing); 
  std_spont(i) = std(spont_firing);
  std_mean(i) = std(ON_OFFmeans(i,:));
  min_fire(i) = min(ON_OFFmeans(i,:));
%   threshold(i) = mean(ON_OFFmeans(i,:));    
 threshold(i) = mean(ON_OFFmeans(i,:)) + std_mean(i); %threshold is avg spont + 2*SD
% threshold(i) = 0;
   thresh = threshold(i);
   
   for j = 1: nCond
       for k = 1: length(rf_cell(1,i).condition(1,j).trial)
           if (rf_cell(1,i).condition(1,j).ONcount{k}- thresh) > 0
              rf_cell(1,i).condition(1,j).ONpass(k) = 1;
           else
              rf_cell(1,i).condition(1,j).ONpass(k) = 0; 
           end
            if (rf_cell(1,i).condition(1,j).OFFcount{k}- thresh) > 0
              rf_cell(1,i).condition(1,j).OFFpass(k) = 1;
           else
              rf_cell(1,i).condition(1,j).OFFpass(k) = 0; 
           end          
             if (rf_cell(1,i).condition(1,j).ONcount{k} + rf_cell(1,i).condition(1,j).OFFcount{k} - (2*thresh)) > 0
              rf_cell(1,i).condition(1,j).OnOffpass(k) = 1;
           else
              rf_cell(1,i).condition(1,j).OnOffpass(k) = 0; 
           end         
       end
       rf_cell(1,i).ONpassthresh(1,j) = (mean(rf_cell(1,i).condition(1,j).ONpass) >= 0.4); %must be above thresh 40% of trials
       rf_cell(1,i).OFFpassthresh(1,j) = (mean(rf_cell(1,i).condition(1,j).OFFpass) >= 0.4);
       rf_cell(1,i).ON_OFFpassthresh(1,j) = (mean(rf_cell(1,i).condition(1,j).OnOffpass) >= 0.4);
       rf_cell(1,i).BOTHpassthresh(1,j) =  ((rf_cell(1,i).ONpassthresh(1,j) + rf_cell(1,i).OFFpassthresh(1,j)) == 2); %ON AND OFF
       rf_cell(1,i).EITHERpassthresh(1,j) = ((rf_cell(1,i).ONpassthresh(1,j) + rf_cell(1,i).OFFpassthresh(1,j)) > 0); % ON or OFF
       rf_cell(1,i).ONmeansThresh(1,j)= rf_cell(1,i).ONpassthresh(1,j) * ONmeans(i,j);
       rf_cell(1,i).OFFmeansThresh(1,j)= rf_cell(1,i).OFFpassthresh(1,j) * OFFmeans(i,j);
       rf_cell(1,i).ONOFFmeansThresh(1,j) = rf_cell(1,i).ON_OFFpassthresh(1,j) * ON_OFFmeans(i,j);
   end
   
   OnOffOverlapRatio1(i) = sum(rf_cell(1,i).BOTHpassthresh)/sum(rf_cell(1,i).EITHERpassthresh);
end
%%
%plot histograms
binranges = linspace(0,1,11); %for 100ms bins= 11, change 11 to 21 for 50ms bins etc
% % plotinx = [1,18,35,52,69,86,103,120,137,154,171,188,205,2,19,36,53,70,87,104,121,138,155,172,189,206,3,20,37,54,71,88,105,122,139,156,173,190,207,4,21,38,55,72,89,106,123,140,157,174,191,208,5,22,39,56,73,90,107,124,141,158,175,192,209,6,23,40,57,74,91,108,125,142,159,176,193,210,7,24,41,58,75,92,109,126,143,160,177,194,211,8,25,42,59,76,93,110,127,144,161,178,195,212,9,26,43,60,77,94,111,128,145,162,179,196,213,10,27,44,61,78,95,112,129,146,163,180,197,214,11,28,45,62,79,96,113,130,147,164,181,198,215,12,29,46,63,80,97,114,131,148,165,182,199,216,13,30,47,64,81,98,115,132,149,166,183,200,217,14,31,48,65,82,99,116,133,150,167,184,201,218,15,32,49,66,83,100,117,134,151,168,185,202,219,16,33,50,67,84,101,118,135,152,169,186,203,220,17,34,51,68,85,102,119,136,153,170,187,204,221];
% plotinx = [1,28,55,82,109,136,163,190,217,244,271,298,325,352,379,2,29,56,83,110,137,164,191,218,245,272,299,326,353,380,3,30,57,84,111,138,165,192,219,246,273,300,327,354,381,4,31,58,85,112,139,166,193,220,247,274,301,328,355,382,5,32,59,86,113,140,167,194,221,248,275,302,329,356,383,6,33,60,87,114,141,168,195,222,249,276,303,330,357,384,7,34,61,88,115,142,169,196,223,250,277,304,331,358,385,8,35,62,89,116,143,170,197,224,251,278,305,332,359,386,9,36,63,90,117,144,171,198,225,252,279,306,333,360,387,10,37,64,91,118,145,172,199,226,253,280,307,334,361,388,11,38,65,92,119,146,173,200,227,254,281,308,335,362,389,12,39,66,93,120,147,174,201,228,255,282,309,336,363,390,13,40,67,94,121,148,175,202,229,256,283,310,337,364,391,14,41,68,95,122,149,176,203,230,257,284,311,338,365,392,15,42,69,96,123,150,177,204,231,258,285,312,339,366,393,16,43,70,97,124,151,178,205,232,259,286,313,340,367,394,17,44,71,98,125,152,179,206,233,260,287,314,341,368,395,18,45,72,99,126,153,180,207,234,261,288,315,342,369,396,19,46,73,100,127,154,181,208,235,262,289,316,343,370,397,20,47,74,101,128,155,182,209,236,263,290,317,344,371,398,21,48,75,102,129,156,183,210,237,264,291,318,345,372,399,22,49,76,103,130,157,184,211,238,265,292,319,346,373,400,23,50,77,104,131,158,185,212,239,266,293,320,347,374,401,24,51,78,105,132,159,186,213,240,267,294,321,348,375,402,25,52,79,106,133,160,187,214,241,268,295,322,349,376,403,26,53,80,107,134,161,188,215,242,269,296,323,350,377,404,27,54,81,108,135,162,189,216,243,270,297,324,351,378,405];
% %eqt = zeros(n_cells,6);
% 
for cell_n = 1:n_cells
    channel_no = cells(cell_n,1);
    clust_no = cells(cell_n,2);
%     figure
%     bincounts= zeros((n_rows*n_col+1),11);%10 bins
%     for j = 1: nCond
%         bincounts(j,:) = histc(rf_cell(1,cell_n).condition(1,j).all_trials,binranges);
%     end
%     ymax = max(max(bincounts));
%     for j = 1:(nCond-1)%dont plot blank cond
%         subplot(n_rows,n_col,plotinx(j))
%         bar(binranges,bincounts(j,:),'histc')
%         axis([0 1 0 ymax])
%         set(gca,'XTicklabel',[])
%         set(gca,'YTicklabel',[])
%      end
% % Save the file as PNG
%     print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '1histogram'],'-dpng','-r300');
%%
%plot heatmap
[xgrid, ygrid]= meshgrid(1:n_col,1:n_rows);
xcontour = xgrid(:)';
ycontour = fliplr(ygrid(:));
figure %contour plot
matrixx = zeros(n_rows,n_col);
    for p = 1: (n_rows*n_col)
        matrixx(ycontour(p),xcontour(p))= rf_cell(1,cell_n).ONOFFmeansThresh(1,p);
        ONmatrixx(ycontour(p),xcontour(p))= rf_cell(1,cell_n).ONmeansThresh(1,p);
        OFFmatrixx(ycontour(p),xcontour(p))= rf_cell(1,cell_n).OFFmeansThresh(1,p);
    end
    contourf(matrixx);
    set(gca,'Ydir', 'reverse')
    colorbar;
    rf_cell(1,cell_n).matrixx = matrixx;
     axis equal;
%%
     % Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '2contourplot'],'-dpng','-r300');
    
[neighb] = neighbormap( matrixx );
og_matrixx = matrixx;
og_ONmatrixx = ONmatrixx;
og_OFFmatrixx = OFFmatrixx;
neighb_matrixx = zeros(n_rows, n_col);
neighb_ONmatrixx = zeros(n_rows, n_col);
neighb_OFFmatrixx = zeros(n_rows, n_col);
for row = 1:n_rows
    for col = 1:n_col
        neighbs = length(neighb(row, col).map);
        clear neighb_avgs
        neighb_avgs = zeros(neighbs,3);
        for nb = 1: neighbs
            neighb_avgs(nb,1) = og_matrixx((neighb(row,col).map(nb,1)), (neighb(row,col).map(nb,2)));
            neighb_avgs(nb,2) = og_ONmatrixx((neighb(row,col).map(nb,1)), (neighb(row,col).map(nb,2)));
            neighb_avgs(nb,3) = og_OFFmatrixx((neighb(row,col).map(nb,1)), (neighb(row,col).map(nb,2)));
        end
        nb_avg = mean(neighb_avgs);
        neighb_matrixx(row, col) = nb_avg(1);
        neighb_ONmatrixx(row,col) = nb_avg(2);
        neighb_OFFmatrixx(row,col) = nb_avg(3);
    end
end

% neighb_thresh= mean(neighb_matrixx(:))+ std(neighb_matrixx(:)*0);
% neighb_threshON= mean(neighb_ONmatrixx(:))+ std(neighb_ONmatrixx(:)*0);
% neighb_threshOFF = mean(neighb_OFFmatrixx(:))+ std(neighb_OFFmatrixx(:)*0);
% new_thresh = neighb_matrixx - neighb_thresh;
% matrixx = og_matrixx.*(new_thresh > 0);
% new_ONthresh = neighb_ONmatrixx - neighb_threshON;
% ONmatrixx = og_ONmatrixx.*(new_ONthresh>0);
% new_OFFthresh = neighb_OFFmatrixx- neighb_threshOFF;
% OFFmatrixx = og_OFFmatrixx.*(new_OFFthresh>0);


ON_flat_thresh = ONmatrixx>0;
OFF_flat_thresh = OFFmatrixx>0;
ON_flat_thresh = ON_flat_thresh(:);
OFF_flat_thresh = OFF_flat_thresh(:);
for i = 1: (n_rows*n_col)
    BOTHpass(i)= (ON_flat_thresh(i)+OFF_flat_thresh(i))>1;
    EITHERpass(i) = (ON_flat_thresh(i)+OFF_flat_thresh(i))>0;
end
ON_matrices(cell_n,:)= ONmatrixx(:);
OFF_matrices(cell_n,:) = OFFmatrixx(:);
OIindex1(cell_n) = sum(BOTHpass)/sum(EITHERpass);
    

        
    
    
%%    
  %ON_OFF ellipse
 
 [rf_cell(1,cell_n).amp.ONOFF, rf_cell(1,cell_n).sigmaX.ONOFF, rf_cell(1,cell_n).sigmaY.ONOFF,...
     rf_cell(1,cell_n).theta.ONOFF, rf_cell(1,cell_n).x0.ONOFF, rf_cell(1,cell_n).y0.ONOFF, ...
     rf_cell(1,cell_n).Earea.ONOFF, rf_cell(1,cell_n).axes_ratio.ONOFF, rf_cell(1,cell_n).azimuth.ONOFF, rf_cell(1,cell_n).elevation.ONOFF,...
     rf_cell(1,cell_n).rf_clust_ratio.ONOFF] = makeEllipse(matrixx, mean_spont(cell_n), threshold(cell_n));
     
%%    
%heatmap ON_OFF 
 R_data = matrixx; 
[Xq,Yq] = meshgrid(.5:.2:(n_col+0.5));
Vq = interp2(R_data,Xq,Yq);
Vq(isnan(Vq)) = 0;
rf_plot = Vq(1:(n_rows*5),1:(n_col*5)); 
    for i=1:round(n_rows*spacing)
        for j=1:round(n_col*spacing)
            index_x=floor((i-0.001)/spacing)+1;
            index_y=floor((j-0.001)/spacing)+1;
            if index_x==0
                index_x=1;
            end
            if index_y==0
                index_y=1;
            end
            rf_plot(i,j)=R_data(index_x,index_y);
        end
    end
    RF_figure = figure; imagesc(rf_plot); axis equal;
    % Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , 'ONOFF_map'],'-dpng','-r300');
%%    
 %plot ON_OFF ellipse
hold on
clear x 
clear y
 t = linspace(0,2*pi,1000);
 h = rf_cell(1,cell_n).x0.ONOFF;
 k = rf_cell(1,cell_n).y0.ONOFF;
 if rf_cell(1,cell_n).sigmaX.ONOFF > rf_cell(1,cell_n).sigmaY.ONOFF
     a = rf_cell(1,cell_n).sigmaX.ONOFF;
     b = rf_cell(1,cell_n).sigmaY.ONOFF;
 else
     b = rf_cell(1,cell_n).sigmaX.ONOFF;
     a = rf_cell(1,cell_n).sigmaY.ONOFF;
 end
 thet = rf_cell(1,cell_n).theta.ONOFF;
 try
     rf_cell(1,cell_n).Exs.ONOFF = h+ a * cos(t) * cos(thet) - b * sin(t) * sin(thet);
     rf_cell(1,cell_n).Eys.ONOFF = k+ b * sin(t) * cos(thet) + a * cos(t) * sin(thet) ;
     plot(rf_cell(1,cell_n).Exs.ONOFF, rf_cell(1,cell_n).Eys.ONOFF, 'w')
     axis([0 135 0 75]);
     set ( gca, 'ydir', 'reverse' )
     axis equal;
 end
 print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , 'ONOFF_ellipse'],'-dpng','-r300');
%% 
 %ON_only
  [rf_cell(1,cell_n).amp.ON, rf_cell(1,cell_n).sigmaX.ON, rf_cell(1,cell_n).sigmaY.ON,...
     rf_cell(1,cell_n).theta.ON, rf_cell(1,cell_n).x0.ON, rf_cell(1,cell_n).y0.ON, ...
     rf_cell(1,cell_n).Earea.ON, rf_cell(1,cell_n).axes_ratio.ON, rf_cell(1,cell_n).azimuth.ON, rf_cell(1,cell_n).elevation.ON,...
     rf_cell(1,cell_n).rf_clust_ratio.ON] = makeEllipse(ONmatrixx, mean_spont(cell_n), threshold(cell_n));
%% 
%heatmap ON
R_data = ONmatrixx; 
[Xq,Yq] = meshgrid(.5:.2:(n_col+0.5));
Vq = interp2(R_data,Xq,Yq);
Vq(isnan(Vq)) = 0;
rf_plot = Vq(1:(n_rows*5),1:(n_col*5)); 
    for i=1:round(n_rows*spacing)
        for j=1:round(n_col*spacing)
            index_x=floor((i-0.001)/spacing)+1;
            index_y=floor((j-0.001)/spacing)+1;
            if index_x==0
                index_x=1;
            end
            if index_y==0
                index_y=1;
            end
            rf_plot(i,j)=R_data(index_x,index_y);
        end
    end
    RF_figure = figure; imagesc(rf_plot); axis equal;
      colormap(RF_figure, redmap);
      colorbar;
      % Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , 'ON_Map'],'-dpng','-r300');
%%      
  %plot ON ellipse
hold on
clear x 
clear y
 t = linspace(0,2*pi,1000);
 h = rf_cell(1,cell_n).x0.ON;
 k = rf_cell(1,cell_n).y0.ON;
 if rf_cell(1,cell_n).sigmaX.ON > rf_cell(1,cell_n).sigmaY.ON
     a = rf_cell(1,cell_n).sigmaX.ON;
     b = rf_cell(1,cell_n).sigmaY.ON;
 else
     b = rf_cell(1,cell_n).sigmaX.ON;
     a = rf_cell(1,cell_n).sigmaY.ON;
 end
 thet = rf_cell(1,cell_n).theta.ON;
 try
     rf_cell(1,cell_n).Exs.ON = h+ a * cos(t) * cos(thet) - b * sin(t) * sin(thet);
     rf_cell(1,cell_n).Eys.ON = k+ b * sin(t) * cos(thet) + a * cos(t) * sin(thet) ;
     plot(rf_cell(1,cell_n).Exs.ON, rf_cell(1,cell_n).Eys.ON, 'w')
     axis([0 135 0 75]);
     set ( gca, 'ydir', 'reverse' )
     axis equal;
 end
 % Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , 'ON_ellipse'],'-dpng','-r300');
 %%
    %OFF_only
 
 [rf_cell(1,cell_n).amp.OFF, rf_cell(1,cell_n).sigmaX.OFF, rf_cell(1,cell_n).sigmaY.OFF,...
     rf_cell(1,cell_n).theta.OFF, rf_cell(1,cell_n).x0.OFF, rf_cell(1,cell_n).y0.OFF, ...
     rf_cell(1,cell_n).Earea.OFF, rf_cell(1,cell_n).axes_ratio.OFF, rf_cell(1,cell_n).azimuth.OFF, rf_cell(1,cell_n).elevation.OFF,... 
     rf_cell(1,cell_n).rf_clust_ratio.OFF] = makeEllipse(OFFmatrixx, mean_spont(cell_n), threshold(cell_n));
%%      
%heatmap OFF
R_data = OFFmatrixx; 
[Xq,Yq] = meshgrid(.5:.2:(n_col+0.5));
Vq = interp2(R_data,Xq,Yq);
Vq(isnan(Vq)) = 0;
rf_plot = Vq(1:(n_rows*5),1:(n_col*5)); 
    for i=1:round(n_rows*spacing)
        for j=1:round(n_col*spacing)
            index_x=floor((i-0.001)/spacing)+1;
            index_y=floor((j-0.001)/spacing)+1;
            if index_x==0
                index_x=1;
            end
            if index_y==0
                index_y=1;
            end
            rf_plot(i,j)=R_data(index_x,index_y);
        end
    end
    RF_figure = figure; imagesc(rf_plot); axis equal; 
      colormap(RF_figure,bluemap);
      colorbar;
 % Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , 'OFF_map'],'-dpng','-r300');
%%    
  %plot OFF ellipse
hold on
clear x 
clear y
 t = linspace(0,2*pi,1000);
 h = rf_cell(1,cell_n).x0.OFF;
 k = rf_cell(1,cell_n).y0.OFF;
 if rf_cell(1,cell_n).sigmaX.OFF > rf_cell(1,cell_n).sigmaY.OFF
     a = rf_cell(1,cell_n).sigmaX.OFF;
     b = rf_cell(1,cell_n).sigmaY.OFF;
 else
     b = rf_cell(1,cell_n).sigmaX.OFF;
     a = rf_cell(1,cell_n).sigmaY.OFF;
 end
 thet = rf_cell(1,cell_n).theta.OFF;
 try
     rf_cell(1,cell_n).Exs.OFF = h+ a * cos(t) * cos(thet) - b * sin(t) * sin(thet);
     rf_cell(1,cell_n).Eys.OFF = k+ b * sin(t) * cos(thet) + a * cos(t) * sin(thet) ;
     plot(rf_cell(1,cell_n).Exs.OFF, rf_cell(1,cell_n).Eys.OFF, 'w')
     axis([0 135 0 75]);
     set ( gca, 'ydir', 'reverse' )
     axis equal;
 end 
  % Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , 'OFF_ellipse'],'-dpng','-r300');
 %%
 %calculate ON OFF overlap 2 (ellipse overlap instead of grid)
 if (rf_cell(1,cell_n).x0.ON + rf_cell(1,cell_n).x0.OFF) > 0;
     
     ONx0 = rf_cell(1,cell_n).x0.ON;  OFFx0 = rf_cell(1,cell_n).x0.OFF;
     ONy0 = rf_cell(1,cell_n).y0.ON; OFFy0 = rf_cell(1,cell_n).y0.OFF;
     c = sqrt((ONx0- OFFx0)^2 + (ONy0- OFFy0)^2); %distance bt subfield centers
     m = (ONy0-OFFy0)/(ONx0-OFFx0);
     b = ONy0-(m*ONx0);
     hiONX = ceil(max(rf_cell(1,cell_n).Exs.ON));
     loONX = floor(min(rf_cell(1,cell_n).Exs.ON));
     hiOFFX = ceil(max(rf_cell(1,cell_n).Exs.OFF));
     loOFFX = floor(min(rf_cell(1,cell_n).Exs.OFF));
     ONlinexs = linspace(loONX,hiONX,500);
     ONlineys = m*ONlinexs + b;
     OFFlinexs = linspace(loOFFX,hiOFFX,500);
     OFFlineys = m*OFFlinexs + b;     
     alldistancesON = zeros(1000,500);
     alldistancesOFF = zeros(1000,500);
     tic
     for e = 1:1000
         for l = 1:500
             alldistancesON(e,l) = sqrt((ONlinexs(l) - rf_cell(1,cell_n).Exs.ON(e))^2 + (ONlineys(l) - rf_cell(1,cell_n).Eys.ON(e))^2);
             alldistancesOFF(e,l) = sqrt((OFFlinexs(l) - rf_cell(1,cell_n).Exs.OFF(e))^2 + (OFFlineys(l) - rf_cell(1,cell_n).Eys.OFF(e))^2);
         end
     end
     toc
     [minD, bestLON]= min(min(alldistancesON));
     %bestEON = find(alldistancesON(:,bestLON)== minD);
     [minD, bestLOFF] = min(min(alldistancesOFF));
     %bestEOFF = find(alldistancesOFF(:,bestLOFF)== minD);
     w1 = sqrt((ONx0- ONlinexs(bestLON))^2 + (ONy0- ONlineys(bestLON))^2);
     w2 = sqrt((OFFx0- OFFlinexs(bestLOFF))^2 + (OFFy0- OFFlineys(bestLOFF))^2);
     OI_index2(cell_n) = (w1+w2 -c)/(w1+w2+c);
 else
     OI_index2(cell_n) = 0;
 end
 
 response_ratio(cell_n) = (rf_cell(1,cell_n).amp.ON-rf_cell(1,cell_n).amp.OFF)/(rf_cell(1,cell_n).amp.ON+rf_cell(1,cell_n).amp.OFF);
 area_ratio(cell_n) = (rf_cell(1,cell_n).Earea.ON- rf_cell(1,cell_n).Earea.OFF)/(rf_cell(1,cell_n).Earea.ON+ rf_cell(1,cell_n).Earea.OFF);
 
    excel_variables{1,cell_n} = channel_no;
    excel_variables{2,cell_n} = clust_no;
    excel_variables{3,cell_n} = rf_cell(1,cell_n).Earea.ONOFF;
    excel_variables{4,cell_n} = rf_cell(1,cell_n).Earea.ON;
    excel_variables{5,cell_n} = rf_cell(1,cell_n).Earea.OFF;
    excel_variables{6,cell_n} = rf_cell(1,cell_n).theta.ONOFF;
    excel_variables{7,cell_n} = rf_cell(1,cell_n).theta.ON;
   excel_variables{8,cell_n} = rf_cell(1,cell_n).theta.OFF;
    excel_variables{9,cell_n} = rf_cell(1,cell_n).x0.ONOFF;
    excel_variables{10,cell_n} = rf_cell(1,cell_n).x0.ON;
    excel_variables{11,cell_n} = rf_cell(1,cell_n).x0.OFF;
    excel_variables{12,cell_n} = rf_cell(1,cell_n).y0.ONOFF;
    excel_variables{13,cell_n} = rf_cell(1,cell_n).y0.ON;
    excel_variables{14,cell_n} = rf_cell(1,cell_n).y0.OFF;
    excel_variables{15,cell_n} = rf_cell(1,cell_n).sigmaX.ONOFF;
    excel_variables{16,cell_n} = rf_cell(1,cell_n).sigmaX.ON;
    excel_variables{17,cell_n} = rf_cell(1,cell_n).sigmaX.OFF;
    excel_variables{18,cell_n} = rf_cell(1,cell_n).sigmaY.ONOFF;
    excel_variables{19,cell_n} = rf_cell(1,cell_n).sigmaY.ON;
    excel_variables{20,cell_n} = rf_cell(1,cell_n).sigmaY.OFF;
    excel_variables{21,cell_n} = rf_cell(1,cell_n).axes_ratio.ONOFF;
    excel_variables{22,cell_n} = rf_cell(1,cell_n).axes_ratio.ON;
    excel_variables{23,cell_n} = rf_cell(1,cell_n).axes_ratio.OFF;
    excel_variables{24,cell_n} = rf_cell(1,cell_n).azimuth.ONOFF;
    excel_variables{25,cell_n} = rf_cell(1,cell_n).azimuth.ON;
    excel_variables{26,cell_n} = rf_cell(1,cell_n).azimuth.OFF;   
    excel_variables{27,cell_n} = rf_cell(1,cell_n).elevation.ONOFF;
    excel_variables{28,cell_n} = rf_cell(1,cell_n).elevation.ON;
    excel_variables{29,cell_n} = rf_cell(1,cell_n).elevation.OFF; 
    excel_variables{30,cell_n} = rf_cell(1,cell_n).rf_clust_ratio.ONOFF;
    excel_variables{31,cell_n} = rf_cell(1,cell_n).rf_clust_ratio.ON;
    excel_variables{32,cell_n} = rf_cell(1,cell_n).rf_clust_ratio.OFF;   
    excel_variables{33,cell_n} = rf_cell(1,cell_n).amp.ONOFF;
    excel_variables{34,cell_n} = rf_cell(1,cell_n).amp.ON;
    excel_variables{35,cell_n} = rf_cell(1,cell_n).amp.OFF;     
    excel_variables{36,cell_n} = mean_spont(cell_n);
    excel_variables{37,cell_n} = OIindex1(cell_n);
    excel_variables{38,cell_n} = OI_index2(cell_n);
    excel_variables{39,cell_n} = correlation_coef(cell_n);
    excel_variables{40,cell_n} = response_ratio(cell_n);
    excel_variables{41,cell_n} = area_ratio(cell_n);
    excel_variables{42,cell_n} = rf_clust_ratio(cell_n); %raw clustering ratio
    excel_variables{43,cell_n} = rf_cell(1,cell_n).azimuth.ONOFF * rf_cell(1,cell_n).elevation.ONOFF * pi;
    excel_variables{44,cell_n} = rf_cell(1,cell_n).azimuth.ON * rf_cell(1,cell_n).elevation.ON * pi;
    excel_variables{45,cell_n} = rf_cell(1,cell_n).azimuth.OFF * rf_cell(1,cell_n).elevation.OFF * pi;
    excel_variables{46,cell_n} = (excel_variables{44,cell_n}- excel_variables{45,cell_n})/(excel_variables{44,cell_n}+ excel_variables{45,cell_n});
end  %%% cell

%save(AnalysisFileName, 'Response')
    

