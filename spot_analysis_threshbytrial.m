%finds receptive fields but thresholds each spike rate by trial
%by Rachel Kay 1/12/16
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

gridX = 27; % number of X positions on grid
gridY = 15; % number of Y positions on grid
n_col = 27;
n_rows = 15;
duration = 0.5; %500 ms
waitInt = 0.2; %not used
blankframe=1; %1 = yes

nCond = gridX*gridY;
if blankframe
    nCond = nCond+1;
end

%xTrg can only hold up to 255, made extra fTrg of two to indicate stimcond+
%255
StimOnset = stimulus.epocs.xTrg.onset;
fTrg = stimulus.epocs.fTrg.data;
high_cond = 0;
if max(fTrg) > 1
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
end

matrixx = zeros(15,27);
ONmatrixx = zeros(15,27);
OFFmatrixx = zeros(15,27);
n_cells = length(spikeT);
excel_variables = cell(12, n_cells);
for m = 1: n_cells
    rf_cell(m).condition(1).trial{1}={};
    rf_cell(m).condition(:).all_trials = [];
    rf_cell(m).condition(:).ONtrials{1} = [];
    thresholded(m).condition(:)= [];
end 
%fill rf.cells with spiketimes by condition, ON OFF
ONgrid = zeros(406,n_cells);
OFFgrid = zeros(406,n_cells);
ORgrid = zeros(406,n_cells);
ONyes_no = zeros(406, n_cells);
OFFyes_no = zeros(406, n_cells);
BOTHyes_no = zeros(406, n_cells);
EITHERyes_no = zeros(406, n_cells);
for i = 1:n_cells
    cellspike = spikeT{i};
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
        thresh1 =  cellspike(cellspike > (startT - 0.2) & cellspike < startT);
        thresh2 =  cellspike(cellspike > stopT & cellspike < (stopT + 0.2));
        thresh = ((length(thresh1) + length(thresh2))/2)* 1.2;
        rf_cell(i).condition(StimCond(j)).ONcount{trial_no} = length(on_spikes);
        rf_cell(i).condition(StimCond(j)).OFFcount{trial_no} = length(off_spikes);
        rf_cell(i).condition(StimCond(j)).thresh{trial_no} = thresh;
        rf_cell(i).condition(StimCond(j)).ONthreshed{trial_no} = length(on_spikes)- thresh;
        rf_cell(i).condition(StimCond(j)).OFFthreshed{trial_no} = length(off_spikes)- thresh;
        rf_cell(i).condition(StimCond(j)).ON_OFFthreshed{trial_no} = ((length(on_spikes)+ length(off_spikes))/2) - thresh;
        rf_cell(i).condition(StimCond(j)).ONpass{trial_no} = (length(on_spikes)- thresh) > 0 ;
        rf_cell(i).condition(StimCond(j)).OFFpass{trial_no} = (length(off_spikes)- thresh) > 0 ;
        rf_cell(i).condition(StimCond(j)).ON_OFFpass{trial_no} = (((length(on_spikes)+ length(off_spikes))/2) - thresh) > 0 ;
        trialpercond(StimCond(j))= trial_no;
    end
    %use all trials for thresholding 
    for k = 1:405 %stim condition
        if ((mean(cell2mat(rf_cell(1,i).condition(1,k).ONpass))) >= 0.4)
            avgONspike = cell2mat(rf_cell(1,i).condition(1,k).ONthreshed);
            avgONspike(avgONspike<0) = 0;
            avgONspike = mean(avgONspike);
            ONyes_no(k,i) = 1;
        else avgONspike = 0; ONyes_no(k,i) = 0;
        end
        if ((mean(cell2mat(rf_cell(1,i).condition(1,k).OFFpass))) >= 0.4)
            avgOFFspike = cell2mat(rf_cell(1,i).condition(1,k).OFFthreshed);
            avgOFFspike(avgOFFspike<0) = 0;
            avgOFFspike = mean(avgOFFspike);
            OFFyes_no(k,i) = 1;
        else avgOFFspike = 0; OFFyes_no(k,i) = 0;
        end
     ONgrid(k,i) = avgONspike;
     OFFgrid(k,i) = avgOFFspike;
     ORgrid(k,i) = (avgONspike + avgOFFspike)/2;
        if ONyes_no(k,i)+ OFFyes_no(k,i) == 2
            BOTHyes_no(k,i)= 1;
        end
        if ONyes_no(k,i) + OFFyes_no(k,i) >0
                EITHERyes_no(k,i) = 1;
        end
     
    end
end


%plot histograms
binranges = linspace(0,1,11); %for 100ms bins= 11, change 11 to 21 for 50ms bins etc
%plotinx = [1,18,35,52,69,86,103,120,137,154,171,188,205,2,19,36,53,70,87,104,121,138,155,172,189,206,3,20,37,54,71,88,105,122,139,156,173,190,207,4,21,38,55,72,89,106,123,140,157,174,191,208,5,22,39,56,73,90,107,124,141,158,175,192,209,6,23,40,57,74,91,108,125,142,159,176,193,210,7,24,41,58,75,92,109,126,143,160,177,194,211,8,25,42,59,76,93,110,127,144,161,178,195,212,9,26,43,60,77,94,111,128,145,162,179,196,213,10,27,44,61,78,95,112,129,146,163,180,197,214,11,28,45,62,79,96,113,130,147,164,181,198,215,12,29,46,63,80,97,114,131,148,165,182,199,216,13,30,47,64,81,98,115,132,149,166,183,200,217,14,31,48,65,82,99,116,133,150,167,184,201,218,15,32,49,66,83,100,117,134,151,168,185,202,219,16,33,50,67,84,101,118,135,152,169,186,203,220,17,34,51,68,85,102,119,136,153,170,187,204,221];
plotinx =  [1,16,31,46,61,76,91,106,121,136,151,166,181,196,211,226,241,256,271,286,301,316,331,346,361,376,391,2,17,32,47,62,77,92,107,122,137,152,167,182,197,212,227,242,257,272,287,302,317,332,347,362,377,392,3,18,33,48,63,78,93,108,123,138,153,168,183,198,213,228,243,258,273,288,303,318,333,348,363,378,393,4,19,34,49,64,79,94,109,124,139,154,169,184,199,214,229,244,259,274,289,304,319,334,349,364,379,394,5,20,35,50,65,80,95,110,125,140,155,170,185,200,215,230,245,260,275,290,305,320,335,350,365,380,395,6,21,36,51,66,81,96,111,126,141,156,171,186,201,216,231,246,261,276,291,306,321,336,351,366,381,396,7,22,37,52,67,82,97,112,127,142,157,172,187,202,217,232,247,262,277,292,307,322,337,352,367,382,397,8,23,38,53,68,83,98,113,128,143,158,173,188,203,218,233,248,263,278,293,308,323,338,353,368,383,398,9,24,39,54,69,84,99,114,129,144,159,174,189,204,219,234,249,264,279,294,309,324,339,354,369,384,399,10,25,40,55,70,85,100,115,130,145,160,175,190,205,220,235,250,265,280,295,310,325,340,355,370,385,400,11,26,41,56,71,86,101,116,131,146,161,176,191,206,221,236,251,266,281,296,311,326,341,356,371,386,401,12,27,42,57,72,87,102,117,132,147,162,177,192,207,222,237,252,267,282,297,312,327,342,357,372,387,402,13,28,43,58,73,88,103,118,133,148,163,178,193,208,223,238,253,268,283,298,313,328,343,358,373,388,403,14,29,44,59,74,89,104,119,134,149,164,179,194,209,224,239,254,269,284,299,314,329,344,359,374,389,404,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,360,375,390,405];
%plotinx2= [213,196,179,162,145,128,111,94,77,60,43,26,9,212,195,178,161,144,127,110,93,76,59,42,25,8,211,194,177,160,143,126,109,92,75,58,41,24,7,210,193,176,159,142,125,108,91,74,57,40,23,6,209,192,175,158,141,124,107,90,73,56,39,22,5,208,191,174,157,140,123,106,89,72,55,38,21,4,207,190,173,156,139,122,105,88,71,54,37,20,3,206,189,172,155,138,121,104,87,70,53,36,19,2,205,188,171,154,137,120,103,86,69,52,35,18,1,221,204,187,170,153,136,119,102,85,68,51,34,17,220,203,186,169,152,135,118,101,84,67,50,33,16,219,202,185,168,151,134,117,100,83,66,49,32,15,218,201,184,167,150,133,116,99,82,65,48,31,14,217,200,183,166,149,132,115,98,81,64,47,30,13,216,199,182,165,148,131,114,97,80,63,46,29,12,215,198,181,164,147,130,113,96,79,62,45,28,11,214,197,180,163,146,129,112,95,78,61,44,27,10];
%eqt = zeros(n_cells,6);

for cell_n = 1:n_cells
    channel_no = cells(cell_n,1);
    clust_no = cells(cell_n,2);
    figure
    bincounts= zeros(406,11);%10 bins
    for j = 1: nCond
        bincounts(j,:) = histc(rf_cell(1,cell_n).condition(1,j).all_trials,binranges);
    end
    ymax = max(max(bincounts));
    for j = 1:(nCond-1)%dont plot blank cond
        subplot(15,27,plotinx(j))
        bar(binranges,bincounts(j,:),'histc')
        axis([0 1 0 ymax])
        set(gca,'XTicklabel',[])
        set(gca,'YTicklabel',[])
    end
% Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '1histogram'],'-dpng','-r300');
%plot heatmap
%xcontour= [1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17];
%ycontour = [13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1,13,12,11,10,9,8,7,6,5,4,3,2,1];
xcontour= [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,22,22,22,22,22,22,22,22,22,22,22,22,22,22,22,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,24,24,24,24,24,24,24,24,24,24,24,24,24,24,24,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27];
ycontour= [15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1];
figure 
    for p = 1: 405
        matrixx(ycontour(p),xcontour(p))=ORgrid(p,cell_n);
        ONmatrixx(ycontour(p),xcontour(p))=ONgrid(p,cell_n);
        OFFmatrixx(ycontour(p),xcontour(p))=OFFgrid(p,cell_n);
    end
    contourf(matrixx);
% Save the file as PNG
    print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '2contourplot'],'-dpng','-r300');

%Either On or Off RF    
    % prepare initial parameters for JC's code PSA_AnalyzeSpots
    temp = reshape(matrixx,405,1);
    rfthresh = matrixx;
    spacing = 5;
    %spontmean = mean(cell2mat(rf_cell(1,cell_n).condition(1,222).ONcount))+ mean(cell2mat(rf_cell(1,cell_n).condition(1,222).OFFcount));
    spontmean = 0;
    RF_figure = figure;
    note_text=sprintf('%s %s Channel%d Cluster%d ',Tank_Name,Block_Name, channel_no, clust_no);
    axis equal; colorbar; title(note_text); 
    rf_data = reshape(temp(1:n_col*n_rows), n_rows, n_col);

    % estimating center and dispersion along x axis
    n= 0 ; nsx = 0; nsx2 = 0;
    for i = 1:n_col
        ns = sum(rfthresh(:, i));
        n = n+ns;
        nsx = nsx+ns*i;
        nsx2 = nsx2+ns*i*i;
    end;
    centroid_x(cell_n) = nsx/n;
    dispersion_x(cell_n) = sqrt((nsx2/n) - (nsx/n)^2);

    % estimating center and dispersion along y axis
    n = 0 ; nsy = 0; nsy2 = 0;
    for i = 1:n_rows
        ns = sum(rfthresh(i, :));
        n = n+ns;
        nsy = nsy+ns*i;
        nsy2 = nsy2+ns*i*i;
    end;
    centroid_y(cell_n) = nsy/n;
    dispersion_y(cell_n) = sqrt((nsy2/n) - (nsy/n)^2);
        
    obs = rf_data(:); %reduce to one dimension
    [xgrid ygrid]= meshgrid(1:size(rf_data,2),1:size(rf_data,1));
    x(:,1) = xgrid(:);
    x(:,2) = ygrid(:);
    try 
        coeff = nlinfit(x,obs,@rf2fit,[spontmean max(obs)-spontmean ...
        centroid_x(cell_n) dispersion_x(cell_n) centroid_y(cell_n) dispersion_y(cell_n)]);
        est = rf2fit(coeff,x);
        figure(RF_figure); imagesc(reshape(est,[size(rf_data,1) size(rf_data,2)]));
        set(gca,'YDir','normal')
        axis equal; colorbar; title('rf fit');
        rf_amp = coeff(2);
        spont = coeff(1);
        x0 = coeff(3)*spacing;
        wx0 = abs(coeff(4))*spacing;
        y0 = coeff(5)*spacing;
        wy0 = abs(coeff(6))*spacing;
        e_area = pi*wx0*wy0*5;
        est = 1000.*est;
        est = (round(est)/1000);
        est(est == mode(est)) = 0;
    catch
        rf_amp = 'error'; spont = 'error'; x0= 'error'; wx0 = 'error'; y0 = 'error'; wy0 = 'error'; e_area = 'error'; est = 'error'; fitted = temp;
    end
 % Save the file as PNG
     print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '3RF_figure'],'-dpng','-r300');  
    
    
    % RF fit for ON Response
    
    
    % prepare initial parameters for JC's code PSA_AnalyzeSpots
    ONtemp = reshape(ONmatrixx,405,1);
    ONrfthresh = ONmatrixx;
    spacing = 5;
    %spontmean = mean(cell2mat(rf_cell(1,cell_n).condition(1,222).ONcount))+ mean(cell2mat(rf_cell(1,cell_n).condition(1,222).OFFcount));
    spontmean = 0.001;
    RF_figure = figure;
    note_text=sprintf('%s %s Channel%d Cluster%d ',Tank_Name,Block_Name, channel_no, clust_no);
    axis equal; colormap(summer); title(note_text); 
    rf_data = reshape(ONtemp(1:n_col*n_rows), n_rows, n_col);

    % estimating center and dispersion along x axis
    n= 0 ; nsx = 0; nsx2 = 0;
    for i = 1:n_col
        ns = sum(ONrfthresh(:, i));
        n = n+ns;
        nsx = nsx+ns*i;
        nsx2 = nsx2+ns*i*i;
    end;
    centroid_x(cell_n) = nsx/n;
    dispersion_x(cell_n) = sqrt((nsx2/n) - (nsx/n)^2);

    % estimating center and dispersion along y axis
    n = 0 ; nsy = 0; nsy2 = 0;
    for i = 1:n_rows
        ns = sum(ONrfthresh(i, :));
        n = n+ns;
        nsy = nsy+ns*i;
        nsy2 = nsy2+ns*i*i;
    end;
    centroid_y(cell_n) = nsy/n;
    dispersion_y(cell_n) = sqrt((nsy2/n) - (nsy/n)^2);
        
    obs = rf_data(:); %reduce to one dimension
    [xgrid ygrid]= meshgrid(1:size(rf_data,2),1:size(rf_data,1));
    x(:,1) = xgrid(:);
    x(:,2) = ygrid(:);
    try
    coeff = nlinfit(x,obs,@rf2fit,[spontmean max(obs)-spontmean ...
        centroid_x(cell_n) dispersion_x(cell_n) centroid_y(cell_n) dispersion_y(cell_n)]);
    ON_est = rf2fit(coeff,x);
    figure(RF_figure); imagesc(reshape(ON_est,[size(rf_data,1) size(rf_data,2)]));
    set(gca,'YDir','normal')
    axis equal; colormap(summer); colorbar; title('rf fit');
    ONrf_amp = coeff(2);
    spont = coeff(1);
    ONx0 = coeff(3)*spacing;
    ONwx0 = abs(coeff(4))*spacing;
    ONy0 = coeff(5)*spacing;
    ONwy0 = abs(coeff(6))*spacing;
    ONe_area = pi*ONwx0*ONwy0*5;
    ON_est = 1000.*ON_est;
    ON_est = (round(ON_est)/1000);
    ON_est(ON_est == mode(ON_est)) = 0;
    catch
        ONrf_amp = 'error'; spont = 'error'; ONx0= 'error'; ONwx0 = 'error'; ...
            ONy0 = 'error'; ONwy0 = 'error'; ONe_area = 'error'; ON_est = 'error'; ...
            ONfitted = ONtemp;
    end
    
   % note_text = sprintf('RF fit: Spont %0.2g, Amp %0.2g, X0 %0.2gdeg, stdX %0.2gdeg, Y0 %0.2gdeg, stdY %0.2gdeg', ...
    %    spont, ONrf_amp, ONx0, ONwx0, ONy0, ONwy0);
%     figure(RF_figure);
   % axes('Position',[0 0 1 1],'Visible','off'); 
    %text(.025,0.025,note_text,'FontSize',10);

% Save the file as PNG
      print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '4RF_on'],'-dpng','-r300');
      
    % RF fit for OFF Response
    

    % prepare initial parameters for JC's code PSA_AnalyzeSpots
    OFFtemp = reshape(OFFmatrixx,405,1);
    OFFrfthresh = OFFmatrixx;
    spacing = 5;
    %spontmean = mean(cell2mat(rf_cell(1,cell_n).condition(1,222).ONcount))+ mean(cell2mat(rf_cell(1,cell_n).condition(1,222).OFFcount));
    spontmean = 0;
    RF_figure = figure;
    note_text=sprintf('%s %s Channel%d Cluster%d ',Tank_Name,Block_Name, channel_no, clust_no);
    axis equal; colormap(autumn); colorbar; title(note_text); 
    rf_data = reshape(OFFtemp(1:n_col*n_rows), n_rows, n_col);

    % estimating center and dispersion along x axis
    n= 0 ; nsx = 0; nsx2 = 0;
    for i = 1:n_col
        ns = sum(OFFrfthresh(:, i));
        n = n+ns;
        nsx = nsx+ns*i;
        nsx2 = nsx2+ns*i*i;
    end;
    centroid_x(cell_n) = nsx/n;
    dispersion_x(cell_n) = sqrt((nsx2/n) - (nsx/n)^2);

    % estimating center and dispersion along y axis
    n = 0 ; nsy = 0; nsy2 = 0;
    for i = 1:n_rows
        ns = sum(OFFrfthresh(i, :));
        n = n+ns;
        nsy = nsy+ns*i;
        nsy2 = nsy2+ns*i*i;
    end;
    centroid_y(cell_n) = nsy/n;
    dispersion_y(cell_n) = sqrt((nsy2/n) - (nsy/n)^2);
        
    obs = rf_data(:); %reduce to one dimension
    [xgrid ygrid]= meshgrid(1:size(rf_data,2),1:size(rf_data,1));
    x(:,1) = xgrid(:);
    x(:,2) = ygrid(:);
    try 
        coeff = nlinfit(x,obs,@rf2fit,[spontmean max(obs)-spontmean ...
        centroid_x(cell_n) dispersion_x(cell_n) centroid_y(cell_n) dispersion_y(cell_n)]);
        OFF_est = rf2fit(coeff,x);
        figure(RF_figure); imagesc(reshape(OFF_est,[size(rf_data,1) size(rf_data,2)]));
        set(gca,'YDir','normal')
        axis equal; colorbar; colormap(autumn); title('rf fit');
        OFFrf_amp = coeff(2);
        spont = coeff(1);
        OFFx0 = coeff(3)*spacing;
        OFFwx0 = abs(coeff(4))*spacing;
        OFFy0 = coeff(5)*spacing;
        OFFwy0 = abs(coeff(6))*spacing;
        OFFe_area = pi*OFFwx0*OFFwy0*5;
        OFF_est = 1000.*OFF_est;
        OFF_est = (round(OFF_est)/1000);
        OFF_est(OFF_est == mode(OFF_est)) = 0;
    catch
        OFFrf_amp = 'error'; spont = 'error'; OFFx0= 'error'; OFFwx0 = 'error'; OFFy0 = 'error'; OFFwy0 = 'error'; OFFe_area = 'error'; OFF_est = 'error'; OFFfitted = OFFtemp;
    end
 % Save the file as PNG
     print([Block_Name,'_ch_no',num2str(channel_no),'_clust_no', num2str(clust_no), '_' , '5RF_off'],'-dpng','-r300');
      

 %ON_OFF overlap
    %ON_OFF_overlap_ratios 2= not fit, 1 = from fit
    ON_OFF_overlap_ratio2 = sum(BOTHyes_no(:,cell_n))/sum(EITHERyes_no(:,cell_n));
    
    ON_fit_grid = ON_est;
    ON_est(ON_est > 0) = 1;
    OFF_fit_grid = OFF_est;
    OFF_est(OFF_est > 0) = 1;
    try
    for k= 1:221
    if ON_est(k,1)+ OFF_est(k,1) == 2
            BOTHy_n_fit(k,cell_n)= 1;
    end
    if ON_est(k,1) + OFF_est(k,1) >0
                EITHERy_n_fit(k,cell_n) = 1;
    end
    end
    ON_OFF_overlap_ratio1 = sum(BOTHy_n_fit(:,cell_n))/sum(EITHERy_n_fit(:,cell_n));
    catch 
     ON_OFF_overlap_ratio1 = 'error';
    end
    
    
    %response ratio
    try response_ratio = (max(ONgrid(:, cell_n))-max(OFFgrid(:, cell_n)))/(max(ONgrid(:, cell_n))+max(OFFgrid(:, cell_n)));
    catch
        response_ratio = 'error';
    end
    %area ratio
    try area_ratio = (ONe_area- OFFe_area)/(ONe_area + OFFe_area);
    catch
        area_ratio = 'error';
    end
    
    
    Response(cell_n) = struct('Tank',Tank_Name, 'Block', Block_Name, 'ch', channel_no, 'cluster', clust_no, ...
        'rf_data', rf_data, 'Spont', spontmean, 'ON_Amp_Fit', ONrf_amp, 'Spont_Fit', spont, ...
        'ON_X0_Fit', ONx0, 'ON_wX0_Fit', ONwx0, 'ON_Y0_Fit', ONy0, 'ON_wY0_Fit', ONwy0, 'ON_est', ON_est,...
        'ON_e_area', ONe_area,  ...
         'OFF_Amp_Fit', OFFrf_amp, 'OFF_X0_Fit', OFFx0, 'OFF_wX0_Fit', OFFwx0, 'OFF_Y0_Fit', OFFy0,...
         'OFF_wY0_Fit', OFFwy0, 'OFF_est', OFF_est, 'OFF_e_area', OFFe_area, 'ON_OFF_overlap_ratio1', ON_OFF_overlap_ratio1, 'matrix', matrixx, 'ONmatrix', ONmatrixx, 'OFFmatrix', OFFmatrixx, ...
         'ON_OFF_overlap_ratio2', ON_OFF_overlap_ratio2, 'response_ratio', response_ratio, 'area_ratio', area_ratio);

    excel_variables{1,cell_n} = channel_no;
    excel_variables{2,cell_n} = clust_no;
    excel_variables{3,cell_n} = e_area;
    excel_variables{4,cell_n} = ONe_area;
    excel_variables{5,cell_n} = OFFe_area;
    excel_variables{6,cell_n} = ON_OFF_overlap_ratio1;
    excel_variables{7,cell_n} = ON_OFF_overlap_ratio2;
    excel_variables{8,cell_n} = rf_amp;
    excel_variables{9,cell_n} = ONrf_amp;
    excel_variables{10,cell_n} = OFFrf_amp;
    excel_variables{11,cell_n} = response_ratio;
    excel_variables{12,cell_n} = area_ratio;

end  %%% cell

save(AnalysisFileName, 'Response')
    
   
%best-fit ellipse
%     indxALL = find(ALLgrid(:,i) > 0);
%     for b = 1:length(indxALL)
%         e_points(1,b)= xcontour(indxALL(b))*5;
%         e_points(2,b)= ycontour(indxALL(b))*5;
%     end
%     [rf_cell(i).ALLequation , rf_cell(i).ALLarea, rf_cell(i).ALLx_center, rf_cell(i).ALLy_center] = EllipseFitRachel(e_points);
%     indxON = find(ONgrid(:,i) > 0);
%     for b = 1:length(indxON)
%         e_points(1,b)= xcontour(indxON(b))*5;
%         e_points(2,b)= ycontour(indxON(b))*5;
%     end
%     [rf_cell(i).ONequation , rf_cell(i).ONarea, rf_cell(i).ONx_center, rf_cell(i).ONy_center] = EllipseFitRachel(e_points);
%     indxOFF = find(OFFgrid(:,i) > 0);
%     for b = 1:length(indxOFF)
%         e_points(1,b)= xcontour(indxOFF(b))*5;
%         e_points(2,b)= ycontour(indxOFF(b))*5;
%     end
%     [rf_cell(i).OFFequation, rf_cell(i).OFFarea, rf_cell(i).OFFx_center, rf_cell(i).OFFy_center] = EllipseFitRachel(e_points);
%     
%     figure  
%     hold on
% p1 = ezplot(rf_cell(i).ALLequation, [0 85 0 65]);
% set(p1,'Color','black','LineWidth', 2)
% p2 = ezplot(rf_cell(i).ONequation, [0 85 0 65]);
% set(p2,'Color','green')
% p3 = ezplot(rf_cell(i).OFFequation, [0 85 0 65]);
% set(p3,'Color','red')
% title(['Cell ',num2str(i)])
% hold off

%     A = EllipseDirectFit(e_points');
%     %Equation 
%     eqt(i,1)= A(1);
%     eqt(i,2)= A(2);
%     eqt(i,3)= A(3);
%     eqt(i,4)= A(4);
%     eqt(i,5)=A(5);
%     eqt(i,6)=A(6);
%     
%     %Convert the A to str 
%     aa = num2str(A(1)); 
%     bb = num2str(A(2)); 
%     cc = num2str(A(3)); 
%     dd = num2str(A(4)); 
%     ee = num2str(A(5)); 
%     ff = num2str(A(6));
%     
%     eqtion= ['(',aa, ')*x^2 + (',bb,')*x*y + (',cc,')*y^2 + (',dd,')*x+ (',ee,')*y + (',ff,')'];
%     xmin=0; 
%     xmax=17; 
%     ezplot(eqtion,[xmin,xmax])
%     %scatter(e_points(1,:),e_points(2,:)) 
%     
%     a = A(1);
%     b = A(2)/2;
%     c = A(3);
%     d = A(4)/2;
%     f = A(5)/2;
%     g = A(6);
%     
%     
%     aprime = sqrt(((2*((a*f^2)+(c*d^2)+(g*b^2)-(2*b*d*f)-(a*c*g)))/((b^2-(a*c))*((sqrt(a-c)^2)+(4*b^2))-(a+c))));
%     bprime = sqrt(((2*((a*f^2)+(c*d^2)+(g*b^2)-(2*b*d*f)-(a*c*g)))/((b^2-(a*c))*((-1)*(sqrt(a-c)^2)+(4*b^2))-(a+c))));
%     areaaaa(i)= pi*(aprime/2)*(bprime/2);
%     
%     %coefficient normalizing factor
%     %q = ((64*f)*(4*a*c - b^2)-(64*((a*e^2) - (b*d*e) - (c*d^2))))/ ((4*a*c-b^2)^2);
%     
%     q = (((64*g)*((4*a*c)-b^2))- (64*((a*f^2)-(b*d*f)-(c*d^2))))/(((4*a*c)-b^2)^2);
%     %distance between center and focal point
%     s = 0.25*(sqrt(abs(q)*sqrt((b^2)+((a-c)^2))));
%     %semi-major axis length
%     %rmax = 0.125*(sqrt(2*abs(q)*sqrt(b^2+((a-c)^2))-(2*q*(a+c))));
%     rmax= 0.125*(sqrt((2*abs(q))*sqrt((b^2)+((a-c)^2))-((2*q)*(a+c))));
%     
%     %semi-minor axis length
%     %rmin = sqrt((rmax^2)-(.25*(sqrt(abs(q)*sqrt(b^2+((a-c)^2)))))^2);
%     rmin = sqrt((rmax^2)-(s^2));
%     %center of ellipse
%     elliX(i) = (b*f-2*c*d)/((4*a*c)-b^2);
%     elliY(i) = (b*d-2*a*f)/((4*a*c)-b^2);
%     %area
%     Area(i) = pi*(rmax/2)*(rmin/2);
% 
 

 
 
 