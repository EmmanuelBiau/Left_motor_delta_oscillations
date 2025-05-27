%% source grandaverage across participants;
clearvars; clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd xxx\headmodel\; load standard_mri;
cd xxx\source_analyses\headmodel\; load('grid.mat');

subjects = [1 2 4:11 13:19 21:26];

%Create a structure with all the MI_source data together;
clear allsubj_source

for s = subjects

    %load participant's source data;
    cd(['xxx\individual_source_data\subj_',num2str(s)]);
    load(['source_data_subj',num2str(s)]);
    
    %Select the time-window and frequency band of interest;
    cfg = [];
    cfg.latency     = [3 9];
    cfg.avgovertime = 'yes';
    cfg.frequency   = [2 3];
    cfg.avgoverfreq = 'yes';
    NMAC_source = ft_selectdata(cfg,source_NMAC);
    NMSC_source = ft_selectdata(cfg,source_NMSC);
    HMAC_source = ft_selectdata(cfg,source_HMAC);
    HMSC_source = ft_selectdata(cfg,source_HMSC);
    
    %Put all participants in a common structure;
    allsubj_source.NMAC{1,s} = NMAC_source;
    allsubj_source.NMSC{1,s} = NMSC_source;
    allsubj_source.HMAC{1,s} = HMAC_source;
    allsubj_source.HMSC{1,s} = HMSC_source;

end

clear allsource_NMCs allsource_HMCs
allsource_NMCs = cell(1,size(allsubj_source.NMAC,2));
allsource_HMCs = cell(1,size(allsubj_source.HMAC,2));

for s = subjects

    diff_NMCs_temp = allsubj_source.NMAC{s};
    diff_NMCs_temp.powspctrm = [];
    diff_NMCs_temp.powspctrm = (allsubj_source.NMAC{s}.powspctrm - allsubj_source.NMSC{s}.powspctrm)./allsubj_source.NMSC{s}.powspctrm;  
    allsource_NMCs{s} = diff_NMCs_temp;
    
    diff_HMCs_temp = allsubj_source.HMAC{s};
    diff_HMCs_temp.powspctrm = [];
    diff_HMCs_temp.powspctrm = (allsubj_source.HMAC{s}.powspctrm - allsubj_source.HMSC{s}.powspctrm)./allsubj_source.HMSC{s}.powspctrm;   
    allsource_HMCs{s} = diff_HMCs_temp;
    
    clear diff_NMCs_temp diff_HMCs_temp
    
end

%% Extract delta power difference Async-Sync at grids of interest then correlat with behaviors;

% load grids within the significant cluster of Left Motor and behavioural scores;
load allsubj_behavior_scores; load index_grid_Left_Motor; 

clear allsubj_grid_delta

for s = subjects

    cfg = [];
    cfg.parameter = 'powspctrm';
    GA_diff_NMCs = ft_freqgrandaverage(cfg,allsource_NMCs{s});
    GA_diff_HMCs = ft_freqgrandaverage(cfg,allsource_HMCs{s});
    %Using template grid;
    sourceTmpl = [];
    sourceTmpl.inside = grid.inside;
    sourceTmpl.dim = grid.dim;
    sourceTmpl.pos = grid.pos;
    sourceTmpl.unit = grid.unit;
    % Apply this source 'template' to sound and movie data and mask the data;
    GA_source_NMCs_temp = sourceTmpl;
    GA_source_NMCs_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMCs_temp.powspctrm(sourceTmpl.inside)= GA_diff_NMCs.powspctrm;
    GA_source_HMCs_temp = sourceTmpl;
    GA_source_HMCs_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMCs_temp.powspctrm(sourceTmpl.inside)= GA_diff_HMCs.powspctrm;
    %Interpolate the parameter 'powspctrm';
    cfg = [];
    cfg.downsample = 2;
    cfg.parameter = 'powspctrm';
    GA_source_NMCs = ft_sourceinterpolate(cfg, GA_source_NMCs_temp, mri); 
    GA_source_HMCs = ft_sourceinterpolate(cfg, GA_source_HMCs_temp, mri);
    
    %store 2-3 Hz power difference at the grids of interests  
    allsubj_grid_delta.NMCs(s,:) = GA_source_NMCs_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMCs(s,:) = GA_source_HMCs_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';

end

% save delta power in the significant cluster;
save allsubj_Left_Motor_delta_diff allsubj_grid_delta

% remove bad subjects;
bad_subj = [3 12 20];
allsubj_grid_delta.NMCs(bad_subj,:) = [];
allsubj_grid_delta.HMCs(bad_subj,:) = [];
allsubj_behavior_scores(bad_subj,:) = [];

% average delta power across grids of LeftMotor cluster;
mean_delta_pwr_NMCs = squeeze(nanmean(allsubj_grid_delta.NMCs,2));
mean_delta_pwr_HMCs = squeeze(nanmean(allsubj_grid_delta.HMCs,2));

methd = 'zscore';
valeurs = 'std';
correl = [];
correl(:,1) = round(normalize(cell2mat(allsubj_behavior_scores.cr_NMA) - cell2mat(allsubj_behavior_scores.cr_NMS),methd,valeurs),5);
correl(:,2) = round(normalize(cell2mat(allsubj_behavior_scores.cr_HMA) - cell2mat(allsubj_behavior_scores.cr_HMS),methd,valeurs),5);
correl(:,3) = round(normalize(mean_delta_pwr_NMCs,methd,valeurs),5);
correl(:,4) = round(normalize(mean_delta_pwr_HMCs,methd,valeurs),5);

% plot the correlation;
close all;
X1 = correl(:,1); Y1 = correl(:,3);
X2 = correl(:,2); Y2 = correl(:,4);

subplot(1,2,1)
s = scatter(X1,Y1);
s.Marker = 'o';
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceColor = [0.5 0.5 0.5];
s.SizeData = 15;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XAxis.Limits = [-3 3];
ax.YAxis.Limits = [-2 3];
ax.XTick = [-3 -2 -1 1 2 3];
ax.YTick = [-2 -1 1 2 3];
ax.YTickLabels = [-2 -1 1 2 3];
title('NMCs');
hold on
line = lsline;
line.LineWidth = 1.5;
line.Color = 'k';
subplot(1,2,2)
s = scatter(X2,Y2);
s.Marker = 'o';
s.MarkerEdgeColor = [0 0 0];
s.MarkerFaceColor = [0.5 0.5 0.5];
s.SizeData = 15;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.XAxis.Limits = [-3 3];
ax.YAxis.Limits = [-2 3];
ax.XTick = [-3 -2 -1 1 2 3];
ax.YTick = [-2 -1 1 2 3];
ax.YTickLabels = [-2 -1 1 2 3];
title('HMCs');
hold on
line = lsline;
line.LineWidth = 1.5;
line.Color = 'k';
set(gcf,'color','w');

% calculate coef + pvalue of the correlation (NMCs and HMCs);
[coef,pval] = corr(correl(:,1),correl(:,3),'type','Pearson','tail','right');

%% Now check the position of the selected ROI grids on the head model;
close all;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd xxx\headmodel\; addpath xxx\source_analyses\headmodel\; 
load standard_bem; 
load standard_mri; 
load('grid.mat');
grid = ft_convert_units(grid,'mm');

%your significant grids of the left motor area;
source_grids = cell2mat(NMCs_significant_grids.grid_id);
grid.pos = grid.pos(source_grids,:);

%plot only the grid positions within the brain and the BEM;
figure;
ft_plot_mesh(grid.pos)
ft_plot_mesh(vol.bnd(1), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
ft_plot_mesh(vol.bnd(2),'edgecolor','none','facealpha',0.4);

%%