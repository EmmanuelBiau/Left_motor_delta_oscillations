%% source grandaverage across participants;
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd C:\toolbox\fieldtrip\template\headmodel; load standard_mri;
cd xxx\source_analyses\headmodel; load('grid.mat');

% subjects = [1];
subjects = [1:2 4:26];

%Create a structure with all the MI_source data together;
clear allsubj_source
scount = 0;
for s = subjects

    %load participant's source data;
    cd(['xxx\subj_',num2str(s)]);
    load(['source_data_subj',num2str(s)]);
    
    %normalise power with pre-stim baseline (-0.7 to -0.2s);
    baseline_type = 'relative'; 
    bs_start = -0.7;
    bs_end = -0.2;
    
    cfg = [];
    cfg.parameter = 'powspctrm';
    cfg.baseline = [bs_start bs_end];
    cfg.baselinetype = baseline_type;
    source_NMAC = ft_freqbaseline(cfg,source_NMAC);
    source_NMSC = ft_freqbaseline(cfg,source_NMSC);
    source_HMAC = ft_freqbaseline(cfg,source_HMAC);
    source_HMSC = ft_freqbaseline(cfg,source_HMSC);
    
    %Select the time-window and frequency band of interest;
    cfg = [];
    cfg.latency = [3 9];
    cfg.avgovertime = 'yes';
    cfg.frequency = [2 3];
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

%Extract delta power at grids of interest;
load index_grid_Left_Motor;

clear allsubj_grid_delta
for s = subjects

    cfg = [];
    cfg.parameter = 'powspctrm';
    GA_NMAC = ft_freqgrandaverage(cfg,allsubj_source.NMAC{s});
    GA_NMSC = ft_freqgrandaverage(cfg,allsubj_source.NMSC{s});
    GA_HMAC = ft_freqgrandaverage(cfg,allsubj_source.HMAC{s});
    GA_HMSC = ft_freqgrandaverage(cfg,allsubj_source.HMSC{s});
    %Using template grid;
    sourceTmpl = [];
    sourceTmpl.inside = grid.inside;
    sourceTmpl.dim = grid.dim;
    sourceTmpl.pos = grid.pos;
    sourceTmpl.unit = grid.unit;
    % Apply this source 'template' to sound and movie data and mask the data;
    GA_source_NMAC_temp = sourceTmpl;
    GA_source_NMAC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMAC_temp.powspctrm(sourceTmpl.inside)= GA_NMAC.powspctrm;
    GA_source_NMSC_temp = sourceTmpl;
    GA_source_NMSC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMSC_temp.powspctrm(sourceTmpl.inside)= GA_NMSC.powspctrm;
    GA_source_HMAC_temp = sourceTmpl;
    GA_source_HMAC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMAC_temp.powspctrm(sourceTmpl.inside)= GA_HMAC.powspctrm;
    GA_source_HMSC_temp = sourceTmpl;
    GA_source_HMSC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMSC_temp.powspctrm(sourceTmpl.inside)= GA_HMSC.powspctrm;
    
    %store 2-3 Hz power difference at the grids of interests
    allsubj_grid_delta.NMAC(s,:) = GA_source_NMAC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.NMSC(s,:) = GA_source_NMSC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMAC(s,:) = GA_source_HMAC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMSC(s,:) = GA_source_HMSC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';

end

bad_subj = [3 12 20];
allsubj_grid_delta.NMAC(bad_subj,:) = [];
allsubj_grid_delta.NMSC(bad_subj,:) = [];
allsubj_grid_delta.HMAC(bad_subj,:) = [];
allsubj_grid_delta.HMSC(bad_subj,:) = [];

%average delta power across grid of ROI;
allsubj_delta_mean_power = [];
allsubj_delta_mean_power(:,1) = squeeze(nanmean(allsubj_grid_delta.NMAC,2));
allsubj_delta_mean_power(:,2) = squeeze(nanmean(allsubj_grid_delta.NMSC,2));
allsubj_delta_mean_power(:,3) = squeeze(nanmean(allsubj_grid_delta.HMAC,2));
allsubj_delta_mean_power(:,4) = squeeze(nanmean(allsubj_grid_delta.HMSC,2));

%save the mean delta power across conditions in the significant cluster;
save allsubj_Left_Motor_delta_mean_power allsubj_delta_mean_power

%%