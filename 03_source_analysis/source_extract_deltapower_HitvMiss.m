%% Extract delta power Hit versus Miss trials;
clearvars; clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd C:\toolbox\fieldtrip\template\headmodel; load standard_mri;
cd xxx\source_analyses\headmodel; load('grid.mat');

subjects = [1 2 4:26];

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
    NMSI_source = ft_selectdata(cfg,source_NMSI);
    HMAI_source = ft_selectdata(cfg,source_HMAI);
    HMSI_source = ft_selectdata(cfg,source_HMSI);
    
    if ~ismember(s,[18 26])
        NMAI_source = ft_selectdata(cfg,source_NMAI);
    end
    
    %Put all participants in a common structure;
    allsubj_source.NMAC{1,s} = NMAC_source;
    allsubj_source.NMSC{1,s} = NMSC_source;
    allsubj_source.HMAC{1,s} = HMAC_source;
    allsubj_source.HMSC{1,s} = HMSC_source;
    allsubj_source.NMSI{1,s} = NMSI_source;
    allsubj_source.HMAI{1,s} = HMAI_source;
    allsubj_source.HMSI{1,s} = HMSI_source;
    
    if ~ismember(s,[18 26])
        allsubj_source.NMAI{1,s} = NMAI_source;
    end

end

for s = subjects

    allsource_NMAC{s} =  allsubj_source.NMAC{s};  
    allsource_NMSC{s} =  allsubj_source.NMSC{s};     
    allsource_NMSI{s} =  allsubj_source.NMSI{s};     
    allsource_HMAC{s} =  allsubj_source.HMAC{s};
    allsource_HMAI{s} =  allsubj_source.HMAI{s};
    allsource_HMSC{s} =  allsubj_source.HMSC{s};     
    allsource_HMSI{s} =  allsubj_source.HMSI{s};     
    
    if ~ismember(s,[18 26])
        allsource_NMAI{s} =  allsubj_source.NMAI{s};
    end
    
    if ~ismember(s,[18 26])
        diff_NMAs_temp = allsubj_source.NMAC{s};
        diff_NMAs_temp.powspctrm = [];
        diff_NMAs_temp.powspctrm = (allsubj_source.NMAC{s}.powspctrm - allsubj_source.NMAI{s}.powspctrm)./allsubj_source.NMAI{s}.powspctrm;  
        allsource_NMA_HvM{s} = diff_NMAs_temp;
    end
    
    diff_NMSs_temp = allsubj_source.NMSC{s};
    diff_NMSs_temp.powspctrm = [];
    diff_NMSs_temp.powspctrm = (allsubj_source.NMSC{s}.powspctrm - allsubj_source.NMSI{s}.powspctrm)./allsubj_source.NMSI{s}.powspctrm;   
    allsource_NMS_HvM{s} = diff_NMSs_temp;
    
    diff_HMAs_temp = allsubj_source.HMAC{s};
    diff_HMAs_temp.powspctrm = [];
    diff_HMAs_temp.powspctrm = (allsubj_source.HMAC{s}.powspctrm - allsubj_source.HMAI{s}.powspctrm)./allsubj_source.HMAI{s}.powspctrm;  
    allsource_HMA_HvM{s} = diff_HMAs_temp;
    
    diff_HMSs_temp = allsubj_source.HMSC{s};
    diff_HMSs_temp.powspctrm = [];
    diff_HMSs_temp.powspctrm = (allsubj_source.HMSC{s}.powspctrm - allsubj_source.HMSI{s}.powspctrm)./allsubj_source.HMSI{s}.powspctrm;  
    allsource_HMS_HvM{s} = diff_HMSs_temp;
  
end

%Extract delta power at grids of interest;
load index_grid_Left_Motor;

clear allsubj_grid_delta

for s = subjects

    cfg = [];
    cfg.parameter = 'powspctrm';
    if ~ismember(s,[18 26])
        GA_NMA_HvM = ft_freqgrandaverage(cfg,allsource_NMA_HvM{s});
        GA_NMAI = ft_freqgrandaverage(cfg,allsubj_source.NMAI{s});
    end
    
    GA_NMAC =  ft_freqgrandaverage(cfg,allsubj_source.NMAC{s});  
    GA_NMSC =  ft_freqgrandaverage(cfg,allsubj_source.NMSC{s});     
    GA_NMSI =  ft_freqgrandaverage(cfg,allsubj_source.NMSI{s});     
    GA_HMAC =  ft_freqgrandaverage(cfg,allsubj_source.HMAC{s});
    GA_HMAI =  ft_freqgrandaverage(cfg,allsubj_source.HMAI{s});
    GA_HMSC =  ft_freqgrandaverage(cfg,allsubj_source.HMSC{s});    
    GA_HMSI =  ft_freqgrandaverage(cfg,allsubj_source.HMSI{s}); 
    GA_NMS_HvM = ft_freqgrandaverage(cfg,allsource_NMS_HvM{s});
    GA_HMA_HvM = ft_freqgrandaverage(cfg,allsource_HMA_HvM{s});
    GA_HMS_HvM = ft_freqgrandaverage(cfg,allsource_HMS_HvM{s});
    
    %Using template grid;
    sourceTmpl = [];
    sourceTmpl.inside = grid.inside;
    sourceTmpl.dim = grid.dim;
    sourceTmpl.pos = grid.pos;
    sourceTmpl.unit = grid.unit;
    % Apply this source 'template' to sound and movie data and mask the data;
    if ~ismember(s,[18 26])
        GA_source_NMA_HvM_temp = sourceTmpl;
        GA_source_NMA_HvM_temp.powspctrm = nan(size(grid.pos,1),1); 
        GA_source_NMA_HvM_temp.powspctrm(sourceTmpl.inside)= GA_NMA_HvM.powspctrm;        
        GA_source_NMAI_temp = sourceTmpl;
        GA_source_NMAI_temp.powspctrm = nan(size(grid.pos,1),1); 
        GA_source_NMAI_temp.powspctrm(sourceTmpl.inside)= GA_NMAI.powspctrm;
    end
    
    GA_source_NMAC_temp = sourceTmpl;
    GA_source_NMAC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMAC_temp.powspctrm(sourceTmpl.inside)= GA_NMAC.powspctrm;
    GA_source_NMSC_temp = sourceTmpl;
    GA_source_NMSC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMSC_temp.powspctrm(sourceTmpl.inside)= GA_NMSC.powspctrm;
    GA_source_NMSI_temp = sourceTmpl;
    GA_source_NMSI_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMSI_temp.powspctrm(sourceTmpl.inside)= GA_NMSI.powspctrm;
    GA_source_HMAC_temp = sourceTmpl;
    GA_source_HMAC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMAC_temp.powspctrm(sourceTmpl.inside)= GA_HMAC.powspctrm;
    GA_source_HMAI_temp = sourceTmpl;
    GA_source_HMAI_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMAI_temp.powspctrm(sourceTmpl.inside)= GA_HMAI.powspctrm;
    GA_source_HMSC_temp = sourceTmpl;
    GA_source_HMSC_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMSC_temp.powspctrm(sourceTmpl.inside)= GA_HMSC.powspctrm;
    GA_source_HMSI_temp = sourceTmpl;
    GA_source_HMSI_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMSI_temp.powspctrm(sourceTmpl.inside)= GA_HMSI.powspctrm;
    GA_source_NMS_HvM_temp = sourceTmpl;
    GA_source_NMS_HvM_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_NMS_HvM_temp.powspctrm(sourceTmpl.inside)= GA_NMS_HvM.powspctrm;
    GA_source_HMA_HvM_temp = sourceTmpl;
    GA_source_HMA_HvM_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMA_HvM_temp.powspctrm(sourceTmpl.inside)= GA_HMA_HvM.powspctrm;
    GA_source_HMS_HvM_temp = sourceTmpl;
    GA_source_HMS_HvM_temp.powspctrm = nan(size(grid.pos,1),1); 
    GA_source_HMS_HvM_temp.powspctrm(sourceTmpl.inside)= GA_HMS_HvM.powspctrm;
    
    %store 2-3 Hz power difference at the grids of interests 
    if ~ismember(s,[18 26])
        allsubj_grid_delta.NMA_HvM(s,:) = GA_source_NMA_HvM_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
        allsubj_grid_delta.NMAI(s,:) = GA_source_NMAI_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    end
    
    allsubj_grid_delta.NMAC(s,:) = GA_source_NMAC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.NMSC(s,:) = GA_source_NMSC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.NMSI(s,:) = GA_source_NMSI_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMAC(s,:) = GA_source_HMAC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMAI(s,:) = GA_source_HMAI_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMSC(s,:) = GA_source_HMSC_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMSI(s,:) = GA_source_HMSI_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.NMS_HvM(s,:) = GA_source_NMS_HvM_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMA_HvM(s,:) = GA_source_HMA_HvM_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';
    allsubj_grid_delta.HMS_HvM(s,:) = GA_source_HMS_HvM_temp.powspctrm(cell2mat(NMCs_significant_grids.grid_id))';

end

%average delta power across grid of ROI;
mean_delta_pwr = [];
max_mean = max([size(squeeze(nanmean(allsubj_grid_delta.NMA_HvM,2)),1) size(squeeze(nanmean(allsubj_grid_delta.NMS_HvM,2)),1) size(squeeze(nanmean(allsubj_grid_delta.HMA_HvM,2)),1) size(squeeze(nanmean(allsubj_grid_delta.HMS_HvM,2)),1)]);
mean_delta_pwr(:,1) = [squeeze(nanmean(allsubj_grid_delta.NMA_HvM,2)); zeros(max_mean-size(squeeze(nanmean(allsubj_grid_delta.NMA_HvM,2)),1))];
mean_delta_pwr(:,2) = [squeeze(nanmean(allsubj_grid_delta.NMS_HvM,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.NMS_HvM,2)),1))];
mean_delta_pwr(:,3) = [squeeze(nanmean(allsubj_grid_delta.HMA_HvM,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.HMA_HvM,2)),1))];
mean_delta_pwr(:,4) = [squeeze(nanmean(allsubj_grid_delta.HMS_HvM,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.HMS_HvM,2)),1))];
mean_delta_pwr(:,5) = [squeeze(nanmean(allsubj_grid_delta.NMAC,2)); zeros(max_mean-size(squeeze(nanmean(allsubj_grid_delta.NMAC,2)),1))];
mean_delta_pwr(:,6) = [squeeze(nanmean(allsubj_grid_delta.NMAI,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.NMAI,2)),1))];
mean_delta_pwr(:,7) = [squeeze(nanmean(allsubj_grid_delta.NMSC,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.NMSC,2)),1))];
mean_delta_pwr(:,8) = [squeeze(nanmean(allsubj_grid_delta.NMSI,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.NMSI,2)),1))];
mean_delta_pwr(:,9) = [squeeze(nanmean(allsubj_grid_delta.HMAC,2)); zeros(max_mean-size(squeeze(nanmean(allsubj_grid_delta.HMAC,2)),1))];
mean_delta_pwr(:,10) = [squeeze(nanmean(allsubj_grid_delta.HMAI,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.HMAI,2)),1))];
mean_delta_pwr(:,11) = [squeeze(nanmean(allsubj_grid_delta.HMSC,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.HMSC,2)),1))];
mean_delta_pwr(:,12) = [squeeze(nanmean(allsubj_grid_delta.HMSI,2)); zeros(max_mean - size(squeeze(nanmean(allsubj_grid_delta.HMSI,2)),1))];
mean_delta_pwr(18,:) = mean(mean_delta_pwr([1 2 4:11 13:17 19 21:25],:));
mean_delta_pwr(26,:) = mean_delta_pwr(18,:);

bad_subjects = [3 12 20];
mean_delta_pwr(bad_subject,:) = [];
allsubj_grid_delta_HvM = mean_delta_pwr; clear mean_delta_pwr;

%save the mean difference of delta power Hit versus Miss in the significant cluster;
save allsubj_Left_Motor_delta_HvM allsubj_grid_delta_HvM

%%
