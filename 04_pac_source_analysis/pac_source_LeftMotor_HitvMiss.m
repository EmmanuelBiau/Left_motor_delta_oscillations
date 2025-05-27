%% Compute delta-beta PAC in the Left Motor cluster for each participant;

% Source reconstruction;
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load Template 128-Easycap headmodel;
addpath xxx\source_analyses\headmodel\;
load('elec_ft.mat');
load('grid.mat');
load('vol.mat');
load index_grid_Left_Motor;

%subjects ID;
subjects = [1 2 4:26];

for s = subjects

rng('shuffle');
       
    %Load preprocessed data;
    load(['data_clean_padded_subj_',num2str(s)]);
    
    %Prepare leadfield;
    cfg = [];
    cfg.grid = grid; 
    cfg.vol = vol;
    cfg.elec = elec_ft;
    grid = ft_prepare_leadfield(cfg,preproc_data);
    
    %Time-lock the data get common source filters;
    cfg = [];
    cfg.channel = {'all', '-AFz', '-TP10'};
    cfg.keeptrials = 'yes';
    cfg.covariance = 'yes';
    cfg.covariancewindow = 'all';
    timelock_data = ft_timelockanalysis(cfg, preproc_data);
    
    %Move the timelocked movie data into source space and create the filter;
    cfg = []; 
    cfg.method = 'lcmv';
    cfg.grid = grid;
    cfg.headmodel = vol;
    cfg.elec = elec_ft;
    cfg.channel = {'all', '-AFz', '-TP10'};
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.fixedori = 'yes';
    cfg.lcmv.lambda = '5%';
    source_tmp = ft_sourceanalysis(cfg, timelock_data);
    filters = cell2mat(source_tmp.avg.filter(source_tmp.inside)); % keep filters
    
    %Create a data structure for the source data;
    data_source = [];
    data_source.time = preproc_data.time;       
    data_source.trialinfo = preproc_data.trialinfo; 
    data_source.fsample = preproc_data.fsample;     
    
    %Create labels for each virtual electrode;
    for c = 1 : sum(grid.inside)
        label{c,1} = ['S' num2str(c)]; 
    end
    
    data_source.label = label;
    %For each trial, apply filters to the recorded data;
    bad_chan = {'AFz', 'TP10'};
    idx = ~ismember(preproc_data.label,bad_chan);
    
    for j = 1 : numel(preproc_data.trial) 
        data_source.trial{1,j} = single(filters*preproc_data.trial{1,j}(idx,:)); %remove AFz and TP10 rows = 127 tot labels;
    end
    
    clear source_tmp;

    %select significant left ROI from source analyses;
    cfg = [];
    cfg.latency = [3 9];
    cfg.channel = data_source.label(cell2mat(NMCs_significant_grids.grid_index));
    data_temp = ft_selectdata(cfg,data_source);
    data_temp = rmfield(data_temp,'time');
    clear data_source
    
    %get trial info;
    trlinfo = zeros(numel(data_temp.trialinfo),1);
    for j = 1 : numel(trlinfo)
        trlinfo(j,1) = data_temp.trialinfo{j}.condition;
        trlinfo(j,2) = data_temp.trialinfo{j}.hit;
    end
    
    min_samp = round(min([sum(trlinfo(:,1) == 1 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 1 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 0)])*0.8);
    
    if min_samp ==0  
     min_samp = round(min([sum(trlinfo(:,1) == 1 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 0)])*0.8);
    end  
    
    clear conc_trial_NMA info_nma_temp conc_trial_NMS info_nms_temp conc_trial_HMA info_hma_temp conc_trial_HMS info_hms_temp conc_trial_time
    
    for itr = 1:50
        
        [idx_nmac,~] = find(trlinfo(:,1) == 1 & trlinfo(:,2) == 1); 
        [idx_nmsc,~] = find(trlinfo(:,1) == 2 & trlinfo(:,2) == 1); 
        [idx_hmac,~] = find(trlinfo(:,1) == 3 & trlinfo(:,2) == 1); 
        [idx_hmsc,~] = find(trlinfo(:,1) == 4 & trlinfo(:,2) == 1); 
        
        %little trick for nami in subjects 18 and 26;
        if round(min([sum(trlinfo(:,1) == 1 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 1 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 0) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 0)])*0.8) == 0;
            [idx_nmai,~] = find(trlinfo(:,1) == 1 & trlinfo(:,2) == 1); 
        else 
            [idx_nmai,~] = find(trlinfo(:,1) == 1 & trlinfo(:,2) == 0); 
        end
            
        [idx_nmsi,~] = find(trlinfo(:,1) == 2 & trlinfo(:,2) == 0); 
        [idx_hmai,~] = find(trlinfo(:,1) == 3 & trlinfo(:,2) == 0); 
        [idx_hmsi,~] = find(trlinfo(:,1) == 4 & trlinfo(:,2) == 0);
        rand_trial_nmac = sortrows(idx_nmac(randperm(length(idx_nmac),min_samp)));
        rand_trial_nmsc = sortrows(idx_nmsc(randperm(length(idx_nmsc),min_samp)));
        rand_trial_hmac = sortrows(idx_hmac(randperm(length(idx_hmac),min_samp)));
        rand_trial_hmsc = sortrows(idx_hmsc(randperm(length(idx_hmsc),min_samp)));
        rand_trial_nmai = sortrows(idx_nmai(randperm(length(idx_nmai),min_samp)));
        rand_trial_nmsi = sortrows(idx_nmsi(randperm(length(idx_nmsi),min_samp)));
        rand_trial_hmai = sortrows(idx_hmai(randperm(length(idx_hmai),min_samp)));
        rand_trial_hmsi = sortrows(idx_hmsi(randperm(length(idx_hmsi),min_samp)));
        trial_nmac_temp = data_temp.trial(1,rand_trial_nmac);
        trial_nmsc_temp = data_temp.trial(1,rand_trial_nmsc);
        trial_hmac_temp = data_temp.trial(1,rand_trial_hmac);
        trial_hmsc_temp = data_temp.trial(1,rand_trial_hmsc);
        trial_nmai_temp = data_temp.trial(1,rand_trial_nmai);
        trial_nmsi_temp = data_temp.trial(1,rand_trial_nmsi);
        trial_hmai_temp = data_temp.trial(1,rand_trial_hmai);
        trial_hmsi_temp = data_temp.trial(1,rand_trial_hmsi);
        
         
        for ttp = 1:length(trial_nmac_temp)
        
            trial_nmac_temp{1,ttp} = trial_nmac_temp{1,ttp}(:,1:end-1); 
            trial_nmsc_temp{1,ttp} = trial_nmsc_temp{1,ttp}(:,1:end-1); 
            trial_hmac_temp{1,ttp} = trial_hmac_temp{1,ttp}(:,1:end-1); 
            trial_hmsc_temp{1,ttp} = trial_hmsc_temp{1,ttp}(:,1:end-1); 
            trial_nmai_temp{1,ttp} = trial_nmai_temp{1,ttp}(:,1:end-1); 
            trial_nmsi_temp{1,ttp} = trial_nmsi_temp{1,ttp}(:,1:end-1); 
            trial_hmai_temp{1,ttp} = trial_hmai_temp{1,ttp}(:,1:end-1); 
            trial_hmsi_temp{1,ttp} = trial_hmsi_temp{1,ttp}(:,1:end-1); 
        
        end    
            
        info_nmac_temp(itr,1) = data_temp.trialinfo(rand_trial_nmac(1),1);
        info_nmsc_temp(itr,1) = data_temp.trialinfo(rand_trial_nmsc(1),1); 
        info_hmac_temp(itr,1) = data_temp.trialinfo(rand_trial_hmac(1),1); 
        info_hmsc_temp(itr,1) = data_temp.trialinfo(rand_trial_hmsc(1),1); 
        info_nmai_temp(itr,1) = data_temp.trialinfo(rand_trial_nmai(1),1);
        info_nmsi_temp(itr,1) = data_temp.trialinfo(rand_trial_nmsi(1),1); 
        info_hmai_temp(itr,1) = data_temp.trialinfo(rand_trial_hmai(1),1); 
        info_hmsi_temp(itr,1) = data_temp.trialinfo(rand_trial_hmsi(1),1); 
        
        conc_trial_NMAC{1,itr} = cat(2,trial_nmac_temp{:});
        conc_trial_NMSC{1,itr} = cat(2,trial_nmsc_temp{:});
        conc_trial_HMAC{1,itr} = cat(2,trial_hmac_temp{:});
        conc_trial_HMSC{1,itr} = cat(2,trial_hmsc_temp{:});
        conc_trial_NMAI{1,itr} = cat(2,trial_nmai_temp{:});
        conc_trial_NMSI{1,itr} = cat(2,trial_nmsi_temp{:});
        conc_trial_HMAI{1,itr} = cat(2,trial_hmai_temp{:});
        conc_trial_HMSI{1,itr} = cat(2,trial_hmsi_temp{:});
        
        clear time_temp
        for conc = 1:min_samp  
            time_temp(conc,:) = ((max(cfg.latency) - min(cfg.latency))*(conc-1) + 0.0005):1/500:((max(cfg.latency) - min(cfg.latency))*conc); 
        end
        
        info_time_temp = reshape(time_temp',[],1)';
        time_nmac_temp{1,itr} = info_time_temp;
        time_nmsc_temp{1,itr} = info_time_temp;
        time_hmac_temp{1,itr} = info_time_temp;
        time_hmsc_temp{1,itr} = info_time_temp;
        time_nmai_temp{1,itr} = info_time_temp;
        time_nmsi_temp{1,itr} = info_time_temp;
        time_hmai_temp{1,itr} = info_time_temp;
        time_hmsi_temp{1,itr} = info_time_temp;

    end

    data_temp.trial = cat(2,conc_trial_NMAC,conc_trial_NMSC,conc_trial_HMAC,conc_trial_HMSC,conc_trial_NMAI,conc_trial_NMSI,conc_trial_HMAI,conc_trial_HMSI);
    data_temp.trialinfo = cat(1,info_nmac_temp,info_nmsc_temp,info_hmac_temp,info_hmsc_temp,info_nmai_temp,info_nmsi_temp,info_hmai_temp,info_hmsi_temp);
    data_temp.time = cat(2,time_nmac_temp,time_nmsc_temp,time_hmac_temp,time_hmsc_temp,time_nmai_temp,time_nmsi_temp,time_hmai_temp,time_hmsi_temp);

    %create the new structure with 50 trials in each condition, coming from the %concanated trials;
    data = data_temp;
    clear data_temp 
    
    %calculate MI;
    cfg = [];
    cfg.trials = 'all';
    cfg.phasefreq = [-0.5 0.5] + allsubj_deltabetapeak_LMC(s,2);
    cfg.powerfreq  = [-5 5] + allsubj_deltabetapeak_LMC(s,3); 
    cfg.nbins = 12;
    cfg.keeptrials = 'yes';
    pac = uni_getBasicPAC(cfg,data);
    clear data
    
    %reformat pac
    %pac = rmfield(pac,{'phase','bin_val','phapow','powpow'});
    pac.powspctrm = permute(pac.powspctrm,[4 1 2 3]);
    pac.dimord = 'rpt_chan_freq_time';
    
    cd(['xxx\HitvsMiss\subj_',num2str(s)]);
    save(['subj_',num2str(s),'_pac_source'],'pac');

end

% Compile delta-beta PAC data across subjects for ANOVA ANALYSIS;
for s = subjects
       
    cd (['xxx\HitvsMiss\subj_',num2str(s)]);
    load(['subj_',num2str(s),'_pac_source'],'pac');
    
    trlinfo = zeros(size(pac_beta.trialinfo));
    for trl = 1 : numel(pac_beta.trialinfo)
        trlinfo(trl,1) = pac_beta.trialinfo{trl,1}.condition;
        trlinfo(trl,2) = pac_beta.trialinfo{trl,1}.hit;
    end

    %sort trials x correctness x mask;
    pac_NMAC_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 1 & trlinfo(:,2)== 1,:,:,:),1),2));
    pac_NMSC_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 2 & trlinfo(:,2)== 1,:,:,:),1),2));
    pac_NMAI_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 1 & trlinfo(:,2)== 0,:,:,:),1),2));
    pac_NMSI_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 2 & trlinfo(:,2)== 0,:,:,:),1),2));  
    pac_HMAC_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 3 & trlinfo(:,2)== 1,:,:,:),1),2));
    pac_HMSC_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 4 & trlinfo(:,2)== 1,:,:,:),1),2));
    pac_HMAI_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 3 & trlinfo(:,2)== 0,:,:,:),1),2));
    pac_HMSI_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 4 & trlinfo(:,2)== 0,:,:,:),1),2));
    
    %compute it for each participant;
    allsubj_pac_beta_NMCs(s,1) = pac_NMAC_beta;
    allsubj_pac_beta_NMCs(s,2) = pac_NMAI_beta;
    allsubj_pac_beta_NMCs(s,3) = pac_NMSC_beta;
    allsubj_pac_beta_NMCs(s,4) = pac_NMSI_beta;
    allsubj_pac_beta_HMCs(s,1) = pac_HMAC_beta;
    allsubj_pac_beta_HMCs(s,2) = pac_HMAI_beta;
    allsubj_pac_beta_HMCs(s,3) = pac_HMSC_beta;
    allsubj_pac_beta_HMCs(s,4) = pac_HMSI_beta;
    pac_NMA_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 1,:,:,:),1),2));
    pac_NMS_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 2,:,:,:),1),2));
    pac_HMA_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 3,:,:,:),1),2));
    pac_HMS_beta = squeeze(nanmean(nanmean(pac_beta.powspctrm(trlinfo(:,1)== 4,:,:,:),1),2));
    allsubj_pac_beta_NMA_NMS(s,1) = pac_NMA_beta;
    allsubj_pac_beta_NMA_NMS(s,2) = pac_NMS_beta;
    allsubj_pac_beta_HMA_HMS(s,1) = pac_HMA_beta;
    allsubj_pac_beta_HMA_HMS(s,2) = pac_HMS_beta;

end

%Reshape data for ANOVA analysis;
%remove nmai data from subj 18 and 26 and replace with mean of the condition;
allsubj_pac_beta_NMCs([18 26],2) = NaN;
allsubj_pac_beta_HMCs([18 26],2) = NaN;
allsubj_pac_beta_NMCs([18 26],2) = nanmean(allsubj_pac_beta_NMCs(:,2));
allsubj_pac_beta_HMCs([18 26],2) = nanmean(allsubj_pac_beta_HMCs(:,2));

%remove bad subjects;
bad_subj = [3 12 20];
allsubj_pac_beta_NMCs(bad_subj,:) = [];
allsubj_pac_beta_HMCs(bad_subj,:) = [];
allsubj_pac_beta_NMA_NMS(bad_subj,:) = [];
allsubj_pac_beta_HMA_HMS(bad_subj,:) = [];

%store delta-beta PAC in NMCs and HMCs for anovas;
anovas_analyses_beta = round([allsubj_pac_beta_NMCs allsubj_pac_beta_HMCs],7);
anova_analyses_freq = [mean(anovas_analyses_beta(:,[1 3]),2) mean(anovas_analyses_beta(:,[5 7]),2)];
anovas_analyses_NMCsvsHMCs = round([allsubj_pac_beta_NMA_NMS allsubj_pac_beta_HMA_HMS],7);
anovas_analyses_beta = anovas_analyses_beta(:, [3 4 1 2 7 8 5 6]);

%save delta-beta pac data from the Left Motor cortex;
save allsubj_LeftMotor_deltabeta_pac_HvM anovas_analyses_beta;

%%