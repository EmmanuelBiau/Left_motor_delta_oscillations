%% Source reconstruction of each subject in the whole brain;
clearvars;clc; close all; 
addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load Template 128-Easycap headmodel;
addpath xxx\source_analysis\headmodel\;
load('elec_ft.mat');
load('grid.mat');
load('vol.mat');
load index_grid_Left_Motor; load Power_peak_info2;

%subjects ID;
subjects = [1 2 4:26];

for s = subjects

    %radnomize rand at each iteration/subj;
    rng('shuffle');    
  
    %Load preprocessed data;
    cd xxx\data_clean_padded\;
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
    filters = cell2mat(source_tmp.avg.filter(source_tmp.inside));

    %Create a data structure for the source data;
    data_source = [];
    data_source.time = preproc_data.time;       
    data_source.trialinfo = preproc_data.trialinfo; 
    data_source.fsample = preproc_data.fsample;     

    %Create labels for each virtual electrode;
    for c = 1 : sum(grid.inside); label{c,1} = ['S' num2str(c)]; end

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
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1) & isequal(x.hit,1), data_source.trialinfo, 'UniformOutput', false)));
    data_temp = ft_selectdata(cfg,data_source);
    data_temp = rmfield(data_temp,'time');
    clear data_source

    %get trial info;
    trlinfo = zeros(numel(data_temp.trialinfo),1);
    for j = 1 : numel(trlinfo)
        trlinfo(j,1) = data_temp.trialinfo{j}.condition;
        trlinfo(j,2) = data_temp.trialinfo{j}.hit;
    end

    min_samp = round(min([sum(trlinfo(:,1) == 1 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 2 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 3 & trlinfo(:,2) == 1) sum(trlinfo(:,1) == 4 & trlinfo(:,2) == 1)])*0.8);
    
    clear conc_trial_NMA info_nma_temp conc_trial_NMS info_nms_temp conc_trial_HMA info_hma_temp conc_trial_HMS info_hms_temp conc_trial_time
    
    for itr = 1:50
        
        [idx_nma,~] = find(trlinfo(:,1) == 1 & trlinfo(:,2) == 1); 
        [idx_nms,~] = find(trlinfo(:,1) == 2 & trlinfo(:,2) == 1); 
        [idx_hma,~] = find(trlinfo(:,1) == 3 & trlinfo(:,2) == 1); 
        [idx_hms,~] = find(trlinfo(:,1) == 4 & trlinfo(:,2) == 1); 
        rand_trial_nma = sortrows(idx_nma(randperm(length(idx_nma),min_samp)));
        rand_trial_nms = sortrows(idx_nms(randperm(length(idx_nms),min_samp)));
        rand_trial_hma = sortrows(idx_hma(randperm(length(idx_hma),min_samp)));
        rand_trial_hms = sortrows(idx_hms(randperm(length(idx_hms),min_samp)));
        
        trial_nma_temp = data_temp.trial(1,rand_trial_nma);
        trial_nms_temp = data_temp.trial(1,rand_trial_nms);
        trial_hma_temp = data_temp.trial(1,rand_trial_hma);
        trial_hms_temp = data_temp.trial(1,rand_trial_hms);
        
        for ttp = 1:length(trial_nma_temp)
        
            trial_nma_temp{1,ttp} = trial_nma_temp{1,ttp}(:,1:end-1); 
            trial_nms_temp{1,ttp} = trial_nms_temp{1,ttp}(:,1:end-1); 
            trial_hma_temp{1,ttp} = trial_hma_temp{1,ttp}(:,1:end-1); 
            trial_hms_temp{1,ttp} = trial_hms_temp{1,ttp}(:,1:end-1); 
            
        end    
            
        info_nma_temp(itr,1) = data_temp.trialinfo(rand_trial_nma(1),1);
        info_nms_temp(itr,1) = data_temp.trialinfo(rand_trial_nms(1),1); 
        info_hma_temp(itr,1) = data_temp.trialinfo(rand_trial_hma(1),1); 
        info_hms_temp(itr,1) = data_temp.trialinfo(rand_trial_hms(1),1); 
        
        conc_trial_NMA{1,itr} = cat(2,trial_nma_temp{:});
        conc_trial_NMS{1,itr} = cat(2,trial_nms_temp{:});
        conc_trial_HMA{1,itr} = cat(2,trial_hma_temp{:});
        conc_trial_HMS{1,itr} = cat(2,trial_hms_temp{:});
    
        clear time_temp
        for conc = 1:min_samp
            time_temp(conc,:) = ((max(cfg.latency) - min(cfg.latency))*(conc-1) + 0.0005):1/500:((max(cfg.latency) - min(cfg.latency))*conc);
        end
        
        info_time_temp = reshape(time_temp',[],1)';
        time_nma_temp{1,itr} = info_time_temp;
        time_nms_temp{1,itr} = info_time_temp;
        time_hma_temp{1,itr} = info_time_temp;
        time_hms_temp{1,itr} = info_time_temp;
    
    end

    data_temp.trial = cat(2,conc_trial_NMA,conc_trial_NMS,conc_trial_HMA,conc_trial_HMS);
    data_temp.trialinfo = cat(1,info_nma_temp,info_nms_temp,info_hma_temp,info_hms_temp);
    data_temp.time = cat(2,time_nma_temp,time_nms_temp,time_hma_temp,time_hms_temp);
    
    %create the new structure with 50 trials in each condition, coming from the %concanated trials;
    data = data_temp;
    clear data_temp 

    %calculate MI BETA;
    cfg = [];
    cfg.trials = 'all';
    cfg.phasefreq = [-0.5 0.5] + allsubj_deltabetapeak_LMC(s,2);
    cfg.powerfreq  = [-5 5] + allsubj_deltabetapeak_LMC(s,3); 
    cfg.nbins = 12;
    cfg.keeptrials = 'yes';
    pac_beta = uni_getBasicPAC(cfg,data);
    %reformat pac
    %pac = rmfield(pac,{'phase','bin_val','phapow','powpow'});
    pac_beta.powspctrm = permute(pac_beta.powspctrm,[4 1 2 3]);
    pac_beta.dimord = 'rpt_chan_freq_time';

    % save whole brain delta-beta pac data of the participant;
    save(['subj_',num2str(s),'_pac_source_wb'],'pac_wb_beta'); clear data

end

%%