%% Source reconstruction;
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%subjects ID;
subjects = [1 2 4:26];

for s = subjects

    %Load Template 128-Easycap headmodel;
    cd xxx\source_analyses\headmodel\;
    load('elec_ft.mat');
    load('grid.mat');
    load('vol.mat');
        
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
    filters = cell2mat(source_tmp.avg.filter(source_tmp.inside)); % keep filters
    
    %save the sound and movie filters;
    mkdir(['subj_',num2str(s)]); save(['subj_',num2str(s),'_filters'],'filters');
    
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
    
    %NMAC condition;
    cfg = [];
    cfg.output = 'pow';
    cfg.method = 'wavelet'; 
    cfg.foi = 1:1:5;
    cfg.toi = -1:0.01:10;
    cfg.width = 5;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    source_NMAC = ft_freqanalysis(cfg,data_source);
    
    if ~ismember(s,[18 26])
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
        source_NMAI = ft_freqanalysis(cfg,data_source);
    end
    
    %NMSC condition;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    source_NMSC = ft_freqanalysis(cfg,data_source);
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
    source_NMSI = ft_freqanalysis(cfg,data_source);
    
    %HMAC condition;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,3) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    source_HMAC = ft_freqanalysis(cfg,data_source);
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,3) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
    source_HMAI = ft_freqanalysis(cfg,data_source);
     
    %HMSC condition;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,4) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    source_HMSC = ft_freqanalysis(cfg,data_source);
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,4) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
    source_HMSI = ft_freqanalysis(cfg,data_source);
    
    %save Data_source for each participant;
    if ~ismember(s,[18 26])
        save(['source_data_subj',num2str(s)],'source_NMAC', 'source_NMSC', 'source_HMAC', 'source_HMSC','source_NMAI', 'source_NMSI', 'source_HMAI', 'source_HMSI');
    elseif ismember(s,[18 26])
        save(['source_data_subj',num2str(s)],'source_NMAC', 'source_NMSC', 'source_HMAC', 'source_HMSC','source_NMSI', 'source_HMAI', 'source_HMSI');
    end
    
    clearvars -except subjects

end

%%