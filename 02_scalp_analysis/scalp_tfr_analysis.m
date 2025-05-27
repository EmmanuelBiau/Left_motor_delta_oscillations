%% Compute tfrs in each condition for all participants;
clear;clc; addpath C:\toolbox\fieldtrip\; 

subjects = [1 2 4:26];

cfg = [];
cfg.channels = 'all';
cfg.output = 'pow';
cfg.method = 'wavelet'; 
cfg.foi = 1:1:30;
cfg.toi = -5:0.04:15;
cfg.width = 5;

for s = subjects
 
    load(['\data_clean_padded_subj_',num2str(s)]);

    %IP NMAC;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_nmac = ft_freqanalysis(cfg,preproc_data);
    %dont compute nmai for 18 and 26 because they dont have misses in NMA condition;
    if ~ismember(s,[18 26]) 
        %IP NMAI;
        cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,1) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
        ip_nmai = ft_freqanalysis(cfg,preproc_data);
    end
    %IP NMSC;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_nmsc = ft_freqanalysis(cfg,preproc_data);
    %IP NMSI;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,2) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_nmsi = ft_freqanalysis(cfg,preproc_data);
    %IP HMAC;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,3) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_hmac = ft_freqanalysis(cfg,preproc_data);
    %IP HMAI;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,3) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_hmai = ft_freqanalysis(cfg,preproc_data);
    %IP HMSC;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,4) & isequal(x.hit,1), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_hmsc = ft_freqanalysis(cfg,preproc_data);
    %IP HMSI;
    cfg.trials = find(cell2mat(cellfun(@(x) isequal(x.condition,4) & isequal(x.hit,0), preproc_data.trialinfo, 'UniformOutput', false)));
    ip_hmsi = ft_freqanalysis(cfg,preproc_data);

    %save TFRs for each participant;
    if ~ismember(s,[18 26])  
        save(['TFRs_',num2str(s),'_perm'],'ip_nmac','ip_nmai','ip_nmsc','ip_nmsi','ip_hmac','ip_hmai','ip_hmsc','ip_hmsi');
    elseif ismember(s,[18 26]) 
        save(['TFRs_',num2str(s),'_perm'],'ip_nmac','ip_nmsc','ip_nmsi','ip_hmac','ip_hmai','ip_hmsc','ip_hmsi');
    end

end

% compile all subjects data;
all_ip_nmac = []; all_ip_nmai = [];
all_ip_nmsc = []; all_ip_nmsi = []; 
all_ip_hmac = []; all_ip_hmai = [];
all_ip_hmsc = []; all_ip_hmsi = [];

for s = subjects

    load(['TFRs_',num2str(s)]);   
    all_ip_nmac{1,s} = ip_nmac;
    all_ip_nmsc{1,s} = ip_nmsc;
    all_ip_nmsi{1,s} = ip_nmsi;
    all_ip_hmac{1,s} = ip_hmac;
    all_ip_hmai{1,s} = ip_hmai;
    all_ip_hmsc{1,s} = ip_hmsc;
    all_ip_hmsi{1,s} = ip_hmsi;
    
    if ~ismember(s,[18 26])
        all_ip_nmai{1,s} = ip_nmai;
    end
    
    keep all_ip_nmac all_ip_nmai all_ip_nmsc all_ip_nmsi all_ip_hmac all_ip_hmai all_ip_hmsc all_ip_hmsi

end

%save all subjects TFRS at scalp level;
% save allsubj_scalp_tfrs all_ip_nmac all_ip_nmai all_ip_nmsc all_ip_nmsi all_ip_hmac all_ip_hmai all_ip_hmsc all_ip_hmsi;

%remove bad subjects;
bad_subjects = [3 12 20];
all_ip_nmac(bad_subjects) = [];
all_ip_nmai(bad_subjects) = [];
all_ip_nmsc(bad_subjects) = [];
all_ip_nmsi(bad_subjects) = [];
all_ip_hmac(bad_subjects) = [];
all_ip_hmai(bad_subjects) = [];
all_ip_hmsc(bad_subjects) = [];
all_ip_hmsi(bad_subjects) = [];

%normalize individual power to the pre-stim baseline;
cfg=[];
cfg.parameter = 'powspctrm';
cfg.baseline = [-0.7 -0.2];
cfg.baselinetype = 'relative';
for s = 1:length(all_ip_nmac)
    all_ip_nmac{1,s} = ft_freqbaseline(cfg,all_ip_nmac{1,s});
    if ~ismember(s,[16 23])
    all_ip_nmai{1,s} = ft_freqbaseline(cfg,all_ip_nmai{1,s}); 
    end
    all_ip_nmsc{1,s} = ft_freqbaseline(cfg,all_ip_nmsc{1,s});
    all_ip_nmsi{1,s} = ft_freqbaseline(cfg,all_ip_nmsi{1,s});
    all_ip_hmac{1,s} = ft_freqbaseline(cfg,all_ip_hmac{1,s});    
    all_ip_hmai{1,s} = ft_freqbaseline(cfg,all_ip_hmai{1,s});
    all_ip_hmsc{1,s} = ft_freqbaseline(cfg,all_ip_hmsc{1,s});
    all_ip_hmsi{1,s} = ft_freqbaseline(cfg,all_ip_hmsi{1,s});   
end
    
all_ip_nmai = all_ip_nmai(~cellfun('isempty', all_ip_nmai));

%average power across subjects;
cfg = [];
cfg.keepindividual = 'yes'; 
avg_NMAC = ft_freqgrandaverage(cfg, all_ip_nmac{:});
avg_NMAI = ft_freqgrandaverage(cfg, all_ip_nmai{:});
avg_NMSC = ft_freqgrandaverage(cfg, all_ip_nmsc{:});
avg_NMSI = ft_freqgrandaverage(cfg, all_ip_nmsi{:});
avg_HMAC = ft_freqgrandaverage(cfg, all_ip_hmac{:});
avg_HMAI = ft_freqgrandaverage(cfg, all_ip_hmai{:});
avg_HMSC = ft_freqgrandaverage(cfg, all_ip_hmsc{:});
avg_HMSI = ft_freqgrandaverage(cfg, all_ip_hmsi{:});

diff_NMCs = avg_NMAC;
diff_NMCs.powspctrm = (avg_NMAC.powspctrm - avg_NMSC.powspctrm);
diff_HMCs = avg_HMAC;
diff_HMCs.powspctrm = (avg_HMAC.powspctrm - avg_HMSC.powspctrm);

% Cluster-based permutation tests;
cfg = [];
cfg_neighb.method = 'distance';
cfg.layout = 'easycapM15.mat';
neighbours = ft_prepare_neighbours(cfg_neighb, all_ip_nmac{1,1});

subj = size(all_ip_nmac,2);
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)= 1;
design(2,subj+1:2*subj) = 2;

cfg = [];
cfg.avgoverchan = 'no';
cfg.avgoverfreq = 'no';
cfg.avgovertime = 'no';
cfg.channel = {'all'};
cfg.latency = [-0.5 10];
cfg.frequency = [1 30];
cfg.design = design;
cfg.uvar = 1;
cfg.ivar = 2;
cfg.neighbours = neighbours;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 3;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.05;
cfg.numrandomization = 5000;
%spectrum [1-30Hz] difference in NMCs and HMs; 
stat_spectrum_diff_NMCs = ft_freqstatistics(cfg,all_ip_nmac{:},all_ip_nmsc{:}); 
stat_spectrum_diff_HMCs = ft_freqstatistics(cfg,all_ip_hmac{:},all_ip_hmsc{:}); 
%topography of delta [2-3Hz] difference in NMCs and HMs;
cfg.latency = [3 9];
cfg.frequency = [2 3];
cfg.avgoverfreq = 'yes';
cfg.avgovertime = 'yes';
stat_delta_diff_NMCs = ft_freqstatistics(cfg,all_ip_nmac{:},all_ip_nmsc{:}); 
stat_delta_diff_HMCs = ft_freqstatistics(cfg,all_ip_hmac{:},all_ip_hmsc{:}); 

%save CBP test outcomes; 
save allsubj_scalp_tfrs_stats stat_spectrum_diff_NMCs stat_spectrum_diff_HMCs stat_delta_diff_NMCs stat_delta_diff_HMCs;

%% Plot TFRs and topographies of NMCs/HMCs contrasts;

close all;
cfg=[];
%ROI channels;
cfg.channel = {'F1','Fz','F2','FFC3h','FFC1h','FFC2h','FFC4h','FC3', 'FC1', 'FCz', 'FC2','FC4','FCC3h', 'FCC1h', 'FCC2h', 'fCC4h','C1','Cz','C2'};
%RONI channels;
% cfg.channel = {'POO1','Oz','Ol2h','POz','POO2','O2','POO10h','PO10','PPO2h','PO4','PPO6h','PO8','Pz'};
cfg.layout = 'easycapM15';
cfg.parameter='powspctrm';	
cfg.xlim =[-0.5 10]; 
cfg.ylim =[1 30];
cfg.zlim = [-0.15 0.15];
cfg.colorbar = 'yes';
cfg.showlabels = 'yes';
figure;ft_singleplotTFR(cfg,diff_NMCs);
set(gcf,'color','w');
figure;ft_singleplotTFR(cfg,diff_HMCs);
set(gcf,'color','w');

cfg = [];
cfg.channel = {'all','-AFz','-TP10'}; 
cfg.xlim = [3 9];          
cfg.ylim = [2 3];                        
cfg.zlim = [-0.15 0.15];        
cfg.layout = 'easycapM15';
cfg.gridscale = 100;  
cfg.colorbar = 'yes'; 
cfg.comment = 'xlim'; 
cfg.marker = 'off'; 
cfg.markersymbol = '*';
cfg.highlight = 'on';
cfg.highlightchannel = stat_delta_diff_NMCs.label(stat_delta_diff_NMCs.mask==1);
figure; ft_topoplotTFR(cfg, diff_NMCs);
cfg.highlightchannel = stat_delta_diff_HMCs.label(stat_delta_diff_HMCs.mask==1);
figure; ft_topoplotTFR(cfg, diff_HMCs);

%%