% Reconstruct source EEG signal to detect the frequency with peak of delta/beta power in the Left Motor cortex; 
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd C:\toolbox\fieldtrip\template\headmodel; load standard_mri;
cd xxx\source_analyses\headmodel; load('grid.mat');

%subjects ID;
subjects = 1:26;

for s = subjects

    %Load Template 128-Easycap headmodel;
    cd xxx\source_analyses\headmodel\;
    load('elec_ft.mat');
    load('grid.mat');
    load('vol.mat');
        
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
    for c = 1 : sum(grid.inside); label{c,1} = ['S' num2str(c)]; end
    data_source.label = label;
    
    %For each trial, apply filters to the recorded data;
    bad_chan = {'AFz', 'TP10'};
    idx = ~ismember(preproc_data.label,bad_chan);
    
    for j = 1 : numel(preproc_data.trial) 
        data_source.trial{1,j} = single(filters*preproc_data.trial{1,j}(idx,:)); %remove AFz and TP10 rows = 127 tot labels;
    end
    
    clear source_tmp; 
    
    % save Data_source for each participant;
    save(['source_data_subj',num2str(s),'_lm'],'data_source');

end

% Delta/Beta peak detection of each participant/trial;
for s = subjects

    for ii = 1 : length(data_source.trial)
        data_tmp = nan(size(grid.pos,1),size(data_source.trial{ii},2));
        data_tmp(grid.inside,:) = data_source.trial{ii};
        data_source.trial{ii} = data_tmp;
    end
    
    for c = 1 : size(grid.pos,1); label{c,1} = ['S' num2str(c)]; end
    data_source.label = label;
    
    cfg = [];
    cfg.latency = [3 9]; 
    cfg.channel = label((cell2mat(NMCs_significant_grids.grid_id)));
    data_temp_LeftMotor = ft_selectdata(cfg,data_source);
        
    peak_power(s,1) = s;    
    
    %FFT_DELTA/BETA;
    cfg = [];
    cfg.keeptrials = 'yes';
    cfg.output = 'pow';
    cfg.method = 'wavelet';
    cfg.pad = 'nextpow2';
    cfg.foi = 2:0.1:4; 
    cfg.toi = 3:0.05:9; 
    cfg.width = 5;
    delta_fft = ft_freqanalysis(cfg,data_temp_LeftMotor);
    cfg.foi = 13:0.1:30; 
    beta_fft = ft_freqanalysis(cfg,data_temp_LeftMotor);

    %---DELTA and BETA Peaks (1-3Hz; 13-30Hz) -----------------------------------------------------------%
    addpath D:\Maastricht\maastricht_AV_speech_2020\PAC_analyses\analysis_ben\
    
    %create fake freq with a FT structure;
    freq{1} = delta_fft; 
    freq{2} = beta_fft;
    
    %create freq structure;
    for trl = 1:size(freq{1}.trialinfo,1)
        
        pow = squeeze(nanmean(nanmean(nanmean(delta_fft.powspctrm(trl,:,:,:),4),2),1));
        delta = struct('time',1,...
                   'freq',freq{1}.freq(freq{1}.freq>=1&freq{1}.freq<=30),...
                   'label',{{'dummy'}},...
                   'dimord','chan_freq_time',...
                   'powspctrm',pow(freq{1}.freq>=1&freq{1}.freq<=30)');
    
        %fractal property of power spectrum;  
        delta = uni_subtract1f(delta);
        delta_peak.delta_peak(trl,1) = uni_getPeak(delta.powspctrm',delta.freq',[2 4],0);
        delta_peak.trialinfo = freq{1}.trialinfo;
    
        pow = squeeze(nanmean(nanmean(nanmean(beta_fft.powspctrm(trl,:,:,:),4),2),1));
        beta = struct('time',1,...
                   'freq',freq{1}.freq(freq{1}.freq>=1&freq{1}.freq<=30),...
                   'label',{{'dummy'}},...
                   'dimord','chan_freq_time',...
                   'powspctrm',pow(freq{1}.freq>=1&freq{1}.freq<=30)');
    
        %fractal property of power spectrum;  
        beta = uni_subtract1f(beta);
        beta_peak.beta_peak(trl,1) = uni_getPeak(beta.powspctrm',beta.freq',[13 30],0);
        beta_peak.trialinfo = freq{1}.trialinfo;

    end
    
    cd(['xxx\subj_',num2str(s)]);
    save(['delta_peak_subj_' num2str(s)], 'delta_peak','beta_peak');

    peak_power(s,2) = mean(delta_peak.delta_peak);
    peak_power(s,3) = mean(delta_peak.beta_peak);
    
    clearvars -except s

end

% save peak information for all subjects;
save([xxx, 'allsubj_deltabetapeak_LMC'],'peak_power');

%% Matching delta peak of participant for each stimulus;
clearvars;clc; close all; 
addpath C:\toolbox\fieldtrip\; ft_defaults;
load list_stim;

%subjects ID;
subjects = [1 2 4:26];

Peak_delta_stim = [];
Peak_delta_stim.NMAC(:,1) = list_stim;
Peak_delta_stim.NMSC(:,1) = list_stim;
Peak_delta_stim.HMAC(:,1) = list_stim;
Peak_delta_stim.HMSC(:,1) = list_stim;
Peak_delta_stim.NMAI(:,1) = list_stim;
Peak_delta_stim.NMSI(:,1) = list_stim;
Peak_delta_stim.HMAI(:,1) = list_stim;
Peak_delta_stim.HMSI(:,1) = list_stim;

for s = subjects

    %load the left_motor data_source;
    cd(['xxx\subj_',num2str(s)]);
    load(['delta_peak_subj_',num2str(s)]); 
    tmp_peaks = [];
    
    for trl = 1:size(delta_peak.delta_peak,1)
        tmp_peaks(trl,1) = delta_peak.trialinfo{trl}.video_num;
        tmp_peaks(trl,2) = delta_peak.trialinfo{trl}.condition;
        tmp_peaks(trl,3) = delta_peak.trialinfo{trl}.hit;
        tmp_peaks(trl,4) = delta_peak.delta_peak(trl,1); 
    end
    
    tmp_peaks = sortrows(tmp_peaks,[1 2], 'ascend');

    for i = 1: length(Peak_delta_stim.NMAC)
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.NMAC(i,1) & tmp_peaks(:,2)== 1 & tmp_peaks(:,3)== 1) == 1 
            Peak_delta_stim.NMAC(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.NMAC(i,1) & tmp_peaks(:,2)== 1 & tmp_peaks(:,3)== 1,4);
        else
            Peak_delta_stim.NMAC(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.NMAC(i,1) & tmp_peaks(:,2)== 1 & tmp_peaks(:,3)== 0) == 1 
            Peak_delta_stim.NMAI(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.NMAI(i,1) & tmp_peaks(:,2)== 1 & tmp_peaks(:,3)== 0,4);
        else
            Peak_delta_stim.NMAI(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.NMSC(i,1) & tmp_peaks(:,2)== 2 & tmp_peaks(:,3)== 1) == 1 
            Peak_delta_stim.NMSC(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.NMSC(i,1) & tmp_peaks(:,2)== 2 & tmp_peaks(:,3)== 1,4);
        else
            Peak_delta_stim.NMSC(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.NMSC(i,1) & tmp_peaks(:,2)== 2 & tmp_peaks(:,3)== 0) == 1 
            Peak_delta_stim.NMSI(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.NMSI(i,1) & tmp_peaks(:,2)== 2 & tmp_peaks(:,3)== 0,4);
        else
            Peak_delta_stim.NMSI(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.HMAC(i,1) & tmp_peaks(:,2)== 3 & tmp_peaks(:,3)== 1) == 1
            Peak_delta_stim.HMAC(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.HMAC(i,1) & tmp_peaks(:,2)== 3 & tmp_peaks(:,3)== 1,4);
        else
            Peak_delta_stim.HMAC(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.HMAC(i,1) & tmp_peaks(:,2)== 3 & tmp_peaks(:,3)== 0) == 1 
            Peak_delta_stim.HMAI(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.HMAI(i,1) & tmp_peaks(:,2)==3 & tmp_peaks(:,3)== 0,4);
        else
            Peak_delta_stim.HMAI(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.HMSC(i,1) & tmp_peaks(:,2)== 4 & tmp_peaks(:,3)== 1) == 1 
            Peak_delta_stim.HMSC(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.HMSC(i,1) & tmp_peaks(:,2)== 4 & tmp_peaks(:,3)== 1,4);
        else
            Peak_delta_stim.HMSC(i,s+1) = NaN; 
        end
    
        if any(tmp_peaks(:,1)== Peak_delta_stim.HMSC(i,1) & tmp_peaks(:,2)== 4 & tmp_peaks(:,3)== 0) == 1 
            Peak_delta_stim.HMSI(i,s+1) = tmp_peaks(tmp_peaks(:,1)== Peak_delta_stim.HMSI(i,1) & tmp_peaks(:,2)== 4 & tmp_peaks(:,3)== 0,4);
        else
            Peak_delta_stim.HMSI(i,s+1) = NaN; 
        end   
        
    end

end

% save delta peak of participants for each stimulus; 
save Peak_delta_stim Peak_delta_stim

%% Compile Delta peak in the stimuli;
clearvars; clc;
load Peak_delta_stim

%Create a table with information;
stim_temp = readtable('Power_Stim.xlsx');
tmp_vid = []; tmp_cond = [];
for i = 1:size(stim_temp,1)  
    
    tmp_vid(i,1) = str2double(stim_temp.Video{i}(5:end-24));

    if strcmp(stim_temp.Video{i}(1:3), 'NMS') == 1
        tmp_cond(i,1) = 2;
    elseif strcmp(stim_temp.Video{i}(1:3), 'HMS') == 1
        tmp_cond(i,1) = 4;
    end

end

stim_temp.Video = [];
stim_temp = addvars(stim_temp,tmp_cond,tmp_vid,'Before','Head','NewVariableNames',{'Condition','Video'});
stim_temp = movevars(stim_temp,'Body','Before','Head'); stim_temp = movevars(stim_temp,'Full','Before','Head');
stim_temp = sortrows(stim_temp,[1 2],'ascend');
clear i tmp_cond tmp_vid
Peak_Freq_Stim = stim_temp;
clear stim_temp

% save delta peak in the stimuli data;
save Peak_freq_stim Peak_freq_stim 

%%