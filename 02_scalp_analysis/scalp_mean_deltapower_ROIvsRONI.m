%% Compute the mean delta power for ANOVA comparisons;
clear;clc; addpath C:\toolbox\fieldtrip\; 
load allsubj_scalp_tfrs; load allsubj_behavior_scores

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

%% Compute Mean delta power across conditions in ROI and RONI pools;

%Chose the time-window of interest;
TWonset  = 3;  
TWoffset = 9;

%Chose the frequencies of interest delta [2-3Hz] or theta [4-8Hz] for control:
Frq_Min = 2;    
Frq_Max = 3; 

pool = 'roi';

%Chose the electrodes of interest:
switch pool  
    case 'roi'
        elec_interest_delta = {'F1','Fz','F2','FFC3h','FFC1h','FFC2h','FFC4h','FC3', 'FC1', 'FCz', 'FC2','FC4','FCC3h', 'FCC1h', 'FCC2h', 'fCC4h','C1','Cz','C2'};
    case 'roni'
        elec_interest_delta = {'POO1','Oz','Ol2h','POz','POO2','O2','POO10h','PO10','PPO2h','PO4','PPO6h','PO8','Pz'};
end

%Extract mean power for all electrodes in window/freq of interest;
clear Power_all_conditions_delta
Power_all_conditions_delta = {};
TWonset = find(round(avg_NMAC.time,2)== TWonset);
TWoffset = find(round(avg_NMAC.time,2)== TWoffset);
Power_all_conditions_delta.pwr_NMAC = squeeze(nanmean(nanmean(avg_NMAC.powspctrm(:,:,Frq_Min:Frq_Max,TWonset:TWoffset),4),3));
Power_all_conditions_delta.pwr_NMSC = squeeze(nanmean(nanmean(avg_NMSC.powspctrm(:,:,Frq_Min:Frq_Max,TWonset:TWoffset),4),3));
Power_all_conditions_delta.pwr_HMAC = squeeze(nanmean(nanmean(avg_HMAC.powspctrm(:,:,Frq_Min:Frq_Max,TWonset:TWoffset),4),3));
Power_all_conditions_delta.pwr_HSMC = squeeze(nanmean(nanmean(avg_HMSC.powspctrm(:,:,Frq_Min:Frq_Max,TWonset:TWoffset),4),3));
Power_all_conditions_delta.power_diff_NMCs = squeeze(nanmean(nanmean(diff_NMCs.powspctrm(:,:,Frq_Min:Frq_Max,TWonset:TWoffset),4),3));
Power_all_conditions_delta.power_diff_HMCs = squeeze(nanmean(nanmean(diff_HMCs.powspctrm(:,:,Frq_Min:Frq_Max,TWonset:TWoffset),4),3)); 
Power_all_conditions_delta.label = avg_NMAC.label';

clear mean_power_pool_delta
idx_elec_delta = ismember(Power_all_conditions_delta.label,elec_interest_delta);
mean_power_pool_delta = zeros(length(Power_all_conditions_delta.pwr_NMAC(:,1)),4);
mean_power_pool_delta(:,1) = mean(Power_all_conditions_delta.pwr_NMAC(:,idx_elec_delta),2);
mean_power_pool_delta(:,2) = mean(Power_all_conditions_delta.pwr_NMSC(:,idx_elec_delta),2);
mean_power_pool_delta(:,3) = mean(Power_all_conditions_delta.pwr_HMAC(:,idx_elec_delta),2);
mean_power_pool_delta(:,4) = mean(Power_all_conditions_delta.pwr_HSMC(:,idx_elec_delta),2);
mean_power_pool_delta(:,5) = mean(Power_all_conditions_delta.power_diff_NMCs(:,idx_elec_delta),2);
mean_power_pool_delta(:,6) = mean(Power_all_conditions_delta.power_diff_HMCs(:,idx_elec_delta),2);

%Create a table to store the delta power in the four conditions;
bad_subj = [3 12 20];
clear allsubj_delta_mean_power
allsubj_delta_mean_power = NaN(size(mean_power_pool_delta,1) + numel(bad_subj),size(mean_power_pool_delta,2)); 
idx = setdiff(1:size(mean_power_pool_delta,1) + numel(bad_subj), bad_subj);
allsubj_delta_mean_power(idx,:) = mean_power_pool_delta;
allsubj_delta_mean_power = [zeros(size(allsubj_delta_mean_power(:,1),1),1) allsubj_delta_mean_power];
allsubj_delta_mean_power = array2table(allsubj_delta_mean_power,'VariableNames',{'participant','pwr_NMA','pwr_NMS','pwr_HMA','pwr_HMS','pwr_Diff_NMCs','pwr_Diff_HMCs'}); 
allsubj_delta_mean_power.participant = allsubj_behavior_scores.participant;
clear idx

bad_subj = [3 12 20];
allsubj_delta_mean_power(bad_subj,:) = [];

if Frq_Min == 4 && Frq_Max == 8
    allsubj_theta_mean_power_ROI = allsubj_delta_mean_power;
    save allsubj_scalp_theta_mean_power allsubj_theta_mean_power_ROI
elseif Frq_Min == 2 && Frq_Max == 3
    switch pool  
        case 'roi'
            allsubj_delta_mean_power_ROI = allsubj_delta_mean_power;
        case 'roni'
            allsubj_delta_mean_power_RONI = allsubj_delta_mean_power;
    end
end

%save mean delta power at scalp level;
save allsubj_scalp_delta_mean_power allsubj_delta_mean_power_ROI allsubj_delta_mean_power_RONI;

%% PLOT the Delta mean power in the four conditions;

%Plot mean delta power from ROI pool in NMA, NMS, HMA, and HMS conditions;
dat = allsubj_theta_mean_power_ROI{:,[3 2 5 4]};           
col = [0 0 0];                              
dotsize = 15;                                
smooth = 0;                               

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

%configure axis;
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.XAxis.Label.Visible ='on';
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'NMS','NMA','HMS','HMA'};
ylabel('Normalized Power');
set(gcf,'color','w');

%Plot differences NMA-NMS and HMA-HMS in ROI and RONI;
dat = [allsubj_delta_mean_power_ROI{:,6} allsubj_delta_mean_power_RONI{:,6} allsubj_delta_mean_power_ROI{:,7} allsubj_delta_mean_power_RONI{:,7}];           
col = [0 0 0];                              
dotsize = 15;                                
smooth = 0;                               

%Generate your plot by calling the function;
close all;
vert_rain_plot(dat,col,dotsize,smooth)
pause(0.01);

%configure axis;
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.XAxis.Label.Visible ='on';
ax.XTick = [1 2 3 4];
ax.XTickLabel = {'ROI-NMCs','RONI-NMCs','ROI-HMCs','RONI-HMCs'};
ylabel('Normalized Power');
set(gcf,'color','w');

%% 