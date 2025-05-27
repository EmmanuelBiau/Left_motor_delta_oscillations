%% Perform permutations on delta power Hit versus Miss trials;
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;

%Load mri template from FT and the grid template from one participant (e.g. grid of headmodel);
cd C:\toolbox\fieldtrip\template\headmodel; load standard_mri;
cd xxx\source_analyses\headmodel; load('grid.mat');

%Extract delta power at grids of interest;
load index_grid_Left_Motor;

subjects = [1 2 4:26];

iterations = 5000;

allsubj_grid_delta_perm.NMA_HvM = zeros(1,iterations);
allsubj_grid_delta_perm.NMS_HvM = zeros(1,iterations);
allsubj_grid_delta_perm.HMA_HvM = zeros(1,iterations);
allsubj_grid_delta_perm.HMS_HvM = zeros(1,iterations);

f = waitbar(0,'Running permutations now','Name','Work in Progress');

for s = subjects
   
    %load participant's source data;
    cd(['xxx\individual_source_data\subj_',num2str(s)]);
    load(['source_data_subj',num2str(s),'_ori']);
    
    NMAC_source.powspctrm = normalize(NMAC_source.powspctrm,'zscore','std');
    NMSC_source.powspctrm = normalize(NMSC_source.powspctrm,'zscore','std');
    HMAC_source.powspctrm = normalize(HMAC_source.powspctrm,'zscore','std');
    HMSC_source.powspctrm = normalize(HMSC_source.powspctrm,'zscore','std');
    NMSI_source.powspctrm = normalize(NMSI_source.powspctrm,'zscore','std');
    HMAI_source.powspctrm = normalize(HMAI_source.powspctrm,'zscore','std');
    HMSI_source.powspctrm = normalize(HMSI_source.powspctrm,'zscore','std');
    
    if ~ismember(s,[18 26])
        NMAI_source.powspctrm = normalize(NMAI_source.powspctrm,'zscore','std');
    end
    
    %Create template grid;
    sourceTmpl = [];
    sourceTmpl.inside = grid.inside;
    sourceTmpl.dim = grid.dim;
    sourceTmpl.pos = grid.pos;
    sourceTmpl.unit = grid.unit;
    
    %Apply this source 'template' to sound and movie data and mask the data;
    NMAC_source2 = sourceTmpl;
    NMAC_source2.trialinfo = NMAC_source.trialinfo;
    NMAC_source2.powspctrm = nan(size(NMAC_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(NMAC_source.powspctrm,1)
        NMAC_source2.powspctrm(i,sourceTmpl.inside)= NMAC_source.powspctrm(i,:);
    end
    NMAC_source2.powspctrm = NMAC_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    
    if ~ismember(s,[18 26])
        NMAI_source2 = sourceTmpl;
        NMAI_source2.trialinfo = NMAI_source.trialinfo;
        NMAI_source2.powspctrm = nan(size(NMAI_source.powspctrm,1),size(grid.pos,1)); 
        for i =1:size(NMAI_source.powspctrm,1)
            NMAI_source2.powspctrm(i,sourceTmpl.inside)= NMAI_source.powspctrm(i,:);
        end
    NMAI_source2.powspctrm = NMAI_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    end
    
    NMSC_source2 = sourceTmpl;
    NMSC_source2.trialinfo = NMSC_source.trialinfo;
    NMSC_source2.powspctrm = nan(size(NMSC_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(NMSC_source.powspctrm,1)
        NMSC_source2.powspctrm(i,sourceTmpl.inside)= NMSC_source.powspctrm(i,:);
    end
    NMSI_source2 = sourceTmpl;
    NMSI_source2.trialinfo = NMSI_source.trialinfo;
    NMSI_source2.powspctrm = nan(size(NMSI_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(NMSI_source.powspctrm,1)
        NMSI_source2.powspctrm(i,sourceTmpl.inside)= NMSI_source.powspctrm(i,:);
    end
    NMSI_source2.powspctrm = NMSI_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    NMSC_source2.powspctrm = NMSC_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    
    HMAC_source2 = sourceTmpl;
    HMAC_source2.trialinfo = HMAC_source.trialinfo;
    HMAC_source2.powspctrm = nan(size(HMAC_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(HMAC_source.powspctrm,1)
        HMAC_source2.powspctrm(i,sourceTmpl.inside)= HMAC_source.powspctrm(i,:);
    end
    HMAI_source2 = sourceTmpl;
    HMAI_source2.trialinfo = HMAI_source.trialinfo;
    HMAI_source2.powspctrm = nan(size(HMAI_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(HMAI_source.powspctrm,1)
        HMAI_source2.powspctrm(i,sourceTmpl.inside)= HMAI_source.powspctrm(i,:);
    end
    HMAI_source2.powspctrm = HMAI_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    HMAC_source2.powspctrm = HMAC_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    
    HMSC_source2 = sourceTmpl;
    HMSC_source2.trialinfo = HMSC_source.trialinfo;
    HMSC_source2.powspctrm = nan(size(HMSC_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(HMSC_source.powspctrm,1)
        HMSC_source2.powspctrm(i,sourceTmpl.inside)= HMSC_source.powspctrm(i,:);
    end
    HMSI_source2 = sourceTmpl;
    HMSI_source2.trialinfo = HMSI_source.trialinfo;
    HMSI_source2.powspctrm = nan(size(HMSI_source.powspctrm,1),size(grid.pos,1)); 
    for i =1:size(HMSI_source.powspctrm,1)
        HMSI_source2.powspctrm(i,sourceTmpl.inside)= HMSI_source.powspctrm(i,:);
    end
    HMSI_source2.powspctrm = HMSI_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    HMSC_source2.powspctrm = HMSC_source2.powspctrm(:,cell2mat(NMCs_significant_grids.grid_id));
    
    %shuffle hit/miss labels;
    NMAs_temp = cat(1, NMAC_source2.powspctrm, NMAI_source2.powspctrm);
    NMAC_source2.powspctrm = datasample(NMAs_temp,size(NMAC_source2.trialinfo,1),'Replace',false);
    idx = ~ismember(NMAs_temp(:,1),NMAC_source2.powspctrm(:,1));
    NMAI_source2.powspctrm = NMAs_temp(idx,:);
    NMSs_temp = cat(1, NMSC_source2.powspctrm, NMSI_source2.powspctrm);
    NMSC_source2.powspctrm = datasample(NMSs_temp,size(NMSC_source2.trialinfo,1),'Replace',false);
    idx = ~ismember(NMSs_temp(:,1),NMSC_source2.powspctrm(:,1));
    NMSI_source2.powspctrm = NMSs_temp(idx,:);
    HMAs_temp = cat(1, HMAC_source2.powspctrm, HMAI_source2.powspctrm);
    HMAC_source2.powspctrm = datasample(HMAs_temp,size(HMAC_source2.trialinfo,1),'Replace',false);
    idx = ~ismember(HMAs_temp(:,1),HMAC_source2.powspctrm(:,1));
    HMAI_source2.powspctrm = HMAs_temp(idx,:);
    HMSs_temp = cat(1, HMSC_source2.powspctrm, HMSI_source2.powspctrm);
    HMSC_source2.powspctrm = datasample(HMSs_temp,size(HMSC_source2.trialinfo,1),'Replace',false);
    idx = ~ismember(HMSs_temp(:,1),HMSC_source2.powspctrm(:,1));
    HMSI_source2.powspctrm = HMSs_temp(idx,:);
    
    clear idx NMAs_temp NMSs_temp HMAs_temp HMSs_temp

    for it = 1:iterations
    
    if ~ismember(s,[18 26])
        samp = min([size(NMAI_source2.trialinfo,1) size(NMAC_source2.trialinfo,1)]);
        trl_tmp = sortrows(datasample(1:size(NMAC_source2.trialinfo,1),samp,'Replace',false)','ascend');
        NMAC_tmp = squeeze(nanmean(nanmean(NMAC_source2.powspctrm(trl_tmp,:),2),1));
        trl_tmp = sortrows(datasample(1:size(NMAI_source2.trialinfo,1),samp,'Replace',false)','ascend');
        NMAI_tmp = squeeze(nanmean(nanmean(NMAI_source2.powspctrm(trl_tmp,:),2),1));
        NMA_HvM = (NMAC_tmp - NMAI_tmp);
        allsubj_grid_delta_perm.NMA_HvM(s,it) = NMA_HvM;
    end
        
    samp = min([size(NMSI_source2.trialinfo,1) size(NMSC_source2.trialinfo,1)]);
    trl_tmp = sortrows(datasample(1:size(NMSC_source2.trialinfo,1),samp,'Replace',false)','ascend');
    NMSC_tmp = squeeze(nanmean(nanmean(NMSC_source2.powspctrm(trl_tmp,:),2),1));
    trl_tmp = sortrows(datasample(1:size(NMSI_source2.trialinfo,1),samp,'Replace',false)','ascend');
    NMSI_tmp = squeeze(nanmean(nanmean(NMSI_source2.powspctrm(trl_tmp,:),2),1));
    NMS_HvM = (NMSC_tmp - NMSI_tmp);
    allsubj_grid_delta_perm.NMS_HvM(s,it) = NMS_HvM;
    
    samp = min([size(HMAI_source2.trialinfo,1) size(HMAC_source2.trialinfo,1)]);
    trl_tmp = sortrows(datasample(1:size(HMAC_source2.trialinfo,1),samp,'Replace',false)','ascend');
    HMAC_tmp = squeeze(nanmean(nanmean(HMAC_source2.powspctrm(trl_tmp,:),2),1));
    trl_tmp = sortrows(datasample(1:size(HMAI_source2.trialinfo,1),samp,'Replace',false)','ascend');
    HMAI_tmp = squeeze(nanmean(nanmean(HMAI_source2.powspctrm(trl_tmp,:),2),1));
    HMA_HvM = (HMAC_tmp - HMAI_tmp);
    allsubj_grid_delta_perm.HMA_HvM(s,it) = HMA_HvM;
    
    samp = min([size(HMSI_source2.trialinfo,1) size(HMSC_source2.trialinfo,1)]);
    trl_tmp = sortrows(datasample(1:size(HMSC_source2.trialinfo,1),samp,'Replace',false)','ascend');
    HMSC_tmp = squeeze(nanmean(nanmean(HMSC_source2.powspctrm(trl_tmp,:),2),1));
    trl_tmp = sortrows(datasample(1:size(HMSI_source2.trialinfo,1),samp,'Replace',false)','ascend');
    HMSI_tmp = squeeze(nanmean(nanmean(HMSI_source2.powspctrm(trl_tmp,:),2),1));
    HMS_HvM = (HMSC_tmp - HMSI_tmp);
    allsubj_grid_delta_perm.HMS_HvM(s,it) = HMS_HvM;
    
    clear NMA_HvM NMS_HvM HMA_HvM HMS_HvM
    
    end

    %Update waitbar and message;
    waitbar(s/length(subjects),f,sprintf('Running permutations now'));

end

delete(f);

%save permuted delta power;
allsubj_grid_delta_perm.NMA_HvM(18,:) = mean(allsubj_grid_delta_perm.NMA_HvM([1 2 4:11 13:17 19 21:25],:));
allsubj_grid_delta_perm.NMA_HvM(26,:) = allsubj_grid_delta_perm.NMA_HvM(18,:);

bad_subjects = [3 12 20];
allsubj_grid_delta_perm.NMA_HvM(bad_subjects,:) = []; 
allsubj_grid_delta_perm.NMS_HvM(bad_subjects,:) = []; 
allsubj_grid_delta_perm.HMA_HvM(bad_subjects,:) = []; 
allsubj_grid_delta_perm.HMS_HvM(bad_subjects,:) = []; 

allsubj_grid_delta_HvM_perm = allsubj_grid_delta_perm; clear allsubj_grid_delta_perm

% save mean difference of delta power Hit versus Miss permutations;
save allsubj_Left_Motor_delta_HvM_perm allsubj_grid_delta_HvM_perm

clc; disp([num2str(iterations) 'PERMUTATIONS DONE AND SAVED!']);

%% source grandaverage across participants
clearvars;clc; close all; addpath C:\toolbox\fieldtrip\; ft_defaults;
load allsubj_Left_Motor_delta_HvM; load allsubj_Left_Motor_delta_HvM_perm;

%1)effect size data ori;
Esize_ori = []; Cohen_d = [];
for i = 1:4
    [~,p,~,stats] = ttest(allsubj_grid_delta_HvM(:,i),0,'tail','both');
    Esize_ori(1,i) = stats.tstat(1,1);
    Esize_ori(2,i) = p;
    clear p stats
    %Cohen's d;
    Cohen_d(1,i) = abs(round(mean(allsubj_grid_delta_HvM(:,i))/std(allsubj_grid_delta_HvM(:,i)),3));
end

%2)effect size data permuted;
Esize_perm = [];

for it = 1:size(allsubj_grid_delta_perm.NMA_HvM,2)
   
    [~,p,~,stats] = ttest(allsubj_grid_delta_perm.NMA_HvM(:,it),0,'tail','both');
    Esize_perm.nma(it,1) = stats.tstat(1,1);
    Esize_perm.nma(it,2) = p;
    clear p stats
    
    [~,p,~,stats] = ttest(allsubj_grid_delta_perm.NMS_HvM(:,it),0,'tail','both');
    Esize_perm.nms(it,1) = stats.tstat(1,1);
    Esize_perm.nms(it,2) = p;
    clear p stats
    
    [~,p,~,stats] = ttest(allsubj_grid_delta_perm.HMA_HvM(:,it),0,'tail','both');
    Esize_perm.hma(it,1) = stats.tstat(1,1);
    Esize_perm.hma(it,2) = p;
    clear p stats
    
    [~,p,~,stats] = ttest(allsubj_grid_delta_perm.HMS_HvM(:,it),0,'tail','both');
    Esize_perm.hms(it,1) = stats.tstat(1,1);
    Esize_perm.hms(it,2) = p;
    clear p stats

end

%sort p-values of permutations;
Esize_perm.nma = sortrows(Esize_perm.nma,1,'descend');
Esize_perm.nms = sortrows(Esize_perm.nms,1,'descend');
Esize_perm.hma = sortrows(Esize_perm.hma,1,'descend');
Esize_perm.hms = sortrows(Esize_perm.hms,1,'descend');

%Now compare the original size effect to the permuted data and get p-value;
P_value_perm_nma = ((sum(abs(Esize_perm.nma(:,1)) > abs(Esize_ori(1,1)))) +1)/(size(Esize_perm.nma,1)+1);
P_value_perm_nms = ((sum(abs(Esize_perm.nms(:,1)) > abs(Esize_ori(1,2)))) +1)/(size(Esize_perm.nms,1)+1);
P_value_perm_hma = ((sum(abs(Esize_perm.hma(:,1)) > abs(Esize_ori(1,3)))) +1)/(size(Esize_perm.hma,1)+1);
P_value_perm_hms = ((sum(abs(Esize_perm.hms(:,1)) > abs(Esize_ori(1,4)))) +1)/(size(Esize_perm.hms,1)+1);

display(['P_value_perm_nms: ' num2str(P_value_perm_nms)]);
display(['P_value_perm_nma: ' num2str(P_value_perm_nma)]);
display(['P_value_perm_hms: ' num2str(P_value_perm_hms)]);
display(['P_value_perm_hma: ' num2str(P_value_perm_hma)]);

%% Interaction of effect size;

inter_ori(:,1) = allsubj_grid_delta_HvM(:,1) - allsubj_grid_delta_HvM(:,2);
inter_ori(:,2) = allsubj_grid_delta_HvM(:,3) - allsubj_grid_delta_HvM(:,4);
inter_perm.NMs = allsubj_grid_delta_perm.NMA_HvM - allsubj_grid_delta_perm.NMS_HvM;
inter_perm.HMs = allsubj_grid_delta_perm.HMA_HvM - allsubj_grid_delta_perm.HMS_HvM;

%1)effect size interaction data ori;
Esize_inter_ori = [];
[~,p,~,stats] = ttest(inter_ori(:,1),inter_ori(:,2),'tail','both');
Esize_inter_ori(1,1) = stats.tstat(1,1);
Esize_inter_ori(2,1) = p;
clear p stats

%2)effect size interaction permuted;
Esize_inter_perm = [];

for it = 1:size(inter_perm.NMs,2)
    [~,p,~,stats] = ttest(inter_perm.NMs(:,it),inter_perm.HMs(:,it),'tail','both');
    Esize_inter_perm(it,1) = stats.tstat(1,1);
    Esize_inter_perm(it,2) = p;
    clear p stats
end

%sort p-values of permutations;
Esize_inter_perm = sortrows(Esize_inter_perm,1,'descend');
%Now compare the original size effect to the permuted data and get p-value;
P_value_inter_perm = ((sum(abs(Esize_inter_perm(:,1)) > abs(Esize_inter_ori(1,1))))+1)/(size(Esize_inter_perm,1)+1);
display(['P_value_inter_perm: ' num2str(P_value_inter_perm)]);

%% Plot data;
mean_diff = [nanmean(allsubj_grid_delta_HvM(:,[2 1])); nanmean(allsubj_grid_delta_HvM(:,[4 3]))];
eblow = [nanstd(allsubj_grid_delta_HvM(:, [2 1]))./sqrt(size(allsubj_grid_delta_HvM(:,1),1)); nanstd(allsubj_grid_delta_HvM(:,[4 3]))./sqrt(size(allsubj_grid_delta_HvM(:,1),1))];

%Plot Hits; 
close all; figure;
x1 = bar(mean_diff,1); 
%compute bar centers;
xBar = (cell2mat(get(x1,'XData')).' + [x1.XOffset]) - 0.04;
hold on
er = errorbar(xBar,mean_diff,eblow,'k.','LineWidth',1);
%configure color bars;
x1(1).FaceColor = [0.4 0.4 0.4];
x1(2).FaceColor = [0.7 0.7 .7];
%configure axis;
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.XAxis.Label.Visible ='on';
ax.XTick = [1 1 2 2] + [[x1.XOffset]';[x1.XOffset]']';
ax.XTickLabel = {'NMS','NMA','HMS','HMA'};
ax.YAxis.Limits = [-0.2 0.2];
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 15;
ylabel('Delta power difference');
hold on
plot(ax.XTick(1:4),allsubj_grid_delta_HvM(:,[2 1 4 3]),'o','Color', [0 0 0],'MarkerSize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0]);
set(gcf,'color','w');

%%