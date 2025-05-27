%% 1- Create the logfile with all the trials together(- the first ten trials);
clearvars; clc; 

% Compute the CR,RTs and Dprimes for all participants together;
allsubj_scores = cell(25,12);
allsubj_scores = array2table(allsubj_scores, 'Variable',{'participant','cr_NMA','cr_NMS','cr_HMA','cr_HMS','rt_NMA','rt_NMS','rt_HMA','rt_HMS','Dprime_NMCs','Dprime_HMCs','Dprime_SyncAsync'}); 

subjects = [1 2 4:26];

for s = subjects
    
    logfile= readtable(['xxx\Logfile_EEG_',num2str(s),'.xls']);

    %count hits/misses (NMS/HMS); 
    HIT_MISS(1,1) = sum(logfile.condition_code ==2 & logfile.hits ==1);
    HIT_MISS(1,2) = sum(logfile.condition_code ==2 & logfile.hits ==0);
    HIT_MISS(2,1) = sum(logfile.condition_code ==4 & logfile.hits ==1);
    HIT_MISS(2,2) = sum(logfile.condition_code ==4 & logfile.hits ==0);
    %count CRj/FA (NMA/HMA);
    CRj_FA(1,1) = sum(logfile.condition_code ==1 & logfile.hits ==1);
    CRj_FA(1,2) = sum(logfile.condition_code ==1 & logfile.hits ==0);
    CRj_FA(2,1) = sum(logfile.condition_code ==3 & logfile.hits ==1);
    CRj_FA(2,2) = sum(logfile.condition_code ==3 & logfile.hits ==0);
    %Correct Responses rates and their RTs;
    CR_NMA(1,1) = sum(logfile.condition_code ==1 & logfile.hits ==1)/54;
    CR_NMS(1,1) = sum(logfile.condition_code ==2 & logfile.hits ==1)/54;
    CR_HMA(1,1) = sum(logfile.condition_code ==3 & logfile.hits ==1)/54;
    CR_HMS(1,1) = sum(logfile.condition_code ==4 & logfile.hits ==1)/54;
    CR_NMA_rt(1,1) = mean(logfile.reaction_time(logfile.condition_code ==1 & logfile.hits ==1 & logfile.reaction_time < (mean(logfile.reaction_time) + 2*std(logfile.reaction_time)) & logfile.reaction_time > (mean(logfile.reaction_time) - 2*std(logfile.reaction_time))));
    CR_NMS_rt(1,1) = mean(logfile.reaction_time(logfile.condition_code ==2 & logfile.hits ==1 & logfile.reaction_time < (mean(logfile.reaction_time) + 2*std(logfile.reaction_time)) & logfile.reaction_time > (mean(logfile.reaction_time) - 2*std(logfile.reaction_time))));
    CR_HMA_rt(1,1) = mean(logfile.reaction_time(logfile.condition_code ==3 & logfile.hits ==1 & logfile.reaction_time < (mean(logfile.reaction_time) + 2*std(logfile.reaction_time)) & logfile.reaction_time > (mean(logfile.reaction_time) - 2*std(logfile.reaction_time))));
    CR_HMS_rt(1,1) = mean(logfile.reaction_time(logfile.condition_code ==4 & logfile.hits ==1 & logfile.reaction_time < (mean(logfile.reaction_time) + 2*std(logfile.reaction_time)) & logfile.reaction_time > (mean(logfile.reaction_time) - 2*std(logfile.reaction_time))));

    %Calculate the Dprimes in the 2 masks NMCs and HMCs;
    rate_hit_NMS = HIT_MISS(1,1)/54;
    rate_fa_NMA = CRj_FA(1,2)/54;
    rate_hit_HMS = HIT_MISS(2,1)/54;
    rate_fa_HMA = CRj_FA(2,2)/54;
    
    %if hits = 1, put 0.999 and fa = 0 put 0.001 to be able to compute the dprime;
    if rate_hit_NMS == 1; rate_hit_NMS = 0.999; end
    if rate_fa_NMA == 0; rate_fa_NMA = 0.001; end
    if rate_hit_HMS == 1; rate_hit_HMS = 0.999; end
    if rate_fa_HMA == 0; rate_fa_HMA = 0.001; end
    Dprime_NMCs = dprime_simple(rate_hit_NMS,rate_fa_NMA);
    Dprime_HMCs = dprime_simple(rate_hit_HMS,rate_fa_HMA);

    %Calculate the Dprimes Synchrony/Asynchrony (to see whether asynchrony bias response in generale);
    %count hits/misses SYNC; 
    HIT_MISS_SYNC(1,1) = sum(logfile.condition_code ==2 | logfile.condition_code ==4 & logfile.hits ==1);
    HIT_MISS_SYNC(1,2) = sum(logfile.condition_code ==2 | logfile.condition_code ==4 & logfile.hits ==0);
    %count CRj/FA ASYNC;
    CRj_FA_ASYNC(1,1) = sum(logfile.condition_code==1 | logfile.condition_code==3 & logfile.hits ==1);
    CRj_FA_ASYNC(1,2) = sum(logfile.condition_code ==1 | logfile.condition_code ==3 & logfile.hits ==0);

    %Dprimes Sync;
    rate_hit_Sync = HIT_MISS_SYNC(1,1)/108;
    rate_fa_Async = CRj_FA_ASYNC(1,2)/108;
    Dprime_SyncAsync = dprime_simple(rate_hit_Sync,rate_fa_Async);
    
    % writetable(logfile, ['Logfile_EEG_', num2str(s), '.xls']);
    save(['xxx\EEG_',num2str(s),'_scores'], 'HIT_MISS', 'CRj_FA','CR_NMA','CR_NMS','CR_HMA','CR_HMS','CR_NMA_rt','CR_NMS_rt','CR_HMA_rt','CR_HMS_rt','Dprime_NMCs','Dprime_HMCs','Dprime_SyncAsync');
       
    allsubj_scores.participant{s} = s;
    allsubj_scores.cr_NMA{s} = CR_NMA;
    allsubj_scores.cr_NMS{s} = CR_NMS;
    allsubj_scores.cr_HMA{s} = CR_HMA;
    allsubj_scores.cr_HMS{s} = CR_HMS;
    allsubj_scores.rt_NMA{s} = CR_NMA_rt;
    allsubj_scores.rt_NMS{s} = CR_NMS_rt;
    allsubj_scores.rt_HMA{s} = CR_HMA_rt;
    allsubj_scores.rt_HMS{s} = CR_HMS_rt;
    allsubj_scores.Dprime_NMCs{s} = Dprime_NMCs;
    allsubj_scores.Dprime_HMCs{s} = Dprime_HMCs;
    allsubj_scores.Dprime_SyncAsync{s} = Dprime_SyncAsync;

end

empty_rows = any(cellfun(@isempty, allsubj_scores.participant), 2);
allsubj_scores(empty_rows,:) = [];

%save allsubject behavioural performances; 
save allsubj_scores allsubj_scores

%% Plot behavioural performances of allsubjects;

%remove the bad subjects before;
bad_subj = [3 12 20];
allsubj_behavior_scores(bad_subj,:) = [];
hits = [nanmean(cell2mat(allsubj_behavior_scores{:,[3 2]})); nanmean(cell2mat(allsubj_behavior_scores{:,[5 4]}))];
sem_hits = [nanstd(cell2mat(allsubj_behavior_scores{:,[3 2]}))/sqrt(size(allsubj_behavior_scores,1));nanstd(cell2mat(allsubj_behavior_scores{:,[5 4]}))/sqrt(size(allsubj_behavior_scores,1))];

%Plot Hits; 
figure;
x1 = bar(hits,1); 
%compute bar centers;
xBar = (cell2mat(get(x1,'XData')).' + [x1.XOffset]) - 0.04;
hold on
er = errorbar(xBar,hits,sem_hits,'k.','LineWidth',1);
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
ax.XTickLabel = {'NMSC','NMAC','HMSC','HMAC'};
ax.XAxis.FontSize = 12;
ax.YAxis.Limits = [0 1];
ax.YTick = 0:0.2:1;
ax.YTickLabel = [0 0.2 0.4 0.6 0.8 1]*100;
ax.YAxis.FontSize = 15;
ylabel('Correct responses (% \pm std)');
hold on
plot(ax.XTick(1:4),cell2mat(allsubj_behavior_scores{:,[3 2 5 4]}),'o','Color', [0 0 0],'MarkerSize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0]);
set(gcf,'color','w');

% PLOT REACTION TIMES;
close all
rt = [nanmean(cell2mat(allsubj_behavior_scores{:,[7 6]})*1000); nanmean(cell2mat(allsubj_behavior_scores{:,[9 8]})*1000)];
sem_rt = [nanstd(cell2mat(allsubj_behavior_scores{:,[7 6]}))/sqrt(size(allsubj_behavior_scores,1)); nanstd(cell2mat(allsubj_behavior_scores{:,[9 8]}))/sqrt(size(allsubj_behavior_scores,1))].*1000;

%Plot Hits; 
figure;
x2 = bar(rt,1); 
%compute bar centers;
xBar = (cell2mat(get(x2,'XData')).' + [x2.XOffset]) - 0.04;
hold on
er = errorbar(xBar,rt,sem_rt,'k.','LineWidth',1);
%configure color bars;
x2(1).FaceColor = [0.4 0.4 0.4];
x2(2).FaceColor = [0.7 0.7 .7];
%configure axis;
ax1 = gca;
ax1.Box = 'off';
ax1.LineWidth = 1.5;
ax1.TickDir = 'out';
ax1.XAxis.Label.Visible ='on';
ax1.XTick = [1 1 2 2] + [[x2.XOffset]';[x2.XOffset]']';
ax1.XTickLabel = {'NMSC','NMAC','HMSC','HMAC'};
ax1.XAxis.FontSize = 12;
ax1.YAxis.Limits = [0 2500];
ax1.YTick = 0:500:2500;
ax1.YAxis.FontSize = 15;
ylabel('Reaction times (ms \pm std)');
hold on
plot(ax1.XTick(1:4),cell2mat(allsubj_behavior_scores{:,[7 6 9 8]})*1000,'o','Color', [0 0 0],'MarkerSize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0]);
hold on
set(gcf,'color','w');

% PLOT DPRIME;
close all
dprime = [nanmean(cell2mat(allsubj_behavior_scores{:,10:11})); nanmean(cell2mat(allsubj_behavior_scores{:,10:11}))];
sem_dprime = [nanstd(cell2mat(allsubj_behavior_scores{:,10:11}))/sqrt(size(allsubj_behavior_scores,1));nanstd(cell2mat(allsubj_behavior_scores{:,10:11}))/sqrt(size(allsubj_behavior_scores,1))];

%Plot Hits; 
figure;
x3 = bar(dprime,1); 
%compute bar centers;
xBar = (cell2mat(get(x3,'XData')).' + [x3.XOffset]) - 0.04;
hold on
er = errorbar(xBar,dprime,sem_dprime,'k.','LineWidth',1);
%configure color bars;
x3(1).FaceColor = [0.4 0.4 0.4];
x3(2).FaceColor = [0.7 0.7 .7];
%configure axis;
ax3 = gca;
ax3.Box = 'off';
ax3.LineWidth = 1.5;
ax3.TickDir = 'out';
ax3.XAxis.Label.Visible ='on';
ax3.XAxis.Limits = [0.5 1.5];
ax3.XTick = [1 1 2 2] + [[x3.XOffset]';[x3.XOffset]']';
ax3.XTickLabel = {'Diff. NMCs','Diff. HMCs'};
ax3.XAxis.FontSize = 12;
ax3.YAxis.Limits = [0 3];
ax3.YTick = 0:0.6:3;
ax3.YAxis.FontSize = 15;
ylabel('Dprime (\pm std)');
hold on
plot(ax3.XTick(1:2),cell2mat(allsubj_behavior_scores{:,10:11}),'o','Color', [0 0 0],'MarkerSize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor', [0 0 0]);
set(gcf,'color','w');

%%