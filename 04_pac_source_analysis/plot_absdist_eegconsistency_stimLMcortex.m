% PLOT THE ABSOLUTE DIFFERENCE BETWEEN DELTA POWER AND STIMULI;
clearvars; clc;

load Peak_delta_stim; load Peak_freq_stim
pk_NMSC = Peak_delta_stim.NMSC(:,2:end);
pk_HMSC = Peak_delta_stim.HMSC(:,2:end);
pk_NMAC = Peak_delta_stim.NMAC(:,2:end);
pk_HMAC = Peak_delta_stim.HMAC(:,2:end);

%compute the absolute difference between stimulus peak and delta power peak of each participant;
diff_peaks_Body = []; diff_peaks_Full = []; diff_peaks_Head = []; diff_peaks_Audio = [];
for ii = 1:size(pk_NMSC,2)
diff_peaks_Body.NMSC(:,ii) = abs(Peak_freq_stim.Body(Peak_freq_stim.Condition==2) - pk_NMSC(:,ii));
diff_peaks_Body.NMAC(:,ii) = abs(Peak_freq_stim.Body(Peak_freq_stim.Condition==2) - pk_NMAC(:,ii));
diff_peaks_Body.HMSC(:,ii) = abs(Peak_freq_stim.Body(Peak_freq_stim.Condition==4) - pk_HMSC(:,ii));
diff_peaks_Body.HMAC(:,ii) = abs(Peak_freq_stim.Body(Peak_freq_stim.Condition==4) - pk_HMAC(:,ii));

diff_peaks_Full.NMSC(:,ii) = abs(Peak_freq_stim.Full(Peak_freq_stim.Condition==2) - pk_NMSC(:,ii));
diff_peaks_Full.NMAC(:,ii) = abs(Peak_freq_stim.Full(Peak_freq_stim.Condition==2) - pk_NMAC(:,ii));
diff_peaks_Full.HMSC(:,ii) = abs(Peak_freq_stim.Full(Peak_freq_stim.Condition==4) - pk_HMSC(:,ii));
diff_peaks_Full.HMAC(:,ii) = abs(Peak_freq_stim.Full(Peak_freq_stim.Condition==4) - pk_HMAC(:,ii));

diff_peaks_Head.NMSC(:,ii) = abs(Peak_freq_stim.Head(Peak_freq_stim.Condition==2) - pk_NMSC(:,ii));
diff_peaks_Head.NMAC(:,ii) = abs(Peak_freq_stim.Head(Peak_freq_stim.Condition==2) - pk_NMAC(:,ii));
diff_peaks_Head.HMSC(:,ii) = abs(Peak_freq_stim.Head(Peak_freq_stim.Condition==4) - pk_HMSC(:,ii));
diff_peaks_Head.HMAC(:,ii) = abs(Peak_freq_stim.Head(Peak_freq_stim.Condition==4) - pk_HMAC(:,ii));

diff_peaks_Audio.NMSC(:,ii) = abs(Peak_freq_stim.Audio(Peak_freq_stim.Condition==2) - pk_NMSC(:,ii));
diff_peaks_Audio.NMAC(:,ii) = abs(Peak_freq_stim.Audio(Peak_freq_stim.Condition==2) - pk_NMAC(:,ii));
diff_peaks_Audio.HMSC(:,ii) = abs(Peak_freq_stim.Audio(Peak_freq_stim.Condition==4) - pk_HMSC(:,ii));
diff_peaks_Audio.HMAC(:,ii) = abs(Peak_freq_stim.Audio(Peak_freq_stim.Condition==4) - pk_HMAC(:,ii));

end

%average difference scores across participants for each stimulus/info;
avg_diff = [];
avg_diff(:,1) = nanmean(diff_peaks_Full.NMSC,2);
avg_diff(:,2) = nanmean(diff_peaks_Full.NMAC,2);
avg_diff(:,3) = nanmean(diff_peaks_Head.NMSC,2);
avg_diff(:,4) = nanmean(diff_peaks_Head.NMAC,2);
avg_diff(:,5) = nanmean(diff_peaks_Body.NMSC,2);
avg_diff(:,6) = nanmean(diff_peaks_Body.NMAC,2);
avg_diff(:,7) = nanmean(diff_peaks_Audio.NMSC,2);
avg_diff(:,8) = nanmean(diff_peaks_Audio.NMAC,2);
avg_diff(:,9) = nanmean(diff_peaks_Full.HMSC,2);
avg_diff(:,10) = nanmean(diff_peaks_Full.HMAC,2);
avg_diff(:,11) = nanmean(diff_peaks_Head.HMSC,2);
avg_diff(:,12) = nanmean(diff_peaks_Head.HMAC,2);
avg_diff(:,13) = nanmean(diff_peaks_Body.HMSC,2);
avg_diff(:,14) = nanmean(diff_peaks_Body.HMAC,2);
avg_diff(:,15) = nanmean(diff_peaks_Audio.HMSC,2);
avg_diff(:,16) = nanmean(diff_peaks_Audio.HMAC,2);

%perform one-sample ttests against zeros;
ttest_diff = [];
for ii = 1:size(avg_diff,2)
    [~,ttest_diff(1,ii),~,stats] = ttest(avg_diff(:,ii),0,'tail','right');
    ttest_diff(2,ii) = round(stats.tstat(1),2);
end

%correct for multiple comparisons;
for ii = 1:size(ttest_diff,2)
    if ttest_diff(1,ii) <= 0.05/16 
        ttest_diff(3,ii) = 1; 
    else
        ttest_diff(3,ii)= 0;
    end
end

% save absolute ditance between stimuli and delta Left Motor;
% save absdist_deltapeak_stimLMC avg_diff ttest_diff

% PLOT Absolute distance;
mean_abs_diff = [mean(avg_diff(:,1:8));mean(avg_diff(:,9:16))];
sem_abs_diff = [std(avg_diff(:,1:8))./sqrt(size(avg_diff,1));std(avg_diff(:,9:16))./sqrt(size(avg_diff,1))];

close all; figure;
x1 = bar(mean_abs_diff,1); 
hold on
xBar = (cell2mat(get(x1,'XData')).' + [x1.XOffset]);
er = errorbar(xBar,mean_abs_diff,sem_abs_diff,'k.','LineWidth',1);

%configure color bars;
x1(1).FaceColor = [0.4 0.4 0.4];
x1(2).FaceColor = [0.8 0.8 0.8];
x1(3).FaceColor = [0.4 0.4 0.4];
x1(4).FaceColor = [0.8 0.8 0.8];
x1(5).FaceColor = [0.4 0.4 0.4];
x1(6).FaceColor = [0.8 0.8 0.8];
x1(7).FaceColor = [0.4 0.4 0.4];
x1(8).FaceColor = [0.8 0.8 0.8];

%configure axis;
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.XTick = [];
ax.YAxis.Limits = [0.6 2.2];
ax.YTick = 0.6:0.2:2.2;
ax.YAxis.FontSize = 12;
ylabel('Absolute difference');
ax.YAxis.Label.FontSize = 14;
set(gcf,'color','w','Position', [329 159.40  724.8  569.6]);

%% Plot Consistency of delta peaks across the ordered stimuli in all four conditions.

body_color = [0 0.2 0.6];
full_color = [0 0.4 0.6];
head_color = [0 0.6 0.6];
audio_color = [0 0.8 0.6];
sync_condition = [0 0 0];
async_condition = [1 0 0];

close all; figure;hold on
subplot(2,1,1)
for ii = 1:size(Peak_freq_stim(Peak_freq_stim.Condition==2,:),1)
    line([ii-0.1 ii-0.1],[nanmean(Peak_delta_stim.NMAC(ii,2:end)')-nanstd(Peak_delta_stim.NMAC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end))) nanmean(Peak_delta_stim.NMAC(ii,2:end)')+ nanstd(Peak_delta_stim.NMAC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end)))],'Color',[0.5 0.5 0.5]);   
    line([ii ii],[nanmean(Peak_delta_stim.NMSC(ii,2:end)')-nanstd(Peak_delta_stim.NMSC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end))) nanmean(Peak_delta_stim.NMSC(ii,2:end)')+nanstd(Peak_delta_stim.NMSC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end)))],'Color',[0.5 0.5 0.5]);
    scatter(ii-0.1,nanmean(Peak_delta_stim.NMAC(ii,2:end)'),20,'s','filled','MarkerFaceColor',async_condition);
    hold on
    scatter(ii,nanmean(Peak_delta_stim.NMSC(ii,2:end)'),20,'s','filled','MarkerFaceColor',sync_condition);
end

%configure axis;
xlim([0 55]); xticks([]); 
ylim([1 3]); yticks([1 2 3]);
ylabel('Neural peak frequency [Hz]');
title('No-mask');
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.YAxis.FontSize = 12;
ax.YAxis.Label.FontSize = 12;

subplot(2,1,2)
for ii = 1:size(Peak_freq_stim(Peak_freq_stim.Condition==4,:),1)
    line([ii-0.1 ii-0.1],[nanmean(Peak_delta_stim.HMAC(ii,2:end)')-nanstd(Peak_delta_stim.HMAC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end))) nanmean(Peak_delta_stim.HMAC(ii,2:end)')+nanstd(Peak_delta_stim.HMAC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end)))],'Color',[0.5 0.5 0.5]);   
    line([ii ii],[nanmean(Peak_delta_stim.HMSC(ii,2:end)')-nanstd(Peak_delta_stim.HMSC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end))) nanmean(Peak_delta_stim.HMSC(ii,2:end)')+nanstd(Peak_delta_stim.HMSC(ii,2:end)')/sqrt(length(Peak_delta_stim.NMAC(ii,2:end)))],'Color',[0.5 0.5 0.5]);
    scatter(ii-0.1,nanmean(Peak_delta_stim.HMAC(ii,2:end)'),20,'s','filled','MarkerFaceColor',async_condition);
    hold on
    scatter(ii,nanmean(Peak_delta_stim.HMSC(ii,2:end)'),20,'s','filled','MarkerFaceColor',sync_condition);
end

%configure axis;
xlim([0 55]); xticks([]); 
ylim([1 3]); yticks([1 2 3]);
ylabel('Neural peak frequency [Hz]');
title('Head-mask');
xlabel('Ordered trials')
ax = gca;
ax.Box = 'off';
ax.LineWidth = 1.5;
ax.TickDir = 'out';
ax.YAxis.FontSize = 12;
ax.YAxis.Label.FontSize = 12;
set(gcf,'color','w','Position', [329 159.40  724.8  569.6]);

%%