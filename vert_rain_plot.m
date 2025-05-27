%% This function creates a vertical raincloud plot which is sooo fancy these days ...

% Input: 
% dat=NxM matrix containing the data to plot; rows=number of
% observations; columns = number of conditions
% col = color in rgb format (default: [0.8 0.2 0.5]);
% dotsize = size of dots to be plotted (default: 50)
% smooth = width of kernel used for smoothing (default: 100)

function vert_rain_plot(dat,col,dotsize,smooth)

ncond = size(dat,2);
if isempty(col)
    col = [0 1 0];
end

if isempty(dotsize)
    dotsize = 50;
end

if isempty(smooth)
    smooth=100;
end



x = linspace(0.4,ncond+0.6);
y = linspace(0,0);
% plot(x,y,'linewidth',1,'color',[0.5 0.5 0.5],'linestyle', '--');
hold on

for k=1:ncond
    tmp=dat(:,k);
    jitt=(rand(1,length(tmp))+0.2).*0.2;
    scatter(k+jitt,tmp,dotsize,col,'filled','MarkerFaceAlpha',.4);
    hold on
end

xlim([0.4 max(x)]);


% Plot Boxes and means ...
xb=[0.865 0.98];
xbr=[0.865 0.98];
for plt=1:ncond
    Mns=nanmean(dat(:,plt));
    STE=nanstd(dat(:,plt))./sqrt(length(dat));% 
    STD=nanstd(dat(:,plt));% 
    
    % draw box for STD
    tmpxb=[xb(1,1) xb(1,1) xb(1,2) xb(1,2) xb(1,1)]+plt-1;
    tmpyb=[Mns-STD Mns+STD Mns+STD Mns-STD Mns-STD];
    plot(tmpxb,tmpyb,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    % draw filled box for STE
    tmpxb=[xb(1,1) xb(1,1) xb(1,2) xb(1,2)]+plt-1;
    tmpyb=[Mns-STE Mns+STE Mns+STE Mns-STE];
    fill(tmpxb,tmpyb,col,'FaceAlpha',.4,'EdgeAlpha',0);
    hold on
    
    % Plot solid Lines for means
    x = linspace(xb(1,1)+plt,xb(1,2)+plt)-1;
    y = linspace(mean(Mns),mean(Mns));

    plot(x,y,'linewidth',1,'color',col,'linestyle', '-','LineWidth',2);
    hold on
    
    % Plot Error bars for 90th percentile
    pct = prctile(dat(:,plt),[5 95]);
    x = linspace(xbr(1,1)+plt,xbr(1,2)+plt)-1;
    y = linspace(pct(1),pct(1));
    % draw lower bound
    plot(x,y,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
    x = linspace(xbr(1,1)+plt,xbr(1,2)+plt)-1;
    y = linspace(pct(2),pct(2));
    % draw upper bound
    plot(x,y,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
    % draw connecting lines to upper bound
    x = [mean(xbr,2),mean(xbr,2)]+plt-1;
    y = [Mns+STD,pct(2)];
    % draw upper bound
    plot(x,y,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
    % draw connecintg lines to lower bound
    x = [mean(xbr,2),mean(xbr,2)]+plt-1;
    y = [Mns-STD,pct(1)];
    % draw upper bound
    plot(x,y,'linewidth',1,'color',col,'linestyle', '-','LineWidth',1.5);
    hold on
    
end

% Let it rain ...
nbins=400;
scal=0.3;
xpos=0.5;
%xb=[0.2 0.2+scal;1.2 1.2+scal;2.2 2.2+scal];

for plt=1:ncond
    histdat=hist(dat(:,plt),nbins);
    dens=conv(histdat,gausswin(smooth),'same')./(smooth/1000);% compute smoothed firing rate using gaussian kernel
    %dens=dens-min(dens);
    dens=(dens./(max(dens))).*-scal;
    tmpxb=[0 dens 0 0]-min(dens)+plt-1+xpos;
    tmpyb=[min(dat(:,plt)) linspace(min(dat(:,plt)),max(dat(:,plt)),nbins) max(dat(:,plt)) min(dat(:,plt))];
    %figure;plot(tmpxb(1:end),tmpyb(1:end));
    fill(tmpxb,tmpyb,col,'FaceAlpha',.4,'EdgeAlpha',0);
    hold on
end

%ax=gca;
%ax.YLim=[-4 6];