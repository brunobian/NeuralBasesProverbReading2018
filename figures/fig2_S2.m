%% Load
clear all
close all
load('../data/EEG') % load EEG pre-analized data
addpath('../functions/')
addpath('/home/brunobian/Documents/Repos/fieldtrip/')
ft_defaults
%% Fig 2 (a-c)     - Results - ERPs - Pred
runStats = 1;
plotTopo = 1;

% ERP
close all
x_lim=[-100 700]; y_lim=[-3.5 3.5];

figure();clf
t= erp.times.clasicos;
set(gcf,'Color','w')
    plot([0 0],y_lim,'k-'); 
    hold on
    plot(x_lim,[0 0],'k-');
    plot([win.N400.time(1)  win.N400.time(1) ],y_lim,'k--'); 
    plot([win.N400.time(2)  win.N400.time(2) ],y_lim,'k--'); 
     
    [~, h1] = niceBars2(t, erp.pred_ROIN400.q1.m, erp.pred_ROIN400.q1.e, [1 0 0], .3,'-');
    [~, h2] = niceBars2(t, erp.pred_ROIN400.q2.m, erp.pred_ROIN400.q2.e, [0 .5 0], .3,'-');
    [~, h3] = niceBars2(t, erp.pred_ROIN400.q3.m, erp.pred_ROIN400.q3.e, [0 0 1], .3,'-');

    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    legend([h1, h2, h3], 'Low Pred','Mid Pred','High Pred')

    % Topoplot hi-lo N400
    times         = win.N400.time;
    limites_color = [-3 3];
    elec_marks    = win.N400.elect;    
    t_ind         = (erp.times.clasicos >= times(1)) & ...
                    (erp.times.clasicos <  times(2));

    amps =  erp.pred.q1.m(:,t_ind) - erp.pred.q3.m(:,t_ind);
    figure();clf;set(gcf,'Color','w')
        topoplot(mean(amps,2), ...
                 CHANS.chanlocs,'maplimits',limites_color);
        
    figure();clf;set(gcf,'Color','w')
        topoplot(zeros(size(amps,1),1), ...
                 CHANS.chanlocs,'maplimits',limites_color, ...
                 'emarker2', {elec_marks,'o','k',10,1});
            
if plotTopo 
    % Topoplot N400
    x = -.2;
    y = .8;
    times         = win.N400.time;
    limites_color = [-3 3];
    elec_marks    = win.N400.elect;    
    t_ind         = (erp.times.clasicos >= times(1)) & ...
                    (erp.times.clasicos <  times(2));

    figure();clf
    set(gcf,'Color','w')
        subplot(1, 3, 1)
            topoplot(mean(erp.pred.q1.m(:,t_ind),2), ...
                     CHANS.chanlocs,'maplimits',limites_color, ...
                     'emarker2', {elec_marks,'o','k',4,1});
            text(x,y, 'Low Pred')

        subplot(1, 3, 2)
            topoplot(mean(erp.pred.q2.m(:,t_ind),2), ...
                     CHANS.chanlocs,'maplimits',limites_color, ...
                     'emarker2', {elec_marks,'o','k',4,1});
            text(x,y, 'Mid Pred')

        subplot(1, 3, 3)
            topoplot(mean(erp.pred.q3.m(:,t_ind),2), ...
                     CHANS.chanlocs,'maplimits',limites_color, ...
                     'emarker2', {elec_marks,'o','k',4,1});
            text(x,y, 'Hi Pred')

end

% Stats
if runStats
    M = [erp.stats.pred.q1' erp.stats.pred.q2' erp.stats.pred.q3'];
    [p,tbl,stats] = kruskalwallis(M);
end
%% Fig 2 (b-d)     - Results - ERPs - SntcType
runStats = 1;

% ERP
close all
x_lim=[-100 700];
y_lim=[-3.5 3.5];

figure();clf
t= erp.times.clasicos;
set(gcf,'Color','w')
    plot([0 0],y_lim,'k-'); 
    hold on
    plot(x_lim,[0 0],'k-');
    plot([win.N400.time(1)  win.N400.time(1) ],y_lim,'k--'); 
    plot([win.N400.time(2)  win.N400.time(2) ],y_lim,'k--'); 
    
    [~, h1] = niceBars2(t, erp.sntcBinary_ROIN400.proverb.m, erp.sntcBinary_ROIN400.proverb.e, [1 0 0], .3,'-');
    [~, h2] = niceBars2(t, erp.sntcBinary_ROIN400.noProverb.m, erp.sntcBinary_ROIN400.noProverb.e, [0 .5 0], .3,'-');
    
    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    legend([h1, h2], 'Proverb','noProverb')


% Topoplot N400
times         = win.N400.time;
limites_color = [-3 3];
elec_marks    = win.N400.elect;    
t_ind         = (erp.times.clasicos >= times(1)) & ...
                (erp.times.clasicos <  times(2));

figure();clf
set(gcf,'Color','w')
    topoplot(mean(erp.sntcBinary.proverb.m(:,t_ind),2), ...
                 CHANS.chanlocs,'maplimits',limites_color);

figure();clf
set(gcf,'Color','w')
    topoplot(mean(erp.sntcBinary.noProverb.m(:,t_ind),2), ...
                 CHANS.chanlocs,'maplimits',limites_color);

% Stats
if runStats

    [p,h,stats1] = ranksum(erp.stats.Sntc.Prov', erp.stats.Sntc.nProv')

end
%% Fig 2 (NotShown)- Results - ERPs - Freq
runStats = 1;
plotTopo = 0;

% ERP
close all
x_lim=[-100 700]; y_lim=[-3.5 3.5];

figure();clf
t= erp.times.clasicos;
set(gcf,'Color','w')
    plot([0 0],y_lim,'k-'); 
    hold on
    plot(x_lim,[0 0],'k-');
    plot([win.N400.time(1)  win.N400.time(1) ],y_lim,'k--'); 
    plot([win.N400.time(2)  win.N400.time(2) ],y_lim,'k--'); 
     
    [~, h1] = niceBars2(t, erp.freq_ROIN400.q1.m, erp.freq_ROIN400.q1.e, [1 0 0], .3,'-');
    [~, h2] = niceBars2(t, erp.freq_ROIN400.q2.m, erp.freq_ROIN400.q2.e, [0 .5 0], .3,'-');
    [~, h3] = niceBars2(t, erp.freq_ROIN400.q3.m, erp.freq_ROIN400.q3.e, [0 0 1], .3,'-');

    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    legend([h1, h2, h3], 'Low freq','Mid freq','High freq')

    % Topoplot resta hi-lo N400
    times         = win.N400.time;
    limites_color = [-3 3];
    elec_marks    = win.N400.elect;    
    t_ind         = (erp.times.clasicos >= times(1)) & ...
                    (erp.times.clasicos <  times(2));

    amps =  erp.freq.q1.m(:,t_ind) - erp.freq.q3.m(:,t_ind);
    figure();clf;set(gcf,'Color','w')
        topoplot(mean(amps,2), ...
                 CHANS.chanlocs,'maplimits',limites_color);
%          colorbar
        
    figure();clf;set(gcf,'Color','w')
        topoplot(zeros(size(amps,1),1), ...
                 CHANS.chanlocs,'maplimits',limites_color, ...
                 'emarker2', {elec_marks,'o','k',10,1});
            
if plotTopo 
    % Topoplot N400
    x = -.2;
    y = .8;
    times         = win.N400.time;
    limites_color = [-3 3];
    elec_marks    = win.N400.elect;    
    t_ind         = (erp.times.clasicos >= times(1)) & ...
                    (erp.times.clasicos <  times(2));

    figure();clf
    set(gcf,'Color','w')
        subplot(1, 3, 1)
            topoplot(mean(erp.freq.q1.m(:,t_ind),2), ...
                     CHANS.chanlocs,'maplimits',limites_color, ...
                     'emarker2', {elec_marks,'o','k',4,1});
            text(x,y, 'Low freq')

        subplot(1, 3, 2)
            topoplot(mean(erp.freq.q2.m(:,t_ind),2), ...
                     CHANS.chanlocs,'maplimits',limites_color, ...
                     'emarker2', {elec_marks,'o','k',4,1});
            text(x,y, 'Mid freq')

        subplot(1, 3, 3)
            topoplot(mean(erp.freq.q3.m(:,t_ind),2), ...
                     CHANS.chanlocs,'maplimits',limites_color, ...
                     'emarker2', {elec_marks,'o','k',4,1});
            text(x,y, 'Hi freq')
end

% Stats
if runStats

    
    M = [erp.stats.freq.q1' erp.stats.freq.q2' erp.stats.freq.q3'];
    [p,tbl,stats] = kruskalwallis(M);

end
%% Fig S2 (a)      - Regression   
close all

Amp  = erp.regresion_pred.m;
Suj  = categorical(erp.regresion_pred.suj);
Pred = erp.regresion_pred.pred;


Amp  = erp.regresion_pred.m;
Suj  = dummyvar(erp.regresion_pred.suj);
Pred = erp.regresion_pred.pred;
b    = regress(Amp, [Pred, Suj]);

Amp2 = b(2) .* Pred + b(1);


figure();
hold on
%scatter(Pred, Amp)
% h=boxplot(Amp, Pred,'PlotStyle','compact', 'Positions', unique(Pred))
h=boxplot(Amp, Pred,'Positions', unique(Pred),'Widths',0.04)
set(h(7,:),'Visible','off')
ylim([-30 30])
xtick = linspace(min(Pred),max(Pred),6)
set(gca, 'xtick', xtick ,'xticklabel', round(xtick,2,'significant') )


plot(Pred, Amp2, 'Color', [0 .5 0], 'LineWidth', 1)
xlim([-1.5 1.5])
ylabel('Voltage')
xlabel('Pred')
%% Fig S2 (b-c)    - CBPT - Pred
v1 = erp.CBPT_pred.q1.data;
v2 = erp.CBPT_pred.q3.data;

ft_defaults;

Y = nan(size(v1 ,1), size(v1,2), 2, size(v1,3));
Y(:,:,1,:) = v1;
Y(:,:,2,:) = v2;

t = erp.times.clasicos;

CONFIG.EEGBLANK         = [];
CONFIG.clusteralpha     = 0.01;     % alpha level of the sample-specific test statistic that will be used for thresholding
CONFIG.minnbchan        = 2;        % minimum number of neighborhood channels that is required for a selected
                                    % sample to be included in the clustering algorithm (default=0).
CONFIG.numrandomization = 500;      % number of draws from the permutation distribution

Nsuj = size(v1,3);
[stat times ftdata grandavg cfg] = my_ftstat(Y,t,[],1:Nsuj,CONFIG,CHANS);

indpos = find([stat.posclusters.prob]<cfg.alpha);
if ~isempty(indpos)
    npos = length(indpos);
    figure;
        z = stat.posclusterslabelmat;
        x = ismember(z,indpos);
        imagesc(t,1:128,x)
end

indneg = find([stat.negclusters.prob]<cfg.alpha);
if ~isempty(indneg)
    nneg = length(indneg);
    figure;
        z = stat.negclusterslabelmat;
        x = ismember(z,indneg);
        imagesc(t,1:128,x)
end
colormap(summer)


% Topoplot of significant cluster in N400 window    
lim = [0 14];
timelim = [200 450];
tval = stat.prob;
ventana = tval(:,t>timelim(1) & t<timelim(2));
meanWindow = nanmean(ventana,2);
meanWindow = sum(ventana<0.05,2);
cmap = flip([repmat(1,14,1), linspace(0,1,14)', linspace(0,1,14)']);

    figure();
        clusnum = 1;
        mask = any(z(:,t>timelim(1) & t<timelim(2))'>=clusnum)';
        window = meanWindow .* mask;
        topoplot(meanWindow .* mask, ...
                 CHANS.chanlocs, ...
                 'maplimits',lim , ...
                 'colormap', colormap(cmap));
        set(gcf,'Color','w')
        set(findobj(gca,'type','patch'),'facecolor',get(gcf,'color'))
        colorbar
        
% Plot histogram with number sig
        figure();
        set(gcf,'Color','w')
        mask = repmat(any(tval<0.05),128,1);
        x = sum(z.*mask, 1);
        plot(t,x, 'LineWidth', 2)
        hold on
            xlim([-105.5 700])
            ylim([0 max(x)+5])
            plot([0 0], [0 max(x)+5], 'k--')
            ylabel('Number of significant electrodes')
            xlabel('Time (ms)')
        hold off
        f=getframe(gca);
%% Fig 2 (NotShown)- CBPT - SntcType

v1 = erp.CBPT_sntcBinary.noProverb.data;
v2 = erp.CBPT_sntcBinary.proverb.data;

Y = nan(size(v1 ,1), size(v1,2), 2, size(v1,3));
Y(:,:,1,:) = v1;
Y(:,:,2,:) = v2;

t = erp.times.clasicos;

CONFIG.EEGBLANK         = [];
CONFIG.clusteralpha     = 0.01;     % alpha level of the sample-specific test statistic that will be used for thresholding
CONFIG.minnbchan        = 2;        % minimum number of neighborhood channels that is required for a selected
                                    % sample to be included in the clustering algorithm (default=0).
CONFIG.numrandomization = 500;      % number of draws from the permutation distribution

Nsuj = size(v1,3);
[stat times ftdata grandavg cfg] = my_ftstat(Y,t,[],1:Nsuj,CONFIG,CHANS);
 

indpos = find([stat.posclusters.prob]<cfg.alpha);
if ~isempty(indpos)
    npos = length(indpos);
    figure;
        z = stat.posclusterslabelmat;
        x = ismember(z,indpos);
        imagesc(t,1:128,x)
end

indneg = find([stat.negclusters.prob]<cfg.alpha);
if ~isempty(indneg)
    nneg = length(indneg);
    figure;
        z = stat.negclusterslabelmat;
        x = ismember(z,indneg);
        imagesc(t,1:128,x)
end












