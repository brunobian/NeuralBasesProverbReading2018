%% Fig 2 (a-c)     - Results - ERPs clásicos - Pred
corro_stats = 1;
topo_todos = 0;
save_path = './figures_paper/Round1/originales/';

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
     
    [a, h1] = niceBars2(t, erp.pred_ROIN400.q1.m, erp.pred_ROIN400.q1.e, [1 0 0], .3,'-');
    [a, h2] = niceBars2(t, erp.pred_ROIN400.q2.m, erp.pred_ROIN400.q2.e, [0 .5 0], .3,'-');
    [a, h3] = niceBars2(t, erp.pred_ROIN400.q3.m, erp.pred_ROIN400.q3.e, [0 0 1], .3,'-');

    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    legend([h1, h2, h3], 'Low Pred','Mid Pred','High Pred')
    saveas(gcf,[save_path 'ERP_pred.pdf'],'pdf');

    % Topoplot resta hi-lo N400
    times         = win.N400.time;
    limites_color = [-3 3];
    elec_marks    = win.N400.elect;    
    t_ind         = (erp.times.clasicos >= times(1)) & ...
                    (erp.times.clasicos <  times(2));

    amps =  erp.pred.q1.m(:,t_ind) - erp.pred.q3.m(:,t_ind);
    figure();clf;set(gcf,'Color','w')
        topoplot(mean(amps,2), ...
                 CHANS.chanlocs,'maplimits',limites_color);
%          colorbar
     saveas(gcf,[save_path 'topo_pred.pdf'],'pdf');
     saveas(gcf,[save_path 'topo_pred.png'],'png');
        
    figure();clf;set(gcf,'Color','w')
        topoplot(zeros(size(amps,1),1), ...
                 CHANS.chanlocs,'maplimits',limites_color, ...
                 'emarker2', {elec_marks,'o','k',10,1});
    saveas(gcf,[save_path 'elect_pred.pdf'],'pdf');
            
if topo_todos 
    % Topoplot N400
    x = -1.3;
    y = 0;
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

        saveas(gcf,[save_path 'topoplotClassic.eps'],'eps2c');
end

% Stats
if corro_stats
    
    tmp_q1 = []; tmp_q2 = []; tmp_q3 = [];

    for iE = 1:length(ERP)
        tmp_q1 = [tmp_q1 ERP(iE).pred_ROI_timeN400.q1.m];
        tmp_q2 = [tmp_q2 ERP(iE).pred_ROI_timeN400.q2.m];
        tmp_q3 = [tmp_q3 ERP(iE).pred_ROI_timeN400.q3.m];
    end
    
    M = [tmp_q1' tmp_q2' tmp_q3'];
    [p,tbl,stats] = kruskalwallis(M);

    % Plot 
    m_400(1) = nanmean(tmp_q1); e_400(1) = nanstd(tmp_q1) / sqrt(length(tmp_q1));
    m_400(2) = nanmean(tmp_q2); e_400(2) = nanstd(tmp_q2) / sqrt(length(tmp_q2));
    m_400(3) = nanmean(tmp_q3); e_400(3) = nanstd(tmp_q3) / sqrt(length(tmp_q3));

    figure(); clf
    colores = {[1 0 0], [0 .5 0], [0 0 1]} ;
    set(gcf,'color','w');
        hold on
        for i=1:3
            errorbar(i, m_400(i), e_400(i),'o','Color',colores{i},'LineWidth',3)
        end
        set(gca, 'Xtick', [1:3], 'xTickLabel', {'Low', 'Mid', 'Hi'})
        xlabel('Predictability tercile')
        ylabel('Mean amplitude in window')
        set(gcf, 'Units', 'Inches', 'Position', [1, 2, 3, 4.7])
        set(gcf,'PaperPositionMode','auto')
        print([save_path 'meanWindows.eps'],'-depsc','-tiff','-r0')
end
%% Fig 2 (b-d)     - Results - ERPs clásicos - SntcType
corro_stats = 1;

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
    
    [a, h1] = niceBars2(t, erp.sntcBinary_ROIN400.proverb.m, erp.sntcBinary_ROIN400.proverb.e, [1 0 0], .3,'-');
    [a, h2] = niceBars2(t, erp.sntcBinary_ROIN400.noProverb.m, erp.sntcBinary_ROIN400.noProverb.e, [0 .5 0], .3,'-');
    
    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    legend([h1, h2], 'Proverb','noProverb')
    saveas(gcf,[save_path 'ERP_SntcType.pdf'],'pdf');


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
    saveas(gcf,[save_path 'topo_SntcType0.png'],'png');

figure();clf
set(gcf,'Color','w')
    topoplot(mean(erp.sntcBinary.noProverb.m(:,t_ind),2), ...
                 CHANS.chanlocs,'maplimits',limites_color);
    saveas(gcf,[save_path 'topo_SntcType.pdf'],'pdf');
    saveas(gcf,[save_path 'topo_SntcType1.png'],'png');

% Stats
if corro_stats
    tmp_proverb = []; tmp_NoProv = []; 

    for iE = 1:length(ERP)
        tmp_proverb = [tmp_proverb ERP(iE).sntcBinary_ROIN400.proverb.m];
        tmp_NoProv  = [tmp_NoProv  ERP(iE).sntcBinary_ROIN400.noProverb.m];
    end

    M = [tmp_proverb', tmp_NoProv'];
    [p1,h,stats1] = ranksum(tmp_proverb', tmp_NoProv')
    [p,tbl,stats] = kruskalwallis(M);


    % Plot
    m_type(1) = nanmean(tmp_proverb); e_type(1) = nanstd(tmp_proverb) / sqrt(length(tmp_proverb));
    m_type(2) = nanmean(tmp_NoProv);  e_type(2) = nanstd(tmp_NoProv) / sqrt(length(tmp_NoProv));

    figure(); clf
    colores = {[1 0 0], [0 .5 0]} ;
    set(gcf,'color','w');
        hold on
        for i=1:2
            errorbar(i, m_type(i), e_type(i),'o','Color',colores{i},'LineWidth',3)
        end
        set(gca, 'Xtick', [1:2], 'xTickLabel', {'Proverb', 'Common'})
        xlabel('Predictability tercile')
        ylabel('Mean amplitude in window')
        set(gcf, 'Units', 'Inches', 'Position', [1, 2, 3, 4.7])
        set(gcf,'PaperPositionMode','auto')
        print([save_path 'meanWindowsSntcType.eps'],'-depsc','-tiff','-r0')

end
%% Fig 2 (NotShown)- Results - ERPs clásicos - Freq
corro_stats = 1;
topo_todos = 0;
save_path = './figures_paper/Round1/originales/';

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
     
    [a, h1] = niceBars2(t, erp.freq_ROIN400.q1.m, erp.freq_ROIN400.q1.e, [1 0 0], .3,'-');
    [a, h2] = niceBars2(t, erp.freq_ROIN400.q2.m, erp.freq_ROIN400.q2.e, [0 .5 0], .3,'-');
    [a, h3] = niceBars2(t, erp.freq_ROIN400.q3.m, erp.freq_ROIN400.q3.e, [0 0 1], .3,'-');

    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    legend([h1, h2, h3], 'Low freq','Mid freq','High freq')
    saveas(gcf,[save_path 'ERP_freq.pdf'],'pdf');

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
     saveas(gcf,[save_path 'topo_freq.pdf'],'pdf');
     saveas(gcf,[save_path 'topo_freq.png'],'png');
        
    figure();clf;set(gcf,'Color','w')
        topoplot(zeros(size(amps,1),1), ...
                 CHANS.chanlocs,'maplimits',limites_color, ...
                 'emarker2', {elec_marks,'o','k',10,1});
    saveas(gcf,[save_path 'elect_freq.pdf'],'pdf');
            
if topo_todos 
    % Topoplot N400
    x = -1.3;
    y = 0;
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

        saveas(gcf,[save_path 'topoplotClassic.eps'],'eps2c');
end

% Stats
if corro_stats

    tmp_q1 = []; tmp_q2 = []; tmp_q3 = [];

    for iE = 1:length(ERP)
        tmp_q1 = [tmp_q1 ERP(iE).freq_ROI_timeN400.q1.m];
        tmp_q2 = [tmp_q2 ERP(iE).freq_ROI_timeN400.q2.m];
        tmp_q3 = [tmp_q3 ERP(iE).freq_ROI_timeN400.q3.m];
    end
    
    M = [tmp_q1' tmp_q2' tmp_q3'];
    [p,tbl,stats] = kruskalwallis(M);

    % Plot
    m_400(1) = nanmean(tmp_q1); e_400(1) = nanstd(tmp_q1) / sqrt(length(tmp_q1));
    m_400(2) = nanmean(tmp_q2); e_400(2) = nanstd(tmp_q2) / sqrt(length(tmp_q2));
    m_400(3) = nanmean(tmp_q3); e_400(3) = nanstd(tmp_q3) / sqrt(length(tmp_q3));

    figure(); clf
    colores = {[1 0 0], [0 .5 0], [0 0 1]} ;
    set(gcf,'color','w');
        hold on
        for i=1:3
            errorbar(i, m_400(i), e_400(i),'o','Color',colores{i},'LineWidth',3)
        end
        set(gca, 'Xtick', [1:3], 'xTickLabel', {'Low', 'Mid', 'Hi'})
        xlabel('freq tercile')
        ylabel('Mean amplitude in window')
        set(gcf, 'Units', 'Inches', 'Position', [1, 2, 3, 4.7])
        set(gcf,'PaperPositionMode','auto')
        print([save_path 'meanWindows.eps'],'-depsc','-tiff','-r0')
end
%% Fig 2 (sup)     - Regression   
close all

Amp  = erp.regresion_pred.m;
Suj  = categorical(erp.regresion_pred.suj);
Typ1  = categorical(erp.regresion_pred.SntcType);
Typ0  = categorical(erp.regresion_pred.SntcType==0);
Pred = erp.regresion_pred.pred;

tbl = table(Amp, Pred, Suj,Typ1, Typ0 , 'VariableNames', {'Amp', 'Pred', 'Suj','Typ1','Typ0'});
lm = fitlm(tbl, 'Amp ~ Pred:Typ1 + Pred:Typ0')
p = lm.Coefficients.Estimate;
Amp2 = p(2) .* Pred + p(1);

Amp  = erp.regresion_pred.m;
Suj  = dummyvar(erp.regresion_pred.suj);
Pred = erp.regresion_pred.pred;
b    = regress(Amp, [Pred, Suj]);

Amp2 = b(2) .* Pred + b(1);

%save('tablePred','tbl')

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
saveas(gcf,[save_path 'regression.pdf'],'pdf');
%% Fig 2 (sup)     - CBPT - Pred
% Primero hay que pasar al formato fieldtrip, eso lo hacia my_ftstat
addpath('./CBPT/')
addpath('/home/brunobian/Documents/Repos/fieldtrip')
save_path = './figures_paper/Round1/originales/';

v1 = erp.CBPT_pred.q1.data;
v2 = erp.CBPT_pred.q3.data;

ft_defaults;

% v1: [128x103x25 single] 
%   Y: (channels, times, conds, subjects)
Y = nan(size(v1 ,1), size(v1,2), 2, size(v1,3));
Y(:,:,1,:) = v1;
Y(:,:,2,:) = v2;

%   t: times (ms)
t = erp.times.clasicos;

CONFIG.EEGBLANK         = [];
CONFIG.clusteralpha     = 0.01;     % alpha level of the sample-specific test statistic that will be used for thresholding
CONFIG.minnbchan        = 2;        % minimum number of neighborhood channels that is required for a selected
                                    % sample to be included in the clustering algorithm (default=0).
CONFIG.numrandomization = 500;      % number of draws from the permutation distribution

Nsuj = size(v1,3);
[stat times ftdata grandavg cfg] = my_ftstat(Y,t,[],1:Nsuj,CONFIG);
%   urevent: ...
%   Nsuj: # of subjects
%   CONFIG: Params of the original FieldTrip function 

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
saveas(gcf,[save_path 'CBPTclasico_pred.pdf'],'pdf');
colormap([1,1,1;1,0,0])
saveas(gcf,[save_path 'CBPTclasico_pred_red.pdf'],'pdf');

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
        saveas(gcf,[save_path 'CBPTclassic_topo.png'],'png')
        
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
        saveas(gcf,[save_path 'CBPTclassic_numberSig.eps'],'eps2c')
%% Fig 2 (sup)     - CBPT - SntcType

v1 = erp.CBPT_sntcBinary.noProverb.data;
v2 = erp.CBPT_sntcBinary.proverb.data;

ft_defaults;

% v1: [128x103x25 single] 
%   Y: (channels, times, conds, subjects)
Y = nan(size(v1 ,1), size(v1,2), 2, size(v1,3));
Y(:,:,1,:) = v1;
Y(:,:,2,:) = v2;

%   t: times (ms)
t = erp.times.clasicos;

CONFIG.EEGBLANK         = [];
CONFIG.clusteralpha     = 0.01;     % alpha level of the sample-specific test statistic that will be used for thresholding
CONFIG.minnbchan        = 2;        % minimum number of neighborhood channels that is required for a selected
                                    % sample to be included in the clustering algorithm (default=0).
CONFIG.numrandomization = 500;      % number of draws from the permutation distribution

Nsuj = size(v1,3);
[stat times ftdata grandavg cfg] = my_ftstat(Y,t,[],1:Nsuj,CONFIG);
%   urevent: ...
%   Nsuj: # of subjects
%   CONFIG: Params of the original FieldTrip function 

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












