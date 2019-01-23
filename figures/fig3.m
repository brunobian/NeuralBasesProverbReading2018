%% Fig 2 and 3: STATS - LMM (permutation) - plot
load coords_LNI_128_toreplaceinEEG
t = load('~/Dropbox/Labo Juan/gerardo_juan_diego/oraciones120/mfiles/tiemposERP_103'); t = t.t;

% basepath = '/media/brunobian/ExtraDrive1/Proverb_WxW/matrices/';
basepath = '/media/brunobian/DataNew/Proverbs_wxw/results/matrices/';
analisis = {'ERPs', 'Freqs','relPos012_new','relPos-101_new'};
band = {{''},{'alpha_simpleModel', 'beta', 'theta_simpleModel'},{''},{''}};

permutations = {'across','within_subjects','within_words','lm'};
electrodes = [1,2,3,4,5,6,19,32,33,34,35,50,51,52,53,54,64,65,66,67,68,75,76,77,84,85,86,87,88,89,90,97,98,99,100,109,110,111,112,113,114,115,124];
electrodes = [33,51,52,53,54,62,63,64,65,66,76,77,78,79,80,81,82,83,84,85,86,87,97,98,99,100,109,110,114,115];
electrodes = [33,51,52,53,54,62,63,64,65,66,67,68,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,97,98,99,107,108,109,100,110,111,113,114,115];

elec = {electrodes};

cfg.tail = 1;
tails = {'neg', 'abs', 'pos'};
tn = tails{cfg.tail + 2};

original          = 0;
butterfly         = 0;
modelo            = 0;
topoplt_clusters  = 0;
topoplt_original  = 0;
plot_numberSig    = 0;
plot_Modelo_numberSig = 0;
guardar           = 0;
pausar            = 0;  
imPath = ['./figures_paper/originales/lmm/' ];
% imPath = ['~/Dropbox/tmp/' ];

colores = [[1 1 1]; ...
           [1 .3 .3]; [0 1 0]; [0 0 1]; ...
           [1 1 0]; [0 1 1]; [1 0 1]; ...
           [.5 1 0]; [0 .5 1]; [1 .5 0]; [.5 0 1]; [1 0 .5]];

iAnalisis = 1;
iBand = 1;

alphaTval     = 2;
alphaClusters = 0.05;


for im = 3 %1:length(permutations)%:5;

perm_type = permutations{im};
load([basepath analisis{iAnalisis} ...
      '/clusters_' perm_type '_' ...
      tn '_' band{iAnalisis}{iBand}])

load([basepath analisis{iAnalisis} ...
      '/pvals_' perm_type '_' ...
      tn '_' band{iAnalisis}{iBand}])

load([basepath analisis{iAnalisis} ...
      '/Original_' perm_type '_' ...
      band{iAnalisis}{iBand}])

fields = fieldnames(clusters);
fprintf('Permutation= %s\n', perm_type)

for iv = [3, 4,6 ]%:length(fields) 
    v      = fields{iv};

    figName = [imPath perm_type '_' num2str(iv) '_' v];
%     close all

    tval   = (values.t.(v)(:, :, 1));
    p      = values.p.(v)(:, :, 1);
    tval_r = (values.t.(v)(:, :, 2:end));

    thisClusters    = clusters.(v).(tn);
    clustersSign    = find(pval.(v).(tn) < alphaClusters);
    clustersFinales = ismember(thisClusters, clustersSign);

    z = thisClusters .* clustersFinales;

    fprintf('Term = %s\n', v)
    fprintf('Sign = %.2f\n', pval.(v).(tn)(clustersSign)')
        

    % plots matrices
    if original
        figure();
        set(gcf,'Color','w')
            col_ax = [-8 8];
            mask = abs(tval)>alphaTval;
            elect = [1:128];
            tmp = tval.*mask;
            imagesc(t, elect, tmp(elect,:), col_ax);
            hold on
                plot([0 0], [1 elect], 'w--')
                set(gca, 'YTickLabel',{})
                set(gca, 'XTickLabel',{})
                box on
            hold off
            title(v)
            f=getframe(gca);
            if guardar; imwrite(f.cdata,[figName '_orig.png'], 'png'); end
    end
    
    if butterfly
        x_lim=[-100 700];
        y_lim=[-10 10];

        figure();
        set(gcf,'Color','w')
            col_ax = [-8 8];
            elect = [1:128];
            hold on
            for i = elect
                plot(t,tval(elect,:),'color',[.8 .8 .8])
            end
            plot(t,mean(tval,1),'color','k', 'LineWidth', 2)
            plot([0 0],y_lim,'k-'); 
            plot(x_lim,[0 0],'k-');
            box on
            hold off
            xlim(x_lim) 
            ylim(y_lim)
%             xlabel('Time [ms]')
%             ylabel('t-value')
%             
            if guardar;saveas(gcf,[figName '_butterfly.pdf'],'pdf');
; end
    end

    if modelo
        figure();
        set(gcf,'Color','w')
            elect = [1:128];
            imagesc(t, elect, z(elect,:),[0 length(unique(z))])
            hold on
                colormap(colores(1:(length(unique(z))+1),:))
                plot([0 0], [1.1 127.90], 'k--')
                set(gca, 'YTickLabel',{})
                set(gca, 'XTickLabel',{})
                box on
            hold off
            f=getframe(gca);
            if guardar; imwrite(f.cdata,[figName '_clusters.png'], 'png'); end
    end
    
    % Topoplots
    lim = [-4 4];
    timelim = [300 500];
    ventana = tval(:,t>timelim(1) & t<timelim(2));
    meanWindow = nanmean(ventana,2);
    if topoplt_clusters
        figure();
        clusnum = 1;
        mask = any(z(:,t>timelim(1) & t<timelim(2))'>=clusnum)';
        
        topoplot(meanWindow .* mask, ...
                 CHANS.chanlocs, ...
                 'maplimits',lim , ...
                 'colormap', colormap('parula'));
        set(gcf,'Color','w')
        
        f=getframe(gca);
        if guardar; saveas(gcf,[figName '_topo_cluster.eps'],'eps2c');end
        if guardar; saveas(gcf,[figName '_topo_cluster.png'],'png');end

    end
    
    if topoplt_original
        figure();
        topoplot(meanWindow, ...
                 CHANS.chanlocs, ...
                 'emarker2', elec, ...
                 'maplimits',lim , ...
                 'colormap', colormap('parula'));
        set(gcf,'Color','w')
        f=getframe(gca);
        if guardar; imwrite(f.cdata,[figName '_topo_orig.png'], 'png'); end
    end
   
    if plot_numberSig
        figure();
        set(gcf,'Color','w')
        x = sum(z,1);
        plot(t,x, 'LineWidth', 2)
        hold on
            xlim([-105.5 700])
            ylim([0 max(x)+5])
            plot([0 0], [0 max(x)+5], 'k--')
            ylabel('Number of significant electrodes')
            xlabel('Time (ms)')
        hold off
        f=getframe(gca);
        if guardar; saveas(gcf,[figName '_numberSig.eps'],'eps2c');end
    end
    
    if plot_Modelo_numberSig
        figure();
        set(gcf,'Color','w')
            elect = [1:128];
            imagesc(t, elect, z(elect,:),[0 length(unique(z))])
            hold on
                plot([0 0], [1.1 127.90], 'k--')
                yyaxis left
                colormap(colores(1:(length(unique(z))+1),:))
                set(gca, 'YTickLabel',{})
                ylabel('Electrode Number')
                yyaxis right
                x = sum(z,1);
                plot(t,x, 'LineWidth', 3)
                xlabel('Time (ms)')
                ylabel('Number of significant electrodes')
                box on
        hold off
        f=getframe(gca);
        if guardar; saveas(gcf,[figName '_numberSigModel.eps'],'eps2c');end
    end
    

    
    if pausar
        pause
    end
end
end
% close all
disp('Fin')

%% (Intercept)
% addpath(genpath('/home/juank/toolbox/eeglab14_1_1b/'),'-end')
iv = 1
    v      = fields{iv};
    fprintf('%s\n', v)

%     figName = [imPath perm_type '_' num2str(iv) '_' v];
%     close all

    tval   = (values.t.(v)(:, :, 1));
    p      = values.p.(v)(:, :, 1);
    tval_r = (values.t.(v)(:, :, 2:end));

    thisClusters    = clusters.(v).(tn);
    clustersSign    = find(pval.(v).(tn) < alphaClusters);
    clustersFinales = ismember(thisClusters, clustersSign);

    z = thisClusters .* clustersFinales;

    figure(1);clf
        set(gcf,'Color','w')
        timelimlist = {[60 110],[150 225],[275 425]};
        for i=1:length(timelimlist)
            subplot(2,3,i)
                topoplot(nanmean(tval(:,t>timelimlist{i}(1) & t<timelimlist{i}(2)),2), ...
                         CHANS.chanlocs, ...
                         'maplimits',lim , ...
                         'colormap', colormap('parula')); %, ...
        end
        subplot(2,1,2)
            elect = [1:128];
            mask = abs(tval)>alphaTval;
            tmp = tval.*mask;
            plot(t, tval(elect,:),'Color',[.7 .7 .7]);
            hold on
                plot(t, std(tval(elect,:)),'k');
                plot(t, -std(tval(elect,:)),'k');
                YLIMI = ylim;
                for i=1:length(timelimlist)
%                     plot([timelimlist{i}(1) timelimlist{i}(1)],ylim, 'b-')
%                     plot([timelimlist{i}(2) timelimlist{i}(2)],ylim, 'b-')
                    patch([timelimlist{i}(1) timelimlist{i}(1) timelimlist{i}(2) timelimlist{i}(2) timelimlist{i}(1)],...
                        [YLIMI YLIMI(2) YLIMI(1) YLIMI(1)], 'b','FaceAlpha',.25)
                end
                plot([min(t) max(t)],[0 0], 'k--')
                plot([0 0], ylim, 'k--')
                set(gca,'XLim',[min(t) max(t)])
                box on
            hold off
            title(v)
%         f=getframe(gca);
%         if guardar; imwrite(f.cdata,[figName '_orig.png'], 'png'); end
    
%% Histograma de repeticion de palabras
% Para justificar por qué el analisis de palabras da diferente al analisis
% de sujetos, mirar esto


a = arrayfun(@(x) x.palabras', DATA, 'UniformOutput', 0);
filt = [DATA.palnum] > 2 & [DATA.palnum] < 8 & ...
        ismember([DATA.catsimple], 'anv') &...
        [DATA.length] > 2;

b = [a{:}];
[C,ia,ic]=unique([b(filt)]);

length(C)
hist(hist(ic,1:257))

