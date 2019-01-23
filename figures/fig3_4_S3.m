%% Load
clear all
close all
load('../data/EEG') % load EEG pre-analized data
addpath('../functions/')

colores = [[1 1 1]; ...
           [1 .3 .3]; [0 1 0]; [0 0 1]; ...
           [1 1 0]; [0 1 1]; [1 0 1]; ...
           [.5 1 0]; [0 .5 1]; [1 .5 0]; [.5 0 1]; [1 0 .5]];

elec = {[33,51,52,53,54,62,63,64,65,66,67,68,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,97,98,99,107,108,109,100,110,111,113,114,115]};
%% Fig 3-4-S3: LMM-CBP plots
t = erp.times.clasicos;

basepath = '../LMM-CBP/csv_out/';

permutations = {'within_subjects','within_words'};

cfg.tail = 1;
tails = {'neg', 'abs', 'pos'};
tn = tails{cfg.tail + 2};

original          = 1;
butterfly         = 1;
model             = 1;
topoplt_clusters  = 1;
topoplt_original  = 1;
plot_numberSig    = 1;
plot_Model_numberSig = 1;

alphaTval     = 2;
alphaClusters = 0.05;

for im = 1:length(permutations);

perm_type = permutations{im};
load([basepath '/clusters_' perm_type '_' tn '_' ])

load([basepath '/pvals_' perm_type '_' tn '_' ])

load([basepath '/Original_' perm_type '_' ])

fields = fieldnames(clusters);
fprintf('Permutation= %s\n', perm_type)

for iv = 1:length(fields) 
    v      = fields{iv};

    tval   = (values.t.(v)(:, :, 1));
    p      = values.p.(v)(:, :, 1);
    tval_r = (values.t.(v)(:, :, 2:end));

    thisClusters    = clusters.(v).(tn);
    clustersSign    = find(pval.(v).(tn) < alphaClusters);
    clustersFinals = ismember(thisClusters, clustersSign);

    z = thisClusters .* clustersFinals;

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
                plot([0 0], [1 max(elect)], 'w--')
                set(gca, 'YTickLabel',{})
                set(gca, 'XTickLabel',{})
                box on
            hold off
            title(v)
            f=getframe(gca);
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
            
    end

    if model
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
    end
    
    if plot_Model_numberSig
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
    end
end
end
disp('Fin')
