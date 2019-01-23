load coords_LNI_128_toreplaceinEEG
t = load('tiemposERP_103'); t = t.t;
imPath = ['~/Dropbox/tmp/' ];

%% imagesc
e = [];
f = []; %(igual que e, pero sin promediar entre sujetos)
fields = {'H_1','H0','H1'};

for i = 1:length(fields)
	fldName = fields{i};
	a = squeeze(ERPH.(fldName)(:,:,1,:));

	c = [];
	for iSuj = 1:max(SujsID.(fldName))
		b = mean(a(:,:,SujsID.(fldName)==iSuj),3)';
		c = cat(3, c, b);
	end

	f = cat(4,f,c);
	d = mean(c,3);
	e = cat(3,e,d);
end

clim = [-1 1];
figure(2);clf
elect = [75:85];
for i = 1:3
    subplot(3,1,i)
    tmp = e(:,elect,i)';% - e(:,elect,i+1)';
    imagesc(t,elect,tmp,clim)
    colorbar
end
name = [imPath 'imagesc1.png'];
print(name,'-dpng')


%% Por sujeto
clim = [-1 1];
figure(2);clf
elect = [70:80];
for i = 1:25
    subplot(5,5,i)
    tmp = f(:,elect,i,1)';% - e(:,elect,i+1)';
    imagesc(t,elect,tmp,clim)
end
name = [imPath 'imagescPorSuj.png'];
print(name,'-dpng')

%% Topo 
figure(2);clf
timelim = [350 500];
name = [imPath 'resta350-500.png'];
e = [];
for i = 1:length(fields)
	fldName = fields{i};
	a = squeeze(ERPH.(fldName)(:,:,1,:));

	c = [];
	for iSuj = 1:max(SujsID.(fldName))
		b = mean(a(:,:,SujsID.(fldName)==iSuj),3)';
		
		ventana = b(t>timelim(1) & t<timelim(2),:);
		meanWindow = nanmean(ventana,1);

		c = cat(1, c, meanWindow);
	end

	d = mean(c,1);
	e = cat(1,e,d);
end


lim = [-2.5 2.5];
lim = clim;
subplot(1,3,1) 
topoplot(e(2,:) - e(1,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
%set(gcf,'Color','w')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test-1.png'], 'png')

subplot(1,3,2) 
topoplot(e(2,:) - e(3,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
%set(gcf,'Color','w')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test0.png'], 'png')
title(['[' num2str(timelim) ']'])

subplot(1,3,3) 
topoplot(e(3,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
set(gcf,'Color','w')
print(name,'-dpng')




%% Topo con mask

load coords_LNI_128_toreplaceinEEG
t = load('tiemposERP_103'); t = t.t;

basepath = '/media/brunobian/DataNew/Proverbs_wxw/results/matrices/';
analisis = {'ERPs', 'Freqs','relPos012_new','relPos-101_new'};
band = {{''},{'alpha_simpleModel', 'beta', 'theta_simpleModel'},{''},{''}};

permutations = {'across','within_subjects','within_words','lm'};

elec = {electrodes};

cfg.tail = 1;
tails = {'neg', 'abs', 'pos'};
tn = tails{cfg.tail + 2};

imPath = ['./figures_paper/' ];
imPath = ['~/Dropbox/tmp/' ];

iAnalisis = 4;
iBand = 1;

alphaTval     = 2;
alphaClusters = 0.05;

im = 2; 

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


iv = 2;
    v      = fields{iv};

    figName = [imPath perm_type '_' num2str(iv) '_' v];
%     close all

    tval   = (values.t.(v)(:, :, 1));
    p      = values.p.(v)(:, :, 1);
    tval_r = (values.t.(v)(:, :, 2:end));

    thisClusters    = clusters.(v).(tn);
    clustersSign    = find(pval.(v).(tn) < alphaClusters);
    clustersFinales = ismember(thisClusters, clustersSign);

    z1 = thisClusters .* clustersFinales;

iv = 3;
    v      = fields{iv};

    figName = [imPath perm_type '_' num2str(iv) '_' v];
%     close all

    tval   = (values.t.(v)(:, :, 1));
    p      = values.p.(v)(:, :, 1);
    tval_r = (values.t.(v)(:, :, 2:end));

    thisClusters    = clusters.(v).(tn);
    clustersSign    = find(pval.(v).(tn) < alphaClusters);
    clustersFinales = ismember(thisClusters, clustersSign);

    z2 = thisClusters .* clustersFinales;

	figure();
	lim = [-4 4];
	timelim = [300 500];
	name = [imPath 'mask350-500.png'];
	ventana = tval(:,t>timelim(1) & t<timelim(2));
	meanWindow = nanmean(ventana,2);
	clusnum = 1;
	mask1 = any(z1(:,t>timelim(1) & t<timelim(2))'>=clusnum)';
	mask2 = any(z2(:,t>timelim(1) & t<timelim(2))'>=clusnum)';

	subplot(1,2,1) 
	topoplot(meanWindow .* mask1, ...
	         CHANS.chanlocs, ...
	         'maplimits',lim , ...
	         'colormap', colormap('parula'));

	topoplot(meanWindow .* mask2, ...
	         CHANS.chanlocs, ...
	         'maplimits',lim , ...
	         'colormap', colormap('parula'));

	set(gcf,'Color','w')
	print(name,'-dpng')

