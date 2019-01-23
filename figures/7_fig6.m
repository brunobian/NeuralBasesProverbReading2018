load coords_LNI_128_toreplaceinEEG
t = load('tiemposERP_103'); t = t.t;

electrodes = [33,51,52,53,54,62,63,64,65,66,67,68,76,77,78,79,80,82,83,84,85,86,87,89,90,91,92,93,97,98,99,107,108,109,100,110,111,113,114,115];
electrodes = [33,51,52,53,54,62,63,64,65,66,76,77,78,79,80,81,82,83,84,85,86,87,97,98,99,100,109,110,114,115];


fields = {'H_1','H0','H1'};

roi1 = [2,3,4,34,112];
roi2 = [85,86,87,75,88];
roi3 = [82,83,84,78,91];

imPath = ['~/Dropbox/tmp/' ];
imName = 'test_roi3.png';
ROI=roi3;
e = [];
for i = 1:length(fields)
	fldName = fields{i};
	a = squeeze(mean(ERPH.(fldName)(ROI,:,1,:),1));

	c = [];
	for iSuj = 1:max(SujsID.(fldName))
		b = mean(a(:,SujsID.(fldName)==iSuj),2)';
		c = cat(1, c, b);
	end

	d = mean(c,1);
	e = cat(1,e,d);
end

% ERPs
figure();clf
hold on
plot(t,e(1,:),'LineWidth',2)
plot(t,e(2,:),'LineWidth',2)
plot(t,e(3,:),'LineWidth',2)
plot([0 0],[-3 3],'k-','LineWidth',2)
plot([-100 700],[0 0],'k-','LineWidth',2)
xlim([-100 700])
legend({'-1','0','1'})
f=getframe(gca);
imwrite(f.cdata,[imPath imName], 'png')
hold off

%% Topoplots
load coords_LNI_128_toreplaceinEEG
fields = {'H_1','H0','H1'};
timelim = [350 450];

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

% e = cat(1, e, e(2,:) - e(1,:));
% e = cat(1, e, e(3,:) - e(2,:));

for i = 1:(length(fields)-1)
	fldName1 = fields{i};
	fldName2 = fields{i+1};

	a1 = squeeze(ERPH.(fldName1)(:,:,1,:));
	a2 = squeeze(ERPH.(fldName2)(:,:,1,:));

	c = [];
	for iSuj = 1:max(SujsID.(fldName))
		b1 = mean(a1(:,:,SujsID.(fldName1)==iSuj),3)';
		b2 = mean(a2(:,:,SujsID.(fldName2)==iSuj),3)';
		
		ventana = b1(t>timelim(1) & t<timelim(2),:) - b2(t>timelim(1) & t<timelim(2),:);
		meanWindow = nanmean(ventana,1);

		c = cat(1, c, meanWindow);
	end

	d = mean(c,1);
	e = cat(1,e,d);
end

lim = [-2.5 2.5];
figure(2);clf
subplot(2,3,1) 
topoplot(e(1,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
%set(gcf,'Color','w')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test-1.png'], 'png')

subplot(2,3,2) 
topoplot(e(2,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
%set(gcf,'Color','w')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test0.png'], 'png')

subplot(2,3,3) 
topoplot(e(3,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
%set(gcf,'Color','w')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test1.png'], 'png')
colorbar

lim = [-1.5 1.5];
subplot(2,3,4) 
topoplot(e(4,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
%set(gcf,'Color','w')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test1.png'], 'png')

subplot(2,3,5) 
topoplot(e(5,:), ...
         CHANS.chanlocs, ...
         'maplimits',lim , ...
         'colormap', colormap('jet'));
     colorbar
set(gcf,'Color','w')
print('~/Dropbox/tmp/testFull.png', '-dpng')
%f=getframe(gca);
%imwrite(f.cdata,[imPath 'test1.png'], 'png')
