%% Fig  - Results - ERPs clásicos - Pred*Type original
load('../data/EEG.mat')
corro_stats = 0;
topo_todos = 0;
t = erp.times.clasicos;

% ERP
close all
x_lim=[-100 700];
y_lim=[-4 4];

figure(1);clf
set(gcf,'Color','w')
    plot([0 0],y_lim,'k-'); 
    hold on
    plot(x_lim,[0 0],'k-');

    tmp = erp.pred_type_ROIN400;

    [~, h1] = niceBars2(t, tmp.LoProv.m, tmp.LoProv.e, [1 0 0], .1,'-');
    [~, h2] = niceBars2(t, tmp.HiProv.m, tmp.HiProv.e, [1 0 0], .1,'--');
    [~, h3] = niceBars2(t, tmp.LoCommon.m, tmp.LoCommon.e, [0 0 1], .1,'-');
    [~, h4] = niceBars2(t, tmp.HiCommon.m, tmp.HiCommon.e, [0 0 1], .1,'--');

    
    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')

    
%% Fig  - Results - ERPs clásicos - Pred*Type remef
% 2/8/18
% Para esta figura:
% Calculamos el ERP para cada trial promediando todos los electrodos N400
% Fitteamos un modelo igual al que usamos en el analisis completo
% Removimos el efecto estimado de posición
% Traemos los datos
% Volvemos a hacer la figura anterior

EEG.all = csvread('../data/EEG_remef.csv');
columnas = {'id','suj_id','n_orac','palnum','tipo','pred','freq','length','bad_epoch','stopword','lngth','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14','t15','t16','t17','t18','t19','t20','t21','t22','t23','t24','t25','t26','t27','t28','t29','t30','t31','t32','t33','t34','t35','t36','t37','t38','t39','t40','t41','t42','t43','t44','t45','t46','t47','t48','t49','t50','t51','t52','t53','t54','t55','t56','t57','t58','t59','t60','t61','t62','t63','t64','t65','t66','t67','t68','t69','t70','t71','t72','t73','t74','t75','t76','t77','t78','t79','t80','t81','t82','t83','t84','t85','t86','t87','t88','t89','t90','t91','t92','t93','t94','t95','t96','t97','t98','t99','t100','t101','t102','t103'};
data = {'id','suj_id','n_orac','palnum','tipo','pred','freq','length','bad_epoch','stopword','lngth'};

%%
x_lim=[-100 700];
y_lim=[-4 4];

indProv = EEG.all(:,5)==0; %tipo proverb
indCommon = EEG.all(:,5)==1; %tipo no-proverb

pred = EEG.all(:,6); % pred
lims = quantile(pred, [0,.33,.66,1]);

predProv = EEG.all(indProv, 6); % pred solo de proverbios
sujProv  = EEG.all(indProv, 2); 
EEG.prov = EEG.all(indProv, length(data)+1 : end)';

predCom = EEG.all(indCommon, 6); % pred solo de no proverbios 
sujCom  = EEG.all(indCommon, 2); 
EEG.com = EEG.all(indCommon, length(data)+1 : end)';

q1p = predProv >= lims(1) & predProv < lims(2);
q3p = predProv >= lims(3) & predProv < lims(4);
q1c = predCom >= lims(1) & predCom < lims(2);
q3c = predCom >= lims(3) & predCom < lims(4);

for i = 1:max(unique(sujCom))
    ERP(i).prov_q1 = nanmean(EEG.prov(:,sujProv == i & q1p),2);
    ERP(i).prov_q3 = nanmean(EEG.prov(:,sujProv == i & q3p),2);

    ERP(i).com_q1  = nanmean(EEG.com(:, sujCom == i & q1c),2);
    ERP(i).com_q3  = nanmean(EEG.com(:, sujCom == i & q3c),2);
end

erp.prov_q1 = nanmean([ERP.prov_q1],2)';
erp.prov_q3 = nanmean([ERP.prov_q3],2)';
erp.prov_e1 = nanstd([ERP.prov_q1],0,2)/sqrt(size([ERP.prov_q1],2));
erp.prov_e3 = nanstd([ERP.prov_q3],0,2)/sqrt(size([ERP.prov_q3],2));

erp.com_q1  = nanmean([ERP.com_q1],2)';
erp.com_q3  = nanmean([ERP.com_q3],2)';
erp.com_e1  = nanstd([ERP.com_q1],0,2)/sqrt(size([ERP.com_q1],2));
erp.com_e3  = nanstd([ERP.com_q3],0,2)/sqrt(size([ERP.com_q3],2));

figure(1);clf
set(gcf,'Color','w')
hold on
    plot([0 0], y_lim, 'k-')
    plot(x_lim, [0 0], 'k-')

    [~, h1] = niceBars2(t, erp.prov_q1, erp.prov_e1', [1 0 0], .1,'-');
    [~, h2] = niceBars2(t, erp.prov_q3, erp.prov_e3', [1 0 0], .1,'--');
    [~, h3] = niceBars2(t, erp.com_q1,  erp.com_e1',  [0 0 1], .1,'-');
    [~, h4] = niceBars2(t, erp.com_q3,  erp.com_e3',  [0 0 1], .1,'--');

    legend([h1,h2,h3,h4],...
        'LowPred - Mem Rel', 'HiPred - Mem Rel', ...
        'LowPred - Common',  'HiPred - Common')
    
    xlim(x_lim) 
    xlabel('Time [ms]')
    ylim(y_lim)
    ylabel('ERP amplitude [microVolts]')
    box on

hold off






