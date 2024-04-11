%% ==========================================================================================================================================
% FIGURE 8 & STATS
% ===========================================================================================================================================
clear all; clc; close all;
cd('\experiment_data\EEG');
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

for u = 1:length(subnums1)
    tic
    subid1 = subnums1{u};
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_PeakAlphaFreq_fooof_0.1Hz_preStim.mat']);
    %for trial info load: 
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_preStim_powerspec_0.1Hzres.mat']);
    datastruc = powerspec;
    clear powerspec;


    tt = length(datastruc.trialinfo);
    difficulty = zeros(1,tt);
    confidence = zeros(1,tt);
    accuracy = zeros(1,tt);
    trialorder = zeros(1,tt);
    rts = zeros(1,tt);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
    end
   
    %remove trials with no alpha found, from behaviour too!
    [~,~,CF] = find(peakCF);
    [~,~,PW] = find(maxAlp);
    confidence = confidence(find(peakCF));
    difficulty = difficulty(find(peakCF));
    accuracy = accuracy(find(peakCF));
    trialorder = trialorder(find(peakCF));
    rts = rts(find(peakCF));

    % everything rank ordered apart from accuracy and confidence!!!!;
    CF_rank = tiedrank(CF);
    PW_rank = tiedrank(PW);
    diff_rank = tiedrank(difficulty);
    trialorder_rank = tiedrank(trialorder);
    rt_rank = tiedrank(rts);

    % create data table
    Regr_table = table(diff_rank',PW_rank',CF_rank',confidence',accuracy',...
        trialorder_rank',rt_rank','VariableNames',{'diff_rank','PW_rank','CF_rank','confidence','accuracy','trialorder_rank','rt_rank'});

    %             CONFIDENCE MODEL CONTROLLING FOR DIFFICULTY AND ACCURACY
    modelspec = 'CF_rank ~ confidence + diff_rank + accuracy + rt_rank + trialorder_rank';
    mdl = fitlm(Regr_table,modelspec);

%     get_conf = mdl.Coefficients(3,1);
%     conf_CF(u) = table2array(get_conf);
%     
%     get_diff = mdl.Coefficients(2,1);
%     diff_CF(u) = table2array(get_diff);
% 
%     get_acc = mdl.Coefficients(4,1);
%     acc_CF(u) = table2array(get_acc);
% 
%     get_rt = mdl.Coefficients(6,1);
%     rt_CF(u) = table2array(get_rt);
% 
%     get_trial = mdl.Coefficients(5,1);
%     trial_CF(u) = table2array(get_trial);

    aic_CF(u) = mdl.ModelCriterion.AIC;
    bic_CF(u) = mdl.ModelCriterion.BIC;


    modelspec2 = 'PW_rank ~ confidence + diff_rank + accuracy + rt_rank + trialorder_rank';
    mdl2 = fitlm(Regr_table,modelspec2);

%     get_conf2 = mdl2.Coefficients(3,1);
%     conf_PW(u) = table2array(get_conf2);
%     
%     get_diff2 = mdl2.Coefficients(2,1);
%     diff_PW(u) = table2array(get_diff2);
% 
%     get_acc2 = mdl2.Coefficients(4,1);
%     acc_PW(u) = table2array(get_acc2);
% 
%     get_rt2 = mdl2.Coefficients(6,1);
%     rt_PW(u) = table2array(get_rt2);
% 
%     get_trial2 = mdl2.Coefficients(5,1);
%     trial_PW(u) = table2array(get_trial2);

    aic_PW(u) = mdl2.ModelCriterion.AIC;
    bic_PW(u) = mdl2.ModelCriterion.BIC;
      
    
    clear peakCF maxAlp CF_rank PW_rank Regr_table diff_rank rt_rank trialorder_rank CF PW;
    toc;

end

saveas = '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_prestim.mat';
    save(saveas, 'aic_CF', 'bic_CF','aic_PW', 'bic_PW');
% %     %Save beta matrices
%     saveas = '\experiment_data\Results\Regression_alphaCF_PW_poststim.mat';
%     save(saveas, 'conf_CF', 'rt_CF','diff_CF','acc_CF','trial_CF','conf_PW', 'rt_PW','diff_PW','acc_PW','trial_PW');%
%     clear datastruc conf_CF rt_CF diff_CF acc_CF trial_CF conf_PW rt_PW diff_PW acc_PW trial_PW; %


%% ==========================================================================================================================================
% REGRESSIONS NOW without time-on-task
% ===========================================================================================================================================
clear all; clc; close all;
cd('\experiment_data\EEG'); %for work laptop 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

for u = 1:length(subnums1)
    tic
    subid1 = subnums1{u};
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_PeakAlphaFreq_fooof_0.1Hz_preStim.mat']);
    %for trial info load: 
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_preStim_powerspec_0.1Hzres.mat']);
    datastruc = powerspec;
    clear powerspec;


    tt = length(datastruc.trialinfo);
    difficulty = zeros(1,tt);
    confidence = zeros(1,tt);
    accuracy = zeros(1,tt);
    trialorder = zeros(1,tt);
    rts = zeros(1,tt);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
    end
   
    %remove trials with no alpha found, from behaviour too!
    [~,~,CF] = find(peakCF);
    [~,~,PW] = find(maxAlp);
    confidence = confidence(find(peakCF));
    difficulty = difficulty(find(peakCF));
    accuracy = accuracy(find(peakCF));
    trialorder = trialorder(find(peakCF));
    rts = rts(find(peakCF));

    % everything rank ordered apart from accuracy and confidence!!!!;
    CF_rank = tiedrank(CF);
    PW_rank = tiedrank(PW);
    diff_rank = tiedrank(difficulty);
    trialorder_rank = tiedrank(trialorder);
    rt_rank = tiedrank(rts);

    % create data table
    Regr_table = table(diff_rank',PW_rank',CF_rank',confidence',accuracy',...
        trialorder_rank',rt_rank','VariableNames',{'diff_rank','PW_rank','CF_rank','confidence','accuracy','trialorder_rank','rt_rank'});

    %             CONFIDENCE MODEL CONTROLLING FOR DIFFICULTY AND ACCURACY
    modelspec = 'CF_rank ~ confidence + diff_rank + accuracy + rt_rank';
    mdl = fitlm(Regr_table,modelspec);

%     get_conf = mdl.Coefficients(3,1);
%     conf_CF(u) = table2array(get_conf);
%     
%     get_diff = mdl.Coefficients(2,1);
%     diff_CF(u) = table2array(get_diff);
% 
%     get_acc = mdl.Coefficients(4,1);
%     acc_CF(u) = table2array(get_acc);
% 
%     get_rt = mdl.Coefficients(6,1);
%     rt_CF(u) = table2array(get_rt);
% 
%     get_trial = mdl.Coefficients(5,1);
%     trial_CF(u) = table2array(get_trial);

    aic_CF2(u) = mdl.ModelCriterion.AIC;
    bic_CF2(u) = mdl.ModelCriterion.BIC;


    modelspec2 = 'PW_rank ~ confidence + diff_rank + accuracy + rt_rank';
    mdl2 = fitlm(Regr_table,modelspec2);

%     get_conf2 = mdl2.Coefficients(3,1);
%     conf_PW(u) = table2array(get_conf2);
%     
%     get_diff2 = mdl2.Coefficients(2,1);
%     diff_PW(u) = table2array(get_diff2);
% 
%     get_acc2 = mdl2.Coefficients(4,1);
%     acc_PW(u) = table2array(get_acc2);
% 
%     get_rt2 = mdl2.Coefficients(6,1);
%     rt_PW(u) = table2array(get_rt2);
% 
%     get_trial2 = mdl2.Coefficients(5,1);
%     trial_PW(u) = table2array(get_trial2);

    aic_PW2(u) = mdl2.ModelCriterion.AIC;
    bic_PW2(u) = mdl2.ModelCriterion.BIC;
      
    
    clear peakCF maxAlp CF_rank PW_rank Regr_table diff_rank rt_rank trialorder_rank CF PW;
    toc;

end

saveas = '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_notot_prestim.mat';
    save(saveas, 'aic_CF2', 'bic_CF2','aic_PW2', 'bic_PW2');
%% poststimulus
clear all;
load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_notot_poststim.mat';
load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_poststim.mat';


n = 1:36;
aic_diff = aic_CF2 - aic_CF;
figure;
scatter(n,aic_diff,30,'o','filled');hold on;
yline(0,'k-'); ylabel('{\Delta} AIC (without t-o-t - with t-o-t)'); xlabel('Participant N^o');
title('Poststimulus \alpha center frequency');

[h,p,ci,stat] = ttest(aic_CF2,aic_CF)

aic_diff_PW = aic_PW2 - aic_PW;
figure;
scatter(n,aic_diff_PW,30,'o','filled');hold on;
yline(0,'k-'); ylabel('{\Delta} AIC (without t-o-t - with t-o-t)'); xlabel('Participant N^o');
title('Poststimulus \alpha peak power');

[h,p,ci,stat] = ttest(aic_PW2,aic_PW)

figure;
 x = ones(1,36); x2 = 2.*ones(1,36); X = [x, x2]; %0.6350 0.0780 0.1840 for dark red
 b = boxplot([aic_diff,aic_diff_PW],X,'Colors',[0 0 0; 0 0 0],'Symbol',''); xticklabels({'Freq','Power'});
 set(b,'LineWidth',1.1); set(gca,'FontSize',14);
 ylabel('{\Delta} AIC','FontSize',14); hold on;
 y = [1,2];
 scatter(1,aic_diff,'k','filled','MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); hold on;
 scatter(2,aic_diff_PW,"filled",'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); 
 hold off;

%% prestimulus now
%%
clear all;
load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_notot_prestim.mat';
load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_prestim.mat';


n = 1:36;
aic_diff = aic_CF2 - aic_CF;
figure;
scatter(n,aic_diff,30,'o','filled');hold on;
yline(0,'k-'); ylabel('{\Delta} AIC (without t-o-t - with t-o-t)'); xlabel('Participant N^o');
title('Prestimulus \alpha center frequency'); hold off;

[h,p,ci,stat] = ttest(aic_CF2,aic_CF)

aic_diff_PW = aic_PW2 - aic_PW;
figure;
scatter(n,aic_diff_PW,30,'o','filled');hold on;
yline(0,'k-'); ylabel('{\Delta} AIC (without t-o-t - with t-o-t)'); xlabel('Participant N^o');
title('Prestimulus \alpha peak power'); hold off;

[h,p,ci,stat] = ttest(aic_PW2,aic_PW)

 figure;
 x = ones(1,36); x2 = 2.*ones(1,36); X = [x, x2]; %0.6350 0.0780 0.1840 for dark red
 b = boxplot([aic_diff,aic_diff_PW],X,'Colors',[0 0 0; 0 0 0],'Symbol',''); xticklabels({'Freq','Power'});
 set(b,'LineWidth',1.1); set(gca,'FontSize',14);
 ylabel('{\Delta} AIC','FontSize',14); hold on;
 y = [1,2];
 scatter(1,aic_diff,'k','filled','MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); hold on;
 scatter(2,aic_diff_PW,"filled",'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); 
 hold off;

%% =============================================================================================
%%==============================================================================================
%===============================================================================================
%                           Figure with regress t-vals and AIC
% ==============================================================================================

%% STREAMLINED VERSION

clear all; clc; close all;
load '\experiment_data\Results\tvals_Regression_alphaCF_PW_prestim.mat';
totpre = stat_all; clear stat_all;
load '\experiment_data\Results\tvals_Regression_alphaCF_PW_poststim.mat';
totpost = stat_all; clear stat_all;
load '\experiment_data\Results\tvals_Regression_noToT_alphaCF_PW_prestim.mat';
notot_pre = stat_all; clear stat_all;
load '\experiment_data\Results\tvals_Regression_noToT_alphaCF_PW_poststim.mat';
notot_post = stat_all; clear stat_all;

% the order of the tvals is as follows:
% for no-tot entrydata = {'conf_CF', 'rt_CF','diff_CF','acc_CF','conf_PW', 'rt_PW','diff_PW','acc_PW'};
% for tot: entrydata = {'conf_CF', 'rt_CF','diff_CF','acc_CF','trial_CF','conf_PW', 'rt_PW','diff_PW','acc_PW','trial_PW'};

subplot(2,3,1)
% figure;
y = [notot_pre(2) totpre(2); notot_pre(1) totpre(1); 0 totpre(5)]; %0 totpre(5)
b = bar(y); b(2).FaceColor = "#0C359E"; %[0.0078431  0.3098  0.12941];
b(1).FaceColor = "#EE99C2"; %[ 0.05098 0.0078431 0.3098]; 
set(gca, 'XTick',[1 2 3],'XTickLabel',{'RT','Confidence','ToT'},'FontSize',14);
title('Peak \alpha frequency'); ylabel('T-value');
ylim([-8 7]); legend('Without ToT included','With ToT included','Location','southwest'); legend('Box','off');

subplot(2,3,2)
% figure;
y = [notot_pre(6) totpre(7); notot_pre(5) totpre(6); 0 totpre(10)]; %0 totpre(10)
b = bar(y); b(2).FaceColor = "#0C359E";%[0.0078431  0.3098  0.12941];
b(1).FaceColor = "#EE99C2";%[ 0.05098 0.0078431 0.3098]; 
set(gca, 'XTick',[1 2 3],'XTickLabel',{'RT','Confidence','ToT'},'FontSize',14);
title('Peak \alpha power'); ylabel('T-value');
ylim([-8 7]);

load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_notot_prestim.mat';
load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_prestim.mat';

n = 1:36;
aic_diff = aic_CF2 - aic_CF;
aic_diff_PW = aic_PW2 - aic_PW;

subplot(2,3,3);
x = ones(1,36); x2 = 2.*ones(1,36); X = [x, x2]; %0.6350 0.0780 0.1840 for dark red
b = boxplot([aic_diff,aic_diff_PW],X,'Colors',[0 0 0; 0 0 0],'Symbol',''); xticklabels({'Freq','Power'});
set(b,'LineWidth',1.1); set(gca,'FontSize',14);
ylabel('{\Delta} AIC','FontSize',14); hold on;
y = [1,2];
scatter(1,aic_diff,'k','filled','MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); hold on;
scatter(2,aic_diff_PW,"filled",'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); 
hold off;

subplot(2,3,4)
% figure;
y = [notot_post(2) totpost(2); notot_post(1) totpost(1); 0 totpost(5)]; %0 totpre(5)
b = bar(y); b(2).FaceColor = "#0C359E";%[0.0078431  0.3098  0.12941];
b(1).FaceColor = "#EE99C2";%[ 0.05098 0.0078431 0.3098]; 
set(gca, 'XTick',[1 2 3],'XTickLabel',{'RT','Confidence','ToT'},'FontSize',14);
title('Peak \alpha frequency'); ylabel('T-value');
ylim([-8 7]);

subplot(2,3,5)
% figure;
y = [ notot_post(6) totpost(7);notot_post(5) totpost(6); 0 totpost(10)]; %0 totpre(10)
b = bar(y); b(2).FaceColor = "#0C359E";%[0.0078431  0.3098  0.12941];
b(1).FaceColor = "#EE99C2";%[ 0.05098 0.0078431 0.3098]; 
set(gca, 'XTick',[1 2 3],'XTickLabel',{'RT','Confidence','ToT'},'FontSize',14);
title('Peak \alpha power'); ylabel('T-value');
ylim([-8 7]);

load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_notot_poststim.mat';
load '\experiment_data\Results\ModelComp_Regression_alphaCF_PW_poststim.mat';


n = 1:36;
aic_diff = aic_CF2 - aic_CF;
aic_diff_PW = aic_PW2 - aic_PW;
subplot(2,3,6);
x = ones(1,36); x2 = 2.*ones(1,36); X = [x, x2]; %0.6350 0.0780 0.1840 for dark red
b = boxplot([aic_diff,aic_diff_PW],X,'Colors',[0 0 0; 0 0 0],'Symbol',''); xticklabels({'Freq','Power'});
set(b,'LineWidth',1.1); set(gca,'FontSize',14);
ylabel('{\Delta} AIC','FontSize',14); hold on;
y = [1,2];
scatter(1,aic_diff,'k','filled','MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); hold on;
scatter(2,aic_diff_PW,"filled",'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.4,'jitter','on','jitterAmount',0.1); 
hold off;

