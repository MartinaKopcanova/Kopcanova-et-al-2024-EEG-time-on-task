%% ========================================================================
%  ==================    FIGURE 6 + regression (not in paper) =============
%  ========================================================================


%% ===================================================================================================================
%                                                    Exponent purely by
%                                                    block (no beh control)
%=====================================================================================================================
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
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_preStim_0.1Hz_ALL.mat']);
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_preStim_powerspec_0.1Hzres.mat']);
    datastruc = powerspec;
    clear powerspec;

    tt = length(datastruc.trialinfo);
    difficulty = zeros(1,tt);
    confidence = zeros(1,tt);
    accuracy = zeros(1,tt);
    trialorder = zeros(1,tt);
    rts = zeros(1,tt);
    block = zeros(1,tt);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        block(1,h) = datastruc.trialinfo(h).block;
    end

    exB1(u,:) = mean(exponent(:,find(block == 1)));
    exB2(u,:) = mean(exponent(:,find(block == 2)));
    exB3(u,:) = mean(exponent(:,find(block == 3)));
    exB4(u,:) = mean(exponent(:,find(block == 4)));
    exB5(u,:) = mean(exponent(:,find(block == 5)));
end

%% boxplots or else for each block (to show variability)ClrMap = colormap('spring');
ClrMap = colormap('cool');
Clr = [42, 84, 126, 168, 210];

figure;
params = [exB1;exB2;exB3;exB4;exB5];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean Exponent','FontSize',12); hold on;
xlabel('Block','FontSize',12); title('Pre-stimulus','FontSize',13);
ylim([0.3 1.9]);
% y = [1,2];
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)

%stats (ANOVA?)
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',exB1,exB2,exB3,exB4,exB5,'VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
c = multcompare(rm,'Blocks')
%% =============================================================================================================
%                                   end
% ==============================================================================================================

%% ===================================================================================================================
%                                                    Exponent purely by
%                                                    block (no beh control)
%=====================================================================================================================
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
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_postStim_0.1Hz_ALL.mat']);
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_postStim_powerspec_0.1Hzres.mat']);
    datastruc = powerspec;
    clear powerspec;

    tt = length(datastruc.trialinfo);
    difficulty = zeros(1,tt);
    confidence = zeros(1,tt);
    accuracy = zeros(1,tt);
    trialorder = zeros(1,tt);
    rts = zeros(1,tt);
    block = zeros(1,tt);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        block(1,h) = datastruc.trialinfo(h).block;
    end

    exB1(u,:) = mean(exponent(:,find(block == 1)));
    exB2(u,:) = mean(exponent(:,find(block == 2)));
    exB3(u,:) = mean(exponent(:,find(block == 3)));
    exB4(u,:) = mean(exponent(:,find(block == 4)));
    exB5(u,:) = mean(exponent(:,find(block == 5)));
end

%% boxplots or else for each block (to show variability)ClrMap = colormap('spring');
ClrMap = colormap('cool');
Clr = [42, 84, 126, 168, 210];

figure;
params = [exB1;exB2;exB3;exB4;exB5];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1); ylim([0.3 1.9]);
ylabel('Mean Exponent','FontSize',12); hold on;
xlabel('Block','FontSize',12); title('Post-stimulus','FontSize',13);
% y = [1,2];
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 
a = get(gca,'YTickLabel'); 
set(gca,'YTickLabel',a,'fontsize',11)

%stats (ANOVA?)
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',exB1,exB2,exB3,exB4,exB5,'VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
c = multcompare(rm,'Blocks')
%% =============================================================================================================
%                                   end
% ==============================================================================================================
%% regression
% clear all; clc; close all;
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

for u = 1:length(subnums1)

    tic
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_preStim_0.1Hz_ALL.mat']);
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
   
    joffset = squeeze(offset);
    offset_rank = tiedrank(joffset);
    jexponent = squeeze(exponent);
    exponent_rank = tiedrank(jexponent);

    diff_rank = tiedrank(difficulty);
    trialorder_rank = tiedrank(trialorder);
    rt_rank = tiedrank(rts);

    % create data table
    Regr_table = table(diff_rank',offset_rank',exponent_rank',confidence',accuracy',...
        trialorder_rank',rt_rank','VariableNames',{'diff_rank','offset_rank','exponent_rank','confidence','accuracy','trialorder_rank','rt_rank'});

    %             CONFIDENCE MODEL CONTROLLING FOR DIFFICULTY AND ACCURACY
    modelspec = 'exponent_rank ~ confidence + diff_rank + accuracy + rt_rank + trialorder_rank';
    mdl = fitlm(Regr_table,modelspec);

    get_conf = mdl.Coefficients(3,1);
    conf_ex(u) = table2array(get_conf);
    
    get_diff = mdl.Coefficients(2,1);
    diff_ex(u) = table2array(get_diff);

    get_acc = mdl.Coefficients(4,1);
    acc_ex(u) = table2array(get_acc);

    get_rt = mdl.Coefficients(6,1);
    rt_ex(u) = table2array(get_rt);

    get_trial = mdl.Coefficients(5,1);
    trial_ex(u) = table2array(get_trial);



%     modelspec2 = 'exponent_rank ~ confidence + diff_rank + accuracy + rt_rank';
%     mdl2 = fitlm(Regr_table,modelspec2);
% 
%     get_conf_ex = mdl2.Coefficients(3,1);
%     conf_EEG_ex = table2array(get_conf_ex);
%     
%     get_diff_ex = mdl2.Coefficients(2,1);
%     diff_EEG_ex = table2array(get_diff_ex);
% 
%     get_acc_ex = mdl2.Coefficients(4,1);
%     acc_EEG_ex = table2array(get_acc_ex);
% 
%     get_rt_ex = mdl2.Coefficients(5,1);
%     rt_EEG_ex = table2array(get_rt_ex);


    toc;
    %     %Save beta matrices
   
end
saveas = '\experiment_data\Results\AperiodicRegression_prestim__exponent_N36.mat';
    save(saveas, 'conf_ex','trial_ex','rt_ex','acc_ex','diff_ex');%
    clear datastruc conf_ex trial_ex rt_ex acc_ex diff_ex; %
%%
clear all; clc; close all;
load('\experiment_data\Results\AperiodicRegression_poststim__exponent_N36.mat');
test = zeros(36,1);
entrydata = {'conf_ex','acc_ex','diff_ex','rt_ex','trial_ex'};
for entd = 1:length(entrydata)
    [h,p,ci,stat] = ttest(eval(entrydata{entd})',test);
    h_all(entd) = h;
    p_all(entd) = p;
    ci_all(entd,:) = ci;
    stat_all(entd) = stat.tstat;
end