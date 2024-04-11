%% =======================================================================================================================
%               PLOTS FOR ToT EFFECTS ON PEAK POWER AND FREQ
% ========================================================================================================================
clear all; clc; close all;
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

for u = 1:length(subnums1)
    subid1 = subnums1{u};
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_PeakAlphaFreq_fooof_0.1Hz_preStim.mat']);
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_preStim_0.1Hz_ALL.mat']);
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
    block = zeros(1,tt);
    flatspec = zeros(tt,371);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        block(1,h) = datastruc.trialinfo(h).block;
        flatspec(h,:) = fooof_results(h).fooofed_spectrum - fooof_results(h).ap_fit;
    end
    %alpha 8-13 = 51:101
    freqs = fooof_results(1).freqs;
    alpha = 51:101;
    flatalpha = flatspec(:,alpha);
    
    for h = 1:tt
        if ~peakCF(1,h) == 0 %if an element has a peak in alpha
            flatalphadetected(h,:) = flatalpha(h,:);
        else
            flatalphadetected(h,:) = nan(1,51);
        end
    end

    a1(u,:) = mean(flatalphadetected(find(block == 1),:),1,'omitnan');
    a2(u,:) = mean(flatalphadetected(find(block == 2),:),1,'omitnan');
    a3(u,:) = mean(flatalphadetected(find(block == 3),:),1,'omitnan');
    a4(u,:) = mean(flatalphadetected(find(block == 4),:),1,'omitnan');
    a5(u,:) = mean(flatalphadetected(find(block == 5),:),1,'omitnan');

    [~,~,fb1] = find(peakCF(1,find(block == 1)));
    [~,~,fb2] = find(peakCF(1,find(block == 2)));
    [~,~,fb3] = find(peakCF(1,find(block == 3)));
    [~,~,fb4] = find(peakCF(1,find(block == 4)));
    [~,~,fb5] = find(peakCF(1,find(block == 5)));

    [~,~,pb1] = find(maxAlp(1,find(block == 1)));
    [~,~,pb2] = find(maxAlp(1,find(block == 2)));
    [~,~,pb3] = find(maxAlp(1,find(block == 3)));
    [~,~,pb4] = find(maxAlp(1,find(block == 4)));
    [~,~,pb5] = find(maxAlp(1,find(block == 5)));
    
    Fb1(u,:) = mean(fb1);
    Fb2(u,:) = mean(fb2);
    Fb3(u,:) = mean(fb3);
    Fb4(u,:) = mean(fb4);
    Fb5(u,:) = mean(fb5);

    Pb1(u,:) = mean(pb1);
    Pb2(u,:) = mean(pb2);
    Pb3(u,:) = mean(pb3);
    Pb4(u,:) = mean(pb4);
    Pb5(u,:) = mean(pb5);
    
    
    clear peakCF maxAlp datastruc fooof_results
end
%% zoomed in spectral plot on alpha + frequency highlighted

ClrMap = colormap('spring');
Clr = [42, 84, 126, 168, 210];
figure;
plot(freqs(1,alpha), mean(a1,1),'Color',ClrMap(Clr(1),:),'LineWidth',1.5); hold on;
plot(freqs(1,alpha), mean(a5,1),'Color',ClrMap(Clr(5),:),'LineWidth',1.5);
plot(freqs(1,alpha), mean(a3,1),'Color',ClrMap(Clr(3),:),'LineWidth',1.5);
plot(freqs(1,alpha), mean(a4,1),'Color',ClrMap(Clr(4),:),'LineWidth',1.5);
plot(freqs(1,alpha), mean(a2,1),'Color',ClrMap(Clr(2),:),'LineWidth',1.5);
xlim([8 13]); ylim([0.2 0.75]);

alpha_freqs = freqs(1,alpha);
[m i] = max(mean(a1,1));
freqb1 = alpha_freqs(i);
[m i] = max(mean(a2,1));
freqb2 = alpha_freqs(i);
[m i] = max(mean(a3,1));
freqb3 = alpha_freqs(i);
[m i] = max(mean(a4,1));
freqb4 = alpha_freqs(i);
[m i] = max(mean(a5,1));
freqb5 = alpha_freqs(i);

xline(freqb1,'--','Color',ClrMap(Clr(1),:),'LineWidth',1.5);
xline(freqb2,'--','Color',ClrMap(Clr(2),:),'LineWidth',1.5);
xline(freqb3,'--','Color',ClrMap(Clr(3),:),'LineWidth',1.5);
xline(freqb4,'--','Color',ClrMap(Clr(4),:),'LineWidth',1.5);
xline(freqb5,'--','Color',ClrMap(Clr(5),:),'LineWidth',1.5);
xlim([8 13]);
ylabel('log(Periodic Power)','FontSize',12);
xlabel('Frequency(Hz)','FontSize',12);
title('Pre-stimulus','FontSize',13);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)
%% boxplots or else for each block (to show variability)
f = figure;
f.Position = [521,254,259,321];
ClrMap = colormap('spring');
Clr = [42, 84, 126, 168, 210];
params = [Fb1;Fb2;Fb3;Fb4;Fb5];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean {\alpha} frequency','FontSize',12); hold on;
xlabel('Block','FontSize',12);title('Center frequency','FontSize',13);
ylim([8.4 13.5]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)

%stats (ANOVA?)
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',Fb1,Fb2,Fb3,Fb4,Fb5,'VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
eta = ranovatbl.SumSq(1)./(ranovatbl.SumSq(1)+ranovatbl.SumSq(2))
c = multcompare(rm,'Blocks')


f = figure;
f.Position = [521,254,259,321];
params = [Pb1;Pb2;Pb3;Pb4;Pb5];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);xlabel('Block','FontSize',12);title('Peak power','FontSize',13)
ylabel('Mean {\alpha} power','FontSize',12); hold on;
ylim([0.5 2]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)

%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',Pb1,Pb2,Pb3,Pb4,Pb5,'VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
eta = ranovatbl.SumSq(1)./(ranovatbl.SumSq(1)+ranovatbl.SumSq(2))
c = multcompare(rm,'Blocks')
%% post stimulus
clear all; clc; close all;
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};


for u = 1:length(subnums1)
    subid1 = subnums1{u};
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_PeakAlphaFreq_fooof_0.1Hz_postStim.mat']);
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_postStim_0.1Hz_ALL.mat']);
    %for trial info load: 
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
    flatspec = zeros(tt,371);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        block(1,h) = datastruc.trialinfo(h).block;
        flatspec(h,:) = fooof_results(h).fooofed_spectrum - fooof_results(h).ap_fit;
    end
    %alpha 8-13 = 51:101
    freqs = fooof_results(1).freqs;
    alpha = 51:101;
    flatalpha = flatspec(:,alpha);
    
    for h = 1:tt
        if ~peakCF(1,h) == 0 %if an element has a peak in alpha
            flatalphadetected(h,:) = flatalpha(h,:);
        else
            flatalphadetected(h,:) = nan(1,51);
        end
    end

    a1(u,:) = mean(flatalphadetected(find(block == 1),:),1,'omitnan');
    a2(u,:) = mean(flatalphadetected(find(block == 2),:),1,'omitnan');
    a3(u,:) = mean(flatalphadetected(find(block == 3),:),1,'omitnan');
    a4(u,:) = mean(flatalphadetected(find(block == 4),:),1,'omitnan');
    a5(u,:) = mean(flatalphadetected(find(block == 5),:),1,'omitnan');

    [~,~,fb1] = find(peakCF(1,find(block == 1)));
    [~,~,fb2] = find(peakCF(1,find(block == 2)));
    [~,~,fb3] = find(peakCF(1,find(block == 3)));
    [~,~,fb4] = find(peakCF(1,find(block == 4)));
    [~,~,fb5] = find(peakCF(1,find(block == 5)));

    [~,~,pb1] = find(maxAlp(1,find(block == 1)));
    [~,~,pb2] = find(maxAlp(1,find(block == 2)));
    [~,~,pb3] = find(maxAlp(1,find(block == 3)));
    [~,~,pb4] = find(maxAlp(1,find(block == 4)));
    [~,~,pb5] = find(maxAlp(1,find(block == 5)));
    
    Fb1(u,:) = mean(fb1);
    Fb2(u,:) = mean(fb2);
    Fb3(u,:) = mean(fb3);
    Fb4(u,:) = mean(fb4);
    Fb5(u,:) = mean(fb5);

    Pb1(u,:) = mean(pb1);
    Pb2(u,:) = mean(pb2);
    Pb3(u,:) = mean(pb3);
    Pb4(u,:) = mean(pb4);
    Pb5(u,:) = mean(pb5);
    
    
    clear peakCF maxAlp datastruc fooof_results
end
%%
ClrMap = colormap('spring');
Clr = [42, 84, 126, 168, 210];
figure;
plot(freqs(1,alpha), mean(a1,1),'Color',ClrMap(Clr(1),:),'LineWidth',1.5); hold on;
plot(freqs(1,alpha), mean(a5,1),'Color',ClrMap(Clr(5),:),'LineWidth',1.5);
plot(freqs(1,alpha), mean(a3,1),'Color',ClrMap(Clr(3),:),'LineWidth',1.5);
plot(freqs(1,alpha), mean(a4,1),'Color',ClrMap(Clr(4),:),'LineWidth',1.5);
plot(freqs(1,alpha), mean(a2,1),'Color',ClrMap(Clr(2),:),'LineWidth',1.5);
xlim([8 13]); ylim([0.2 0.75]);

alpha_freqs = freqs(1,alpha);
[m i] = max(mean(a1,1));
freqb1 = alpha_freqs(i);
[m i] = max(mean(a2,1));
freqb2 = alpha_freqs(i);
[m i] = max(mean(a3,1));
freqb3 = alpha_freqs(i);
[m i] = max(mean(a4,1));
freqb4 = alpha_freqs(i);
[m i] = max(mean(a5,1));
freqb5 = alpha_freqs(i);

xline(freqb1,'--','Color',ClrMap(Clr(1),:),'LineWidth',1.5);
xline(freqb2,'--','Color',ClrMap(Clr(2),:),'LineWidth',1.5);
xline(freqb3,'--','Color',ClrMap(Clr(3),:),'LineWidth',1.5);
xline(freqb4,'--','Color',ClrMap(Clr(4),:),'LineWidth',1.5);
xline(freqb5,'--','Color',ClrMap(Clr(5),:),'LineWidth',1.5);
xlim([8 13]);
ylabel('log(Periodic Power)','FontSize',12);
xlabel('Frequency(Hz)','FontSize',12);
title('Post-stimulus','FontSize',13);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)
%% boxplots or else for each block (to show variability)
ClrMap = colormap('spring');
Clr = [42, 84, 126, 168, 210];
f = figure;
f.Position = [521,254,259,321];
params = [Fb1;Fb2;Fb3;Fb4;Fb5];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean {\alpha} frequency','FontSize',12); hold on;
xlabel('Block','FontSize',12);title('Center frequency','FontSize',13);
ylim([8.4 13.5]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)

%stats (ANOVA?)
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',Fb1,Fb2,Fb3,Fb4,Fb5,'VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
eta = ranovatbl.SumSq(1)./(ranovatbl.SumSq(1)+ranovatbl.SumSq(2))
c = multcompare(rm,'Blocks')


f = figure;
f.Position = [521,254,259,321];
params = [Pb1;Pb2;Pb3;Pb4;Pb5];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);xlabel('Block','FontSize',12);title('Peak power','FontSize',13)
ylabel('Mean {\alpha} power','FontSize',12); hold on;
ylim([0.5 2]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',11)

%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',Pb1,Pb2,Pb3,Pb4,Pb5,'VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
eta = ranovatbl.SumSq(1)./(ranovatbl.SumSq(1)+ranovatbl.SumSq(2))
c = multcompare(rm,'Blocks')
