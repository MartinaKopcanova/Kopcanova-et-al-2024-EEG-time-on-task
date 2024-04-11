%% ============================================================================ 
% power spectra plots by behaviour/time on task -- PLOTS!!
% =============================================================================
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
    power_spec = zeros(tt,371);
    block = zeros(1,tt);
%     power_spec = squeeze(mean(datastruc.powspctrm,2));
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        power_spec(h,:) = 10.^fooof_results(h).power_spectrum;
        block(1,h) = datastruc.trialinfo(h).block;
    end

    
    %aperiodic

    first_ap = squeeze(mean(aperiodic_spec(find( block == 1),:),1));
    first_ap_all(u,:) = first_ap;

    second_ap = squeeze(mean(aperiodic_spec(find( block == 2),:),1));
    second_ap_all(u,:) = second_ap;

    third_ap = squeeze(mean(aperiodic_spec(find( block == 3),:),1));
    third_ap_all(u,:) = third_ap;

    fourth_ap = squeeze(mean(aperiodic_spec(find( block == 4),:),1));
    fourth_ap_all(u,:) = fourth_ap;

    fifth_ap = squeeze(mean(aperiodic_spec(find( block == 5),:),1));
    fifth_ap_all(u,:) = fifth_ap;

    % periodic power

    first_pow = squeeze(mean(power_spec(find( block == 1),:),1));
    first_pow_all(u,:) = first_pow;

    second_pow = squeeze(mean(power_spec(find( block == 2),:),1));
    second_pow_all(u,:) = second_pow;

    third_pow = squeeze(mean(power_spec(find( block == 3),:),1));
    third_pow_all(u,:) = third_pow;

    fourth_pow = squeeze(mean(power_spec(find( block == 4),:),1));
    fourth_pow_all(u,:) = fourth_pow;

    fifth_pow = squeeze(mean(power_spec(find( block == 5),:),1));
    fifth_pow_all(u,:) = fifth_pow;

    clear power_spec
end


first_aperiodic = squeeze(mean(first_ap_all,1));
first_periodic = squeeze(mean(first_pow_all,1));
second_aperiodic = squeeze(mean(second_ap_all,1));
second_periodic = squeeze(mean(second_pow_all,1));
third_aperiodic = squeeze(mean(third_ap_all,1));
third_periodic = squeeze(mean(third_pow_all,1));
fourth_aperiodic = squeeze(mean(fourth_ap_all,1));
fourth_periodic = squeeze(mean(fourth_pow_all,1));
fifth_aperiodic = squeeze(mean(fifth_ap_all,1));
fifth_periodic = squeeze(mean(fifth_pow_all,1));


% TIME-ON-TASK
ClrMap = colormap('spring');
Clr = [42, 84, 126, 168, 210];

% this one is for when using the fooof output poewr specs and ap components
% to make the plots
figure;
plot(fooof_results(1).freqs,first_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(1),:)); hold on;
plot(fooof_results(1).freqs,second_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(2),:)); hold on;
plot(fooof_results(1).freqs,third_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(3),:));
plot(fooof_results(1).freqs,fourth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(4),:));
plot(fooof_results(1).freqs,fifth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(5),:));
hold on;
b = get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'fontsize',11)
xlabel('Frequency(Hz)','FontSize',12); ylabel('Power ({\mu}V^2)','FontSize',12);title('Pre-stimulus','FontSize',13);

%highlight the freq range of the significant cluster from CBPT f-tests
load '\experiment_data\Results\block_Ftest_n36_0.1Hz_prestim_ToT_fitlm.mat';
maskhold1 = double(statsstruc.posclusterslabelmat == 1);
maskhold2 = double(statsstruc.posclusterslabelmat == 2);
freq1 = mean(maskhold1,1);
freq2 = mean(maskhold2,1);
freqsC1 = statsstruc.freq(1,find(freq1));
freqsC2 = statsstruc.freq(1,find(freq2));
scatter(freqsC1, ones(1,length(freqsC1))*0.015,'ko','filled');
scatter(freqsC2, ones(1,length(freqsC2))*0.015,'o','filled','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
legend('Block 1','Block 2','Block 3', 'Block 4','Block 5','','','FontSize',10,'Location','southeast'); legend('box','off');
%add an inset to show the beta frequencies separately on a smaller scale
%from 14-30Hz maybe
axes('Position',[.5 .55 .35 .35]); box off;
set(gca,'FontSize',14);
plot(fooof_results(1).freqs(:,96:171),first_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(1),:)); hold on;
plot(fooof_results(1).freqs(:,96:171),second_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(2),:)); hold on;
plot(fooof_results(1).freqs(:,96:171),third_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(3),:));
plot(fooof_results(1).freqs(:,96:171),fourth_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(4),:));
plot(fooof_results(1).freqs(:,96:171),fifth_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(5),:));
xlim([12.5 20]); xlabel('Frequency(Hz)','FontSize',12); ylabel('Power ({\mu}V^2)','FontSize',12);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',11)
% % this one is if using power_spec gotten from the FFT files
% freqs = datastruc.freq;
% figure;
% plot(freqs,first_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(1),:)); hold on;
% plot(freqs,second_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(2),:)); hold on;
% plot(freqs,third_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(3),:));
% plot(freqs,fourth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(4),:));
% plot(freqs,fifth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(5),:));
% hold off;
% legend('Block 1','Block 2','Block 3', 'Block 4','Block 5'); legend('box','off');
% xlabel('Frequency(Hz)'); ylabel('Power ({\mu}V^2)');title('Pre-stimulus');

%% the same as above but this time with pre-stimulus period
% power spectra plots by behaviour/time on task
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
    power_spec = zeros(tt,371);
    block = zeros(1,tt);
%     power_spec = squeeze(mean(datastruc.powspctrm,2));
%     power_spec = power_spec(:,11:end);
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        power_spec(h,:) = 10.^fooof_results(h).power_spectrum;
        block(1,h) = datastruc.trialinfo(h).block;
    end

    
    %aperiodic

    first_ap = squeeze(mean(aperiodic_spec(find( block == 1),:),1));
    first_ap_all(u,:) = first_ap;

    second_ap = squeeze(mean(aperiodic_spec(find( block == 2),:),1));
    second_ap_all(u,:) = second_ap;

    third_ap = squeeze(mean(aperiodic_spec(find( block == 3),:),1));
    third_ap_all(u,:) = third_ap;

    fourth_ap = squeeze(mean(aperiodic_spec(find( block == 4),:),1));
    fourth_ap_all(u,:) = fourth_ap;

    fifth_ap = squeeze(mean(aperiodic_spec(find( block == 5),:),1));
    fifth_ap_all(u,:) = fifth_ap;

    % periodic power

    first_pow = squeeze(mean(power_spec(find( block == 1),:),1));
    first_pow_all(u,:) = first_pow;

    second_pow = squeeze(mean(power_spec(find( block == 2),:),1));
    second_pow_all(u,:) = second_pow;

    third_pow = squeeze(mean(power_spec(find( block == 3),:),1));
    third_pow_all(u,:) = third_pow;

    fourth_pow = squeeze(mean(power_spec(find( block == 4),:),1));
    fourth_pow_all(u,:) = fourth_pow;

    fifth_pow = squeeze(mean(power_spec(find( block == 5),:),1));
    fifth_pow_all(u,:) = fifth_pow;

    clear power spec;
end


first_aperiodic = squeeze(mean(first_ap_all,1));
first_periodic = squeeze(mean(first_pow_all,1));
second_aperiodic = squeeze(mean(second_ap_all,1));
second_periodic = squeeze(mean(second_pow_all,1));
third_aperiodic = squeeze(mean(third_ap_all,1));
third_periodic = squeeze(mean(third_pow_all,1));
fourth_aperiodic = squeeze(mean(fourth_ap_all,1));
fourth_periodic = squeeze(mean(fourth_pow_all,1));
fifth_aperiodic = squeeze(mean(fifth_ap_all,1));
fifth_periodic = squeeze(mean(fifth_pow_all,1));


% TIME-ON-TASK
ClrMap = colormap('spring');
Clr = [42, 84, 126, 168, 210];
% this one is for when using the fooof output poewr specs and ap components
% to make the plots
figure;
plot(fooof_results(1).freqs,first_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(1),:)); hold on;
plot(fooof_results(1).freqs,second_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(2),:)); hold on;
plot(fooof_results(1).freqs,third_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(3),:));
plot(fooof_results(1).freqs,fourth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(4),:));
plot(fooof_results(1).freqs,fifth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(5),:));
hold on;
b = get(gca,'XTickLabel');
set(gca,'XTickLabel',b,'fontsize',11)
xlabel('Frequency(Hz)','FontSize',12); ylabel('Power ({\mu}V^2)','FontSize',12);title('Post-stimulus','FontSize',13);

%highlight the freq range of the significant cluster from CBPT f-tests
load '\experiment_data\Results\block_Ftest_n36_0.1Hz_poststim_ToT_fitlm.mat';
maskhold1 = double(statsstruc.posclusterslabelmat == 1);
maskhold2 = double(statsstruc.posclusterslabelmat == 2);
freq1 = mean(maskhold1,1);
freq2 = mean(maskhold2,1);
freqsC1 = statsstruc.freq(1,find(freq1));
freqsC2 = statsstruc.freq(1,find(freq2));
scatter(freqsC1, ones(1,length(freqsC1))*0.01,'ko','filled');
scatter(freqsC2, ones(1,length(freqsC2))*0.01,'o','filled','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
% legend('Block 1','Block 2','Block 3', 'Block 4','Block 5','','','FontSize',10,'Location','southeast'); legend('box','off');
%add an inset to show the beta frequencies separately on a smaller scale
%from 14-30Hz maybe
axes('Position',[.5 .55 .35 .35]); box off;
set(gca,'FontSize',14);
plot(fooof_results(1).freqs(:,96:171),first_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(1),:)); hold on;
plot(fooof_results(1).freqs(:,96:171),second_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(2),:)); hold on;
plot(fooof_results(1).freqs(:,96:171),third_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(3),:));
plot(fooof_results(1).freqs(:,96:171),fourth_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(4),:));
plot(fooof_results(1).freqs(:,96:171),fifth_periodic(:,96:171),'-','LineWidth',1.5,'Color',ClrMap(Clr(5),:));
xlim([12.5 20]); xlabel('Frequency(Hz)','FontSize',12); ylabel('Power ({\mu}V^2)','FontSize',12);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',11)

% % this one is if using power_spec gotten from the FFT files
% freqs = datastruc.freq(:,11:end);
% figure;
% plot(freqs,first_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(1),:)); hold on;
% plot(freqs,second_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(2),:)); hold on;
% plot(freqs,third_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(3),:));
% plot(freqs,fourth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(4),:));
% plot(freqs,fifth_periodic,'-','LineWidth',1.5,'Color',ClrMap(Clr(5),:));
% hold off;
% legend('Block 1','Block 2','Block 3', 'Block 4','Block 5'); legend('box','off');
% xlabel('Frequency(Hz)'); ylabel('Power ({\mu}V^2)');title('Post-stimulus');

% =======================================================================================================
%                               END PLOTS 
% =======================================================================================================
%% ======================================================================================================
%                           STATS
% =======================================================================================================
clear all; clc; close all;
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};



for u = 1:length(subnums1)

    subid1 = subnums1{u};
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
    
    %get single elec averages for each block or early/late set of trials
    %(e.g., with block 1&2 = early and block 4&5 = late)

    power_spec = datastruc.powspctrm;
    for el = 1:32
        first(u,el,:) = squeeze(mean(power_spec(find(block == 1),el,:)));
        second(u,el,:) = squeeze(mean(power_spec(find(block == 2),el,:)));
        third(u,el,:) = squeeze(mean(power_spec(find(block == 3),el,:)));
        fourth(u,el,:) = squeeze(mean(power_spec(find(block == 4),el,:)));
        fifth(u,el,:) = squeeze(mean(power_spec(find(block == 5),el,:)));
        early(u,el,:) = squeeze(mean(power_spec(find(block == 1 | block == 2),el,:)));
        late(u,el,:) = squeeze(mean(power_spec(find(block == 4 | block == 5),el,:)));
    end
    clear power_spec datastruc
end

saveas = '\experiment_data\Results\FFT_ToT_poststim_0.1Hz_input_all_n36.mat';
save(saveas,'first','second','third','fourth','fifth','early','late');
%% =================================================================================================
%                       CBPT IN FIELDTRIP - ttest
% -=================================================================================================
clear all; close all; clc;
addpath('C:\Program Files\fieldtrip-20220113'); 
load('\experiment_data\EEG\157061\157061_timefreq.mat');
cfgN = [];
cfgN.method = 'template';
cfgN.template = 'biosemi32_neighb.mat';
neighbours = ft_prepare_neighbours(cfgN, timefreq);
labels = timefreq.elec.label;
clear timefreq;
load('\experiment_data\Results\FFT_ToT_prestim_0.1Hz_input_all_n36.mat');

entrydata = [early; late];
savdat = 'earlyVSlate';

    powerent.label = labels;
    powerent.dimord = 'subj_chan_freq';
    powerent.freq = 1:0.1:40;
    powerent.powspctrm = entrydata;
    % Perform t-test
    cfg = [];
    cfg.statistic = 'ft_statfun_depsamplesT';
    desig1 = [ones(36,1)', 2*ones(36,1)'];
    desig2 = [1:36,1:36];
    desig = [desig1;desig2];
%     cfg.latency = 'all';%
    cfg.frequency = [1 40];
    cfg.ivar = 1;
    cfg.uvar = 2;
    cfg.design = desig;
    cfg.keeptrials = 'no';
%     cfg.avgovertime = 'no';
    cfg.neighbours = neighbours;
    % montecarlo CB only for whole scalp analysis
    cfg.method = 'montecarlo';%'montecarloCB3D2'
    cfg.correctm = 'cluster';%'cluster';
    % cfg.clusteralpha is threshold p-value for cluster inclusion at single
    % datapoint level
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.clustertail = 0;
    % 0.1 for NEWneighbours and whole scalp analysis
    cfg.minnbchan = 1;%0.1;%0.35;
    % cfg.alpha is superfluous other than as a guide for cfg.correcttail.
    cfg.alpha = 0.05;
    cfg.tail = 0;
    cfg.correcttail = 'alpha';
    cfg.computestat = 'yes';
     cfg.computecritval = 'yes';
%     %
     cfg.computeprob = 'no'; 
    cfg.numrandomization = 2000;%
    %means
%     cfg.connectivity = testconmat;
    statsstruc = ft_freqstatistics(cfg,powerent);
    
    savename = ['\experiment_data\Results\',savdat,'_n36_0.1Hz_prestim_ToT_fitlm'];
    save(savename, 'statsstruc');
%     clear statsstruc;


%% =================================================================================================
%                       CBPT IN FIELDTRIP - ANOVA
% -=================================================================================================
clear all; close all; clc;
addpath('C:\Program Files\fieldtrip-20220113'); 
load('\experiment_data\EEG\157061\157061_timefreq.mat');
cfgN = [];
cfgN.method = 'template';
cfgN.template = 'biosemi32_neighb.mat';
neighbours = ft_prepare_neighbours(cfgN, timefreq);
labels = timefreq.elec.label;
clear timefreq;
load('\experiment_data\Results\FFT_ToT_poststim_0.1Hz_input_all_n36.mat');
entrydata = [first;second;third;fourth;fifth];
savdat = 'block_Ftest';

    powerent.label = labels;
    powerent.dimord = 'subj_chan_freq';
    powerent.freq = 1:0.1:40;
    powerent.powspctrm = entrydata;
    % Perform F-test
    cfg = [];
    cfg.statistic = 'ft_statfun_depsamplesFunivariate';
    desig1 = [ones(36,1)', 2*ones(36,1)', 3*ones(36,1)', 4*ones(36,1)', 5*ones(36,1)'];
    desig2 = [1:36, 1:36, 1:36, 1:36, 1:36];
    desig = [desig1;desig2];
%     cfg.latency = 'all';%
    cfg.frequency = [1 40];
    cfg.ivar = 1;
    cfg.uvar = 2;
    cfg.design = desig;
    cfg.keeptrials = 'no';
%     cfg.avgovertime = 'no';
    cfg.neighbours = neighbours;
    % montecarlo CB only for whole scalp analysis
    cfg.method = 'montecarlo';%'montecarloCB3D2'
    cfg.correctm = 'cluster';%'cluster';
    % cfg.clusteralpha is threshold p-value for cluster inclusion at single
    % datapoint level
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.clustertail = 1;
    % 0.1 for NEWneighbours and whole scalp analysis
    cfg.minnbchan = 1;%0.1;%0.35;
    % cfg.alpha is superfluous other than as a guide for cfg.correcttail.
    cfg.alpha = 0.05;
    cfg.tail = 1;
    cfg.correcttail = 'alpha';
    cfg.computestat = 'yes';
     cfg.computecritval = 'yes';
%     %
     cfg.computeprob = 'no'; 
    cfg.numrandomization = 2000;%
    %means
%     cfg.connectivity = testconmat;
    statsstruc = ft_freqstatistics(cfg,powerent);
    
    savename = ['\experiment_data\Results\',savdat,'_n36_0.1Hz_poststim_ToT_fitlm.mat'];
    save(savename, 'statsstruc');
%     clear statsstruc;
%% =================================================================================================
%               TOPOGRAPHIES OF MEAN F VALUES ACROSS ELECTRODES
% ==================================================================================================
clear all; close all; clc;
% addpath('C:\Users\OSS\eeglab2021.1\');
addpath('C:\Program Files\eeglab2022.1');
[ALLEEG, EEG, ~, ALLCOM] = eeglab;

allconds = {'block_Ftest_n36_0.1Hz_prestim_ToT_fitlm','block_Ftest_n36_0.1Hz_poststim_ToT_fitlm'};

for ccc = 2%:length(allconds)
    %         for t = 1:length(testtypes);
    load(['\experiment_data\Results\',allconds{ccc},'.mat']);
    comtype = allconds{ccc};
    %         testty = testtypes{t};
    sta = statsstruc.stat;%staaa(:,201:750);

    maskhold1 = double(statsstruc.posclusterslabelmat == 1);
    maskhold2 = double(statsstruc.posclusterslabelmat == 2);

                
    test1 = statsstruc.posclusterslabelmat.*maskhold1;
    
    %to get sig elec locations
    elec_1 = mean(maskhold1,2);
    els = 1:32
    els_c1 = els(find(elec_1));

    elec_2 = mean(maskhold2,2);
    els_c2 = els(find(elec_2));


    % to get an F at each elecs but only in the sig freqs:
    freq1 = mean(maskhold1,1);
    C1_f = sta(:,find(freq1));
    meanF_1 = mean(C1_f,2);

    freq2 = mean(maskhold2,1);
    C2_f = sta(:,find(freq2));
    meanF_2 = mean(C2_f,2);

    figure;
    topoplot(meanF_1,'\experiment_data\channlocs32.ced','emarker2',{[els_c1], '.','w',18},'maplimits',[0, 7]);
    % change to cbrewer colourscheme?
%     colors = cbrewer( 'RdBu', 64);
    colors = colormap('spring');
    colors = flipud(colors);
    colormap(colors)

     figure;
    topoplot(meanF_2,'\experiment_data\channlocs32.ced','emarker2',{[els_c2], '.','w',18},'maplimits',[0, 7]);
    colors = colormap('spring');
    colors = flipud(colors);
    colormap(colors)
%        
    h = colorbar;
    h.Label.String = 'F-Value';
    h.FontSize = 15;
    %                             h.Label.FontSize = 20;
    % %                             h.Position = [0.87 0.25 0.05 0.55];
    %                             h.Label.Position = [2 0 0];
    %                             h.FontSize = 16;
    %                             title('Cluster Mean','FontSize',18);

           
           
end
%% =================================================================================================
%               TOPOGRAPHIES OF MEAN power for block 1 and 5
% ==================================================================================================
clear all; close all; clc;
% addpath('C:\Users\OSS\eeglab2021.1\');
addpath('C:\Program Files\eeglab2022.1');
[ALLEEG, EEG, ~, ALLCOM] = eeglab;

allconds = {'block_Ftest_n36_0.1Hz_prestim_ToT_fitlm','block_Ftest_n36_0.1Hz_poststim_ToT_fitlm'};

for ccc = 1%:length(allconds)
    %         for t = 1:length(testtypes);
    load(['\experiment_data\Results\',allconds{ccc},'.mat']);
    comtype = allconds{ccc};
    sta = statsstruc.stat;

    maskhold1 = double(statsstruc.posclusterslabelmat == 1);
    maskhold2 = double(statsstruc.posclusterslabelmat == 2);

                
    test1 = statsstruc.posclusterslabelmat.*maskhold1;
    
    %to get sig elec locations
    elec_1 = mean(maskhold1,2);
    els = 1:32;
    els_c1 = els(find(elec_1));

    elec_2 = mean(maskhold2,2);
    els_c2 = els(find(elec_2));



    % to get an F at each elecs but only in the sig freqs:
    freq1 = mean(maskhold1,1);
    C1_f = sta(:,find(freq1));
    meanF_1 = mean(C1_f,2);

    freq2 = mean(maskhold2,1);
    C2_f = sta(:,find(freq2));
    meanF_2 = mean(C2_f,2);
    
    if ccc == 1
        load('\experiment_data\Results\FFT_ToT_prestim_0.1Hz_input_all_n36.mat');
    elseif ccc ==2
        load('\experiment_data\Results\FFT_ToT_poststim_0.1Hz_input_all_n36.mat');
    end

    b1 = squeeze(mean(first,1));
    b5 = squeeze(mean(fifth,1));

    b1_c1 = mean(b1(:,find(freq1)),2);
    b1_c2 = mean(b1(:,find(freq2)),2);

    b5_c1 = mean(b5(:,find(freq1)),2);
    b5_c2 = mean(b5(:,find(freq2)),2);
   
    figure;
    topoplot(b1_c1,'\experiment_data\channlocs32.ced','emarker2',{[els_c1], '.','w',18},'maplimits',[0, 0.3]);
    % change to cbrewer colourscheme?
%     colors = cbrewer( 'RdBu', 64);
    colors = colormap('spring');
    colors = flipud(colors);
    colormap(colors)

     figure;
    topoplot(b5_c1,'\experiment_data\channlocs32.ced','emarker2',{[els_c1], '.','w',18},'maplimits',[0, 0.3]);
    % change to cbrewer colourscheme?
%     colors = cbrewer( 'RdBu', 64);
    colors = colormap('spring');
    colors = flipud(colors);
    colormap(colors)

    h = colorbar;
    h.Label.String = '{\mu}V^2';
    h.FontSize = 10;
    h.Label.FontSize = 12;
%                             h.Position = [0.87 0.25 0.05 0.55];
    h.Label.Position = [2 0 0];
    h.FontSize = 16;
                        
    figure;
    topoplot(b1_c2,'\experiment_data\channlocs32.ced','emarker2',{[els_c2], '.','w',18},'maplimits',[0, 0.07]);
    % change to cbrewer colourscheme?
%     colors = cbrewer( 'RdBu', 64);
    colors = colormap('spring');
    colors = flipud(colors);
    colormap(colors)

     figure;
    topoplot(b5_c2,'\experiment_data\channlocs32.ced','emarker2',{[els_c2], '.','w',18},'maplimits',[0, 0.07]);
    % change to cbrewer colourscheme?
%     colors = cbrewer( 'RdBu', 64);
    colors = colormap('spring');
    colors = flipud(colors);
    colormap(colors)

    h = colorbar;
    h.Label.String = '{\mu}V^2';
    h.FontSize = 10;
    h.Label.FontSize = 12;
%                             h.Position = [0.87 0.25 0.05 0.55];
    h.Label.Position = [2 0 0];
    h.FontSize = 16;

           
           
end
