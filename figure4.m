%% ========================================================================
%  ============  FIGURE 4   ===============================================
%  ========================================================================

clear all; clc; close all;
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

load '\experiment_data\Results\block_Ftest_n36_0.1Hz_poststim_ToT_fitlm.mat';
    maskhold1 = double(statsstruc.posclusterslabelmat == 1);
    maskhold2 = double(statsstruc.posclusterslabelmat == 2);
    freq1 = mean(maskhold1,1);
    freq2 = mean(maskhold2,1);
    freqsC1 = statsstruc.freq(1,find(freq1));
    freqsC2 = statsstruc.freq(1,find(freq2));
    els1 = find(mean(maskhold1,2));
    els2 = find(mean(maskhold2,2));

for u = 1:length(subnums1)

    tic
    subid1 = subnums1{u};
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_postStim_powerspec_0.1Hzres.mat']);
    datastruc = powerspec;
    clear powerspec;
    
    alpha_freq_start = find(datastruc.freq == round(freqsC1(1,1),1));
    alpha_freq_end = find(datastruc.freq == round(freqsC1(1,end),1));

    beta_freq_start = find(datastruc.freq == round(freqsC2(1,1),1));
    beta_freq_end = find(datastruc.freq == round(freqsC2(1,end),1));

    tt = length(datastruc.trialinfo);
    difficulty = zeros(1,tt);
    confidence = zeros(1,tt);
    accuracy = zeros(1,tt);
    trialorder = zeros(1,tt);
    rts = zeros(1,tt);
    block = zeros(1,tt);
  
    power_spec = datastruc.powspctrm;
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        block(1,h) = datastruc.trialinfo(h).block;
    end

    blength = [length(find(block == 1)), length(find(block == 2)), length(find(block == 3)), length(find(block == 4)), length(find(block == 5))];
    bl_div = blength./15;
    bl_div_round = round(bl_div);
    bl_div_round_all(u,:) = bl_div_round;


    alpha_pows = squeeze(mean(squeeze(mean(power_spec(:,els1,[alpha_freq_start : alpha_freq_end]),2)),2));
    beta_pows = squeeze(mean(squeeze(mean(power_spec(:,els2,[beta_freq_start : beta_freq_end]),2)),2));

    zalpha_pows = zscore(alpha_pows);
    zbeta_pows = zscore(beta_pows);
    
    for blc = 1:5
    bstart = 1; 
    bend = [];
        for bin = 1:14
            bstart = [bstart;bl_div_round(blc)*bin + 1];
        end
        binstart(:,blc) = bstart;
        for bin = 1:15
            bend = [bend;bl_div_round(blc)*bin];
        end
        binend(:,blc) = bend;
    end

    for blc = 1:5
        for bin = 1:15
            zalpha1 = zalpha_pows(find(block == blc),:);
            zbeta1 = zbeta_pows(find(block == blc),:);
            if bin == 15
               meanpow_alpha(blc,bin) = mean(zalpha1(binstart(bin,blc):end,:),1);
               meanpow_beta(blc,bin) = mean(zbeta1(binstart(bin,blc):end,:),1);
            else
               meanpow_alpha(blc,bin) = mean(zalpha1(binstart(bin,blc):binend(bin,blc),:),1);
               meanpow_beta(blc,bin) = mean(zbeta1(binstart(bin,blc):end,:),1);
            end
        end
    end

    mpow_alpha_all(u,:,:) = meanpow_alpha;
    mpow_beta_all(u,:,:) = meanpow_beta;
    clear power_spec;
    toc
end

%% to plot

malpha = squeeze(mean(mpow_alpha_all,1));
mbeta = squeeze(mean(mpow_beta_all,1));

A = reshape(malpha',[],1);
B = reshape(mbeta',[],1);

x = 1:75;
f = figure;
f.Position = [713,302,592,443];
scatter(x,A,30,'o','filled','MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[1 0 1]);hold on;
scatter(x,B,30,'bo','filled');
ylabel('z(Band power)','FontSize',12); xlabel('Trial bin','FontSize',12);
xlim([0 76]); ylim([-0.5 0.4]);
title('Post-stimulus','FontSize',13);
ax = gca;
ax.FontSize = 12; 
xline(15.5); xline(30.5); xline(45.5);xline(60.5);
legend('\alpha power','\beta power','','','','','Location','southeast'); legend('Box','off');


%%
clear all; clc; close all;
cd('C:\Users\MKopcanova\OneDrive - University of Dundee\MScDiss_data\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

 load '\experiment_data\Results\block_Ftest_n36_0.1Hz_prestim_ToT_fitlm.mat';
    maskhold1 = double(statsstruc.posclusterslabelmat == 1);
    maskhold2 = double(statsstruc.posclusterslabelmat == 2);
    freq1 = mean(maskhold1,1);
    freq2 = mean(maskhold2,1);
    freqsC1 = statsstruc.freq(1,find(freq1));
    freqsC2 = statsstruc.freq(1,find(freq2));
    els1 = find(mean(maskhold1,2));
    els2 = find(mean(maskhold2,2));

for u = 1:length(subnums1)

    tic
    subid1 = subnums1{u};
    load(['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_preStim_powerspec_0.1Hzres.mat']);
    datastruc = powerspec;
    clear powerspec;
    
    alpha_freq_start = find(datastruc.freq == round(freqsC1(1,1),1));
    alpha_freq_end = find(datastruc.freq == round(freqsC1(1,end),1));

    beta_freq_start = find(datastruc.freq == round(freqsC2(1,1),1));
    beta_freq_end = find(datastruc.freq == round(freqsC2(1,end),1));

    tt = length(datastruc.trialinfo);
    difficulty = zeros(1,tt);
    confidence = zeros(1,tt);
    accuracy = zeros(1,tt);
    trialorder = zeros(1,tt);
    rts = zeros(1,tt);
    block = zeros(1,tt);
  
    power_spec = datastruc.powspctrm;
    % Get difficulty, accuracy and confidence ratings for each trial

    for h = 1:tt
        confidence(1,h) = datastruc.trialinfo(h).confidence2;
        difficulty(1,h) = datastruc.trialinfo(h).SDTcond;
        accuracy(1,h) = datastruc.trialinfo(h).accuracy;
        trialorder(1,h) = datastruc.trialinfo(h).OriginalTrialNo;
        rts(1,h) = datastruc.trialinfo(h).RTresp;
        block(1,h) = datastruc.trialinfo(h).block;
    end

    blength = [length(find(block == 1)), length(find(block == 2)), length(find(block == 3)), length(find(block == 4)), length(find(block == 5))];
    bl_div = blength./15;
    bl_div_round = round(bl_div);


    alpha_pows = squeeze(mean(squeeze(mean(power_spec(:,els1,[alpha_freq_start : alpha_freq_end]),2)),2));
    beta_pows = squeeze(mean(squeeze(mean(power_spec(:,els2,[beta_freq_start : beta_freq_end]),2)),2));

    zalpha_pows = zscore(alpha_pows);
    zbeta_pows = zscore(beta_pows);
    
    for blc = 1:5
    bstart = 1; 
    bend = [];
        for bin = 1:14
            bstart = [bstart;bl_div_round(blc)*bin + 1];
        end
        binstart(:,blc) = bstart;
        for bin = 1:15
            bend = [bend;bl_div_round(blc)*bin];
        end
        binend(:,blc) = bend;
    end

    for blc = 1:5
        for bin = 1:15
            zalpha1 = zalpha_pows(find(block == blc),:);
            zbeta1 = zbeta_pows(find(block == blc),:);
            if bin == 15
               meanpow_alpha(blc,bin) = mean(zalpha1(binstart(bin,blc):end,:),1);
               meanpow_beta(blc,bin) = mean(zbeta1(binstart(bin,blc):end,:),1);
            else
               meanpow_alpha(blc,bin) = mean(zalpha1(binstart(bin,blc):binend(bin,blc),:),1);
               meanpow_beta(blc,bin) = mean(zbeta1(binstart(bin,blc):end,:),1);
            end
        end
    end

    mpow_alpha_all(u,:,:) = meanpow_alpha;
    mpow_beta_all(u,:,:) = meanpow_beta;
    clear power_spec;
    toc
end

%% to plot

malpha = squeeze(mean(mpow_alpha_all,1));
mbeta = squeeze(mean(mpow_beta_all,1));

A = reshape(malpha',[],1);
B = reshape(mbeta',[],1);

x = 1:75;
f = figure;
f.Position = [713,302,592,443];
scatter(x,A,30,'o','filled','MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[1 0 1]);hold on;
scatter(x,B,30,'bo','filled');
ylabel('z(Band power)','FontSize',12); xlabel('Trial bin','FontSize',12);
title('Pre-stimulus','FontSize',13);
xlim([0 76]); ylim([-0.5 0.4]);
ax = gca;
ax.FontSize = 12; 
xline(15.5); xline(30.5); xline(45.5);xline(60.5);
legend('\alpha power','\beta power','','','','','Location','southeast'); legend('Box','off');