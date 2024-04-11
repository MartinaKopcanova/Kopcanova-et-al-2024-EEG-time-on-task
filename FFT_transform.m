
%% ========================================================================================================
%  ========================================================================================================
%  ============         higher frequency resolution for peak freq analysis
%  ========================================================================================================

%% FFT on pre-stimulus data -1-0sec, 1-40Hz, 0.1Hz freq res
clear all; close all; clc;
addpath('C:\Program Files\fieldtrip-20220113');
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','762201','751325','792475','804448'};


for u = 1:length(subnums1)
    
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_S12_fieldNBR.mat']);
    
   for x = 1:length(data.trialinfo)
       data.time{x} = data.time{x}(1,769:1280);%pre-stim = -1-0s [769 1280]
       data.trial{x} = data.trial{x}(:,769:1280);
   end
   
    cfg = [];
    cfg.method= 'mtmfft';
    cfg.channel  = 'all';
    cfg.output = 'pow';
    cfg.keeptrials='yes';
    cfg.pad = 10;
    cfg.padtype = 'zero';
    cfg.foi = [1:0.1:40];
    cfg.taper = 'hanning';
    powerspec = ft_freqanalysis(cfg,data);
    % Attach behavioural data to data structure
    powerspec.trialinfo = data.trialinfo;
    savename = ['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_preStim_powerspec_0.1Hzres.mat'];
    save(savename, 'powerspec');
    
end

%% FFT on post-stimulus data 0-1s, 1-40Hz, 0.1Hz freq res
clear all; close all; clc;
addpath('C:\Program Files\fieldtrip-20220113');
cd('\experiment_data\EEG');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','762201','751325','792475','804448'};

for u = 1:length(subnums1)
    
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_S12_fieldNBR.mat']);
  
   for x = 1:length(data.trialinfo)
       data.time{x} = data.time{x}(1,1281:1793);%post-stim = 0-1s [1281 1793]
       data.trial{x} = data.trial{x}(:,1281:1793);
   end
   
    cfg = [];
    cfg.method= 'mtmfft';
    cfg.channel  = 'all';
    cfg.output = 'pow';
    cfg.keeptrials='yes';
    cfg.pad = 10;
    cfg.padtype = 'zero';
    cfg.foi = [1:0.1:40];
    cfg.taper = 'hanning';
    powerspec = ft_freqanalysis(cfg,data);
    % Attach behavioural data to data structure
    powerspec.trialinfo = data.trialinfo;
    savename = ['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_postStim_powerspec_0.1Hzres.mat'];
    save(savename, 'powerspec');
    
end

%%
%% FOOOF transfroms - full scalp
clear all; close all; clc;
cd('\experiment_data\EEG');
addpath('\fooof_mat-main\fooof_mat\');

%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','762201','751325','792475','804448'};

for u = 1:length(subnums1)
    tic
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FFT_preStim_powerspec_0.1Hzres.mat']);  
    for tt = 1:size(powerspec.powspctrm,1)
   
         % have a look at instructions for 'fooof.m' function
         pow = squeeze(mean(powerspec.powspctrm,2))';
         power = squeeze(pow(:,tt));
         freq = powerspec.freq;
    
         %settings
         settings = struct(); %default
         aperiodic_mode = ['fixed']; settings.aperiodic_mode = aperiodic_mode;  
         peak_threshold = [2]; settings.peak_threshold = peak_threshold; %in SD def=2
         %%setting the relative threshold for stopping the search for peaks
         %%down to 1SD of the flattend spectrum
         peak_width_limits = [2,15]; settings.peak_width_limits = peak_width_limits; %to fix the lower bounds warning
         max_n_peaks = [4]; settings.max_n_peaks = max_n_peaks; %defines max
         %number of peaks to be fitted
         %min_peak_height = []; settings.min_peak_height = min_peak_height;
         %%absolute threshold for stopping peak search as power over and above
         %%aperiodic component
         %fitting range ? 
         f_range = [3,40];
         fooof_results(tt) = fooof(freq, power, f_range,settings,true);
         offset(tt) = fooof_results(tt).aperiodic_params(1,1);
         exponent(tt) = fooof_results(tt).aperiodic_params(1,2);
         aperiodic_spec(tt,:) = fooof_results(tt).ap_fit;
         org_spec(tt,:) = fooof_results(tt).power_spectrum;
         error(tt) = fooof_results(tt).error;
         r_sqr(tt) = fooof_results(tt).r_squared;
    
    end
clear powerspec;
saveas = ['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_preStim_0.1Hz_ALL.mat'];
save(saveas, 'fooof_results','offset','exponent','org_spec','aperiodic_spec','error','r_sqr');
clear fooof_results offset exponent aperiodic_spec org_spec error r_sqr;
toc
end

%% Calculate aperiodic removed spectra & goodness of fit descriptives
%PRE/post-STIMULUS
clear all; close all; clc;
cd('\experiment_data\EEG');
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};


for u = 1:length(subnums1)
    tic
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_postStim_0.1Hz_ALL.mat']);
    merror = mean(error,2); mrsqr = mean(r_sqr,2);
    mmerr = mean(merror); mmrsqr = mean(mrsqr);
    minRS(u,:) = min(r_sqr);
    maxRS(u,:) = max(r_sqr);

    %dimord: trial-elec-freq
    for  aa = 1:size(org_spec,1)
%         for bb = 1:size(org_spec,2)
            spect_flat_log(aa,:) = org_spec(aa,:) - aperiodic_spec(aa,:); 
%         end
    end
    spect_flat = 10.^spect_flat_log;
    
    MeanErr(u) = mmerr;
    MeanRsqr(u) = mmrsqr;
    saveas = ['\experiment_data\EEG\',subid1,'\',subid1,'_FS_aperCorr_0.1Hz_postStim.mat'];
    save(saveas,'spect_flat','spect_flat_log');
    clear spect_flat_log spect_flat
end
save('\experiment_data\Results\fooof_FS_postStim_0.1Hz_GoF.mat', 'MeanRsqr','MeanErr','minRS','maxRS');

%% ========================================================================================================
%  ===============      get alpha peak freq from fooof output =============
% =========================================================================================================

clear all; close all; clc;
cd('\experiment_data\EEG');
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

for u = 1:length(subnums1)
    tic
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_FOOOF_FS_preStim_0.1Hz_ALL.mat']);

    maxAlp = zeros(1,length(fooof_results));
    peakCF = zeros(1,length(fooof_results));
    indAlp = zeros(1,length(fooof_results));
    for tt = 1:length(fooof_results)
        %get peak CF (with highest PW) falling into extended alpha range:
        %5-15Hz?
        for p = 1:length(fooof_results(tt).peak_params(:,1))            
            if (6 < fooof_results(tt).peak_params(p,1))&&(fooof_results(tt).peak_params(p,1) < 15)
                Peaks(tt).alpha(p,:) = fooof_results(tt).peak_params(p,:);  
                [maxAlp(tt), indAlp(tt)] = max(Peaks(tt).alpha(:,2));
                Peaks(tt).maxcf = Peaks(tt).alpha(indAlp(1,tt),1);
                peakCF(tt) = Peaks(tt).maxcf;
            end
        end
    end
    %number of trials with alpha peak detected
    nopeak_trials = find(peakCF == 0); earlyexclusion = find(nopeak_trials <= 450); prop_earlyex(u,:) = length(earlyexclusion)/length(nopeak_trials);
    alphatrials(u,:) = nnz(peakCF);
    excludedtrials(u,:) = length(fooof_results) - alphatrials(u,:);
    
    saveas = ['\experiment_data\EEG\',subid1,'\',subid1,'_PeakAlphaFreq_fooof_0.1Hz_preStim.mat'];
    save(saveas,"maxAlp","peakCF");
    clear Peaks maxAlp indAlp peakCF nopeak_trials earlyexclusion;
end

