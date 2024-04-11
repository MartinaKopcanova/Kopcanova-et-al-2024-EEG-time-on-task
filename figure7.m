%% ========================================================================
%  ==================  FIGURE 7 - TIMEFREQUENCY REGRESSION ================
%  ========================================================================
%%
%% =========================================================================
%============ Perform time-frequency transformation ======================
%=========================================================================
clear all; close all; clc;
cd('\experiment_data\EEG');
addpath('C:\Program Files\fieldtrip-20220113'); 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

for u = 1:length(subnums1)
    
    subid1 = subnums1{u};
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_S12_fieldNBR.mat']);
    cfg = [];
    cfg.method= 'mtmconvol';
    cfg.channel  = 'all';
    cfg.output = 'fourier';
    cfg.keeptrials='yes';
    cfg.foi = [1:40];
    cfg.t_ftimwin = ones(length(cfg.foi),1)*0.5;
    cfg.toi = -2.5:0.02:1.499;
    cfg.taper = 'hanning';
    timefreq = ft_freqanalysis(cfg,data);
    % Attach behavioural data to data structure
    timefreq.trialinfo = data.trialinfo;
    savename = ['\experiment_data\EEG\',subid1,'\',subid1,'_timefreq.mat'];
    save(savename, 'timefreq','-v7.3' );
    
end

%% regressions (with and without time-on-task included)
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
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_timefreq.mat']);
    datastruc = timefreq;
    clear timefreq;
    tim = datastruc.time;

    PreCut = datastruc.time(76); %-1sec
    PostCut = datastruc.time(176); %1sec
    PreCInd = 76;
    PostCInd = 176;
    %
    time = length(datastruc.time(1,PreCInd:PostCInd));

    datastruc.fourierspctrm = datastruc.fourierspctrm(:,:,:,PreCInd:PostCInd);


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
   
    % Loop through electrodes and times performing regression
    % (Power versus accuracy/ratings)
    powertransform = abs(datastruc.fourierspctrm).*abs(datastruc.fourierspctrm);
    powertransformdb = pow2db(powertransform);

    for mm = 1:length(datastruc.label)
        for bb = 1:length(datastruc.freq)
            for uu = 1:time

                % everything rank ordered apart from accuracy and confidence!!!!
                jackit = squeeze(powertransformdb(:,mm,bb,uu))';
                jackit_rank = tiedrank(jackit);
%                 ratings_rank = tiedrank(confidence);
                diff_rank = tiedrank(difficulty);
                trialorder_rank = tiedrank(trialorder);
                rt_rank = tiedrank(rts);

                % create data table
                Regr_table = table(diff_rank',jackit_rank',confidence',accuracy',...
                    trialorder_rank',rt_rank','VariableNames',{'diff_rank','jackit_rank','confidence','accuracy','trialorder_rank','rt_rank'});

                %             CONFIDENCE MODEL CONTROLLING FOR DIFFICULTY AND ACCURACY
                modelspec = 'jackit_rank ~ confidence + diff_rank + accuracy + rt_rank';
                mdl = fitlm(Regr_table,modelspec);
      
                get_conf = mdl.Coefficients(3,1);
                conf_EEG(mm,bb,uu) = table2array(get_conf);
                
                get_diff = mdl.Coefficients(2,1);
                diff_EEG(mm,bb,uu) = table2array(get_diff);

                get_acc = mdl.Coefficients(4,1);
                acc_EEG(mm,bb,uu) = table2array(get_acc);

                get_rt = mdl.Coefficients(5,1);
                rt_EEG(mm,bb,uu) = table2array(get_rt);

%                 get_trial = mdl.Coefficients(5,1);
%                 trial_EEG(mm,bb,uu) = table2array(get_trial);

          
        
                clear jackit jackit_rank Regr_table diff_rank rt_rank trialorder_rank;
            end
        end
    end
    toc;
%     %Save beta matrices
    saveas = ['\experiment_data\Results\',subid1,'_RegressionTF_noToT_1s_LowFreqs.mat'];
    save(saveas, 'conf_EEG', 'rt_EEG','diff_EEG','acc_EEG');%
    clear datastruc conf_EEG rt_EEG diff_EEG acc_EEG ; %
end

%% 
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
    load (['\experiment_data\EEG\',subid1,'\',subid1,'_timefreq.mat']);
    datastruc = timefreq;
    clear timefreq;
    tim = datastruc.time;

    %edit this based on what I want the post-stim cut off time to be!! - or
    %could even add a pre-stim one if don't shorten the epochs to -1s
    %beforehand
    PreCut = datastruc.time(76); %-1sec
    PostCut = datastruc.time(176); %1sec
    PreCInd = 76;
    PostCInd = 176;
    %
    time = length(datastruc.time(1,PreCInd:PostCInd));

    datastruc.fourierspctrm = datastruc.fourierspctrm(:,:,:,PreCInd:PostCInd);


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
   
    % Loop through electrodes and times performing regression
    % (Power versus accuracy/ratings)
    powertransform = abs(datastruc.fourierspctrm).*abs(datastruc.fourierspctrm);
    powertransformdb = pow2db(powertransform);

    for mm = 1:length(datastruc.label)
        for bb = 1:length(datastruc.freq)
            for uu = 1:time

                % everything rank ordered apart from accuracy and confidence!!!!
                jackit = squeeze(powertransformdb(:,mm,bb,uu))';
                jackit_rank = tiedrank(jackit);
%                 ratings_rank = tiedrank(confidence);
                diff_rank = tiedrank(difficulty);
                trialorder_rank = tiedrank(trialorder);
                rt_rank = tiedrank(rts);

                % create data table
                Regr_table = table(diff_rank',jackit_rank',confidence',accuracy',...
                    trialorder_rank',rt_rank','VariableNames',{'diff_rank','jackit_rank','confidence','accuracy','trialorder_rank','rt_rank'});

                %             CONFIDENCE MODEL CONTROLLING FOR DIFFICULTY AND ACCURACY
                modelspec = 'jackit_rank ~ confidence + diff_rank + accuracy + rt_rank + trialorder_rank';
                mdl = fitlm(Regr_table,modelspec);
      
                get_conf = mdl.Coefficients(3,1);
                conf_EEG(mm,bb,uu) = table2array(get_conf);
                
                get_diff = mdl.Coefficients(2,1);
                diff_EEG(mm,bb,uu) = table2array(get_diff);

                get_acc = mdl.Coefficients(4,1);
                acc_EEG(mm,bb,uu) = table2array(get_acc);

                get_rt = mdl.Coefficients(6,1);
                rt_EEG(mm,bb,uu) = table2array(get_rt);

                get_trial = mdl.Coefficients(5,1);
                trial_EEG(mm,bb,uu) = table2array(get_trial);

          
        
                clear jackit jackit_rank Regr_table diff_rank rt_rank trialorder_rank;
            end
        end
    end
    toc;
%     %Save beta matrices
    saveas = ['\experiment_data\Results\',subid1,'_RegressionTF_1s_LowFreqs.mat'];
    save(saveas, 'conf_EEG', 'rt_EEG','diff_EEG','acc_EEG','trial_EEG');%
    clear datastruc conf_EEG rt_EEG diff_EEG acc_EEG trial_EEG; %
end

%%
clear all; clc; close all;
cd('\experiment_data\Results');
 
%define subject;
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};


for I = 1:length(subnums1)
    
    subid1 = subnums1{I};
    load (['\experiment_data\Results\',subid1,'_RegressionTF_noToT_1s_LowFreqs.mat']);
    
    for mm = 1:32
        for bb = 1:size(conf_EEG,2)
            for uu = 1:size(conf_EEG,3)

                beta_conf(I,mm,bb,uu) = conf_EEG(mm,bb,uu);
                beta_acc(I,mm,bb,uu) = acc_EEG(mm,bb,uu);
                beta_diff(I,mm,bb,uu) = diff_EEG(mm,bb,uu);
                beta_rt(I,mm,bb,uu) = rt_EEG(mm,bb,uu);
%                 beta_trial(I,mm,bb,uu) = trial_EEG(mm,bb,uu);
            end
        end
    end
    clear  conf_EEG acc_EEG diff_EEG rt_EEG trial_EEG;
end
save('\experiment_data\Results\RegressionTF_Power_LowFreqs_allbeh_noToT_n36.mat',  "beta_rt", "beta_diff","beta_acc","beta_conf");%,"beta_trial");

%%
clear all; close all; clc;
addpath('C:\Program Files\fieldtrip-20220113'); 
load('\experiment_data\EEG\157061\157061_timefreq.mat');
cfgN = [];
cfgN.method = 'template';
cfgN.template = 'biosemi32_neighb.mat';
neighbours = ft_prepare_neighbours(cfgN, timefreq);
%save('D:\Dissertation_Data\EEG\neighbours', 'neighbours');
time = timefreq.time;
labels = timefreq.elec.label;
clear timefreq;
load('\experiment_data\Results\RegressionTF_Power_LowFreqs_allbeh_ToT_n36.mat');
jackzeros = zeros(36,32,40,size(beta_conf,4));
entrydata = { 'beta_conf','beta_acc','beta_diff','beta_rt','beta_trial'};
savdat = {'confidence_stats','accuracy_stats','difficulty_stats','rt_stats','trial_stats'};

for entd = 1:length(entrydata)
    powerent.label = labels;
    powerent.dimord = 'subj_chan_freq_time';
    powerent.freq = 1:40;
    powerent.time = time(76:176);
    % Perform t-test
    cfg = [];
    cfg.statistic = 'ft_statfun_depsamplesT';
    desig1 = [ones(36,1)', 2*ones(36,1)'];
    desig2 = [1:36,1:36];
    desig = [desig1;desig2];
    cfg.latency = 'all';%
    cfg.frequency = [1 40];
    cfg.ivar = 1;
    cfg.uvar = 2;
    cfg.design = desig;
    cfg.keeptrials = 'no';
    cfg.avgovertime = 'no';
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
    
    savename = ['\experiment_data\Results\',savdat{entd},'_n36_1Hz_1000ms_LowFreqs_ToT_fitlm.mat'];
    save(savename, 'statsstruc');
    clear statsstruc;
end

%% plots
%% with time-on-task results
clear all; close all; clc;
cd('\experiment_data\Results\');

allconds = {'confidence_stats_n36_1Hz_1000ms_LowFreqs_ToT_fitlm','accuracy_stats_n36_1Hz_1000ms_LowFreqs_ToT_fitlm',...
    'difficulty_stats_n36_1Hz_1000ms_LowFreqs_ToT_fitlm','rt_stats_n36_1Hz_1000ms_LowFreqs_ToT_fitlm',...
    'trial_stats_n36_1Hz_1000ms_LowFreqs_ToT_fitlm'};
titles = {'Confidence','Accuracy','Difficulty','Response RT','Time-on-Task'};
for ccc = 1:length(allconds)

        load(['\experiment_data\Results\',allconds{ccc},'.mat']);
        comtype = allconds{ccc};
        
        sta = statsstruc.stat(:,:,:);%
        stap = statsstruc.mask(:,:,:);%
        
        %for significant stats
        if ccc == 2
            stapaware = statsstruc.mask(:,:,1:50);
        else
        end
        
        subplot(2,3,ccc)
%         fig_h=figure;
%         pos = get(gcf, 'Position');
%         width = 3.15;     % Width in inches (try to match journal requirements?)
%         height = 3.15;    % Height in inches
%         set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
%         subplot('Position',[0.29 0.2500 0.55 0.55])
        rast = sta;
        colors = cbrewer( 'RdBu',64); %colors = cbrewer('div', 'RdBu',64);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        negpos = {'Positive','Negative'};
        allrast = squeeze(mean(rast,1));
        multin = squeeze(mean(stap,1));
        imagesc(allrast); hold on
        caxis([-3.0 3]);
        % Highlight significant clusters on t-map
        [B,L,N,A] = bwboundaries(multin);
        %                 colors=[0.5 0.5 0.5];
        for k=1:length(B)
            boundary = B{k};
           % cidx = mod(k,length(colors))+1;
            plot(boundary(:,2), boundary(:,1),...
                'Color',[0 0 0],'LineWidth',0.75);
        end
%         test = 0;
        %             end
        hold on;
       set(gca,'Ytick',[1 5 10 15 20 25 30 35 40],'YTickLabel',[1 5 10 15 20 25 30 35 40]);
        set(gca, 'XTick',[1 26 51 76 101],'XTickLabel',[-1,-0.5,0,0.5,1]);%,'Position',coor);
       set(gca, 'fontsize',10);
       set(gca, 'ydir', 'normal');
       if ccc == 2
       ylabel('Frequency (Hz)');
       xlabel('Time (sec)');
       end
       xline(51,'k-');
       %xtickangle(70);
        hold on;
        title(titles{ccc});
        hold on;
%         if ccc == 4
%         handles = colorbar;
%         handles.Limits = [-3.0 3];
%         handles.TickDirection = 'out';
%         handles.Box = 'on';
%         handles.Label.String = 'T-Value';
%         handles.FontSize = 8;
%         handles.Label.FontSize = 9;
%         handles.Position = [0.87 0.25 0.05 0.55];
%         handles.Label.Position = [1.7 0 0];
%         end
        drawnow;
end
%% without time-on-task results
clear all; close all; clc;
cd('C:\Users\7mart\OneDrive - University of Dundee\MScDiss_data\experiment_data\Results\');

allconds = {'confidence_stats_noToT_n36_1Hz_1000ms_LowFreqs_ToT_fitlm','accuracy_stats_noToT_n36_1Hz_1000ms_LowFreqs_ToT_fitlm',...
    'difficulty_stats_noToT_n36_1Hz_1000ms_LowFreqs_ToT_fitlm','rt_stats_noToT_n36_1Hz_1000ms_LowFreqs_ToT_fitlm'};
titles = {'Confidence','Accuracy','Difficulty','Response RT'};
for ccc = 1:length(allconds)

        load(['\experiment_data\Results\',allconds{ccc},'.mat']);
        comtype = allconds{ccc};
        
        sta = statsstruc.stat(:,:,:);%
        stap = statsstruc.mask(:,:,:);%
        
        %for significant stats
        if ccc == 2
            stapaware = statsstruc.mask(:,:,1:50);
        else
        end
        
        subplot(2,3,ccc)
%         fig_h=figure;
%         pos = get(gcf, 'Position');
%         width = 3.15;     % Width in inches (try to match journal requirements?)
%         height = 3.15;    % Height in inches
%         set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
%         subplot('Position',[0.29 0.2500 0.55 0.55])
        rast = sta;
        colors = cbrewer( 'RdBu',64); %colors = cbrewer('div', 'RdBu',64);
        colors = flipud(colors); % puts red on top, blue at the bottom
        colormap(colors);
        negpos = {'Positive','Negative'};
        allrast = squeeze(mean(rast,1));
        multin = squeeze(mean(stap,1));
        imagesc(allrast); hold on
        caxis([-3.0 3]);
        % Highlight significant clusters on t-map
        [B,L,N,A] = bwboundaries(multin);
        %                 colors=[0.5 0.5 0.5];
        for k=1:length(B)
            boundary = B{k};
           % cidx = mod(k,length(colors))+1;
            plot(boundary(:,2), boundary(:,1),...
                'Color',[0 0 0],'LineWidth',0.75);
        end
%         test = 0;
        %             end
        hold on;
       set(gca,'Ytick',[1 5 10 15 20 25 30 35 40],'YTickLabel',[1 5 10 15 20 25 30 35 40]);
        set(gca, 'XTick',[1 26 51 76 101],'XTickLabel',[-1,-0.5,0,0.5,1]);%,'Position',coor);
       set(gca, 'fontsize',10);
       set(gca, 'ydir', 'normal');
       if ccc == 2
       ylabel('Frequency (Hz)');
       xlabel('Time (sec)');
       end
       xline(51,'k-');
       %xtickangle(70);
        hold on;
        title(titles{ccc});
        hold on;
%         if ccc == 4
%         handles = colorbar;
%         handles.Limits = [-3.0 3];
%         handles.TickDirection = 'out';
%         handles.Box = 'on';
%         handles.Label.String = 'T-Value';
%         handles.FontSize = 8;
%         handles.Label.FontSize = 9;
%         handles.Position = [0.87 0.25 0.05 0.55];
%         handles.Label.Position = [1.7 0 0];
%         end
        drawnow;
end

