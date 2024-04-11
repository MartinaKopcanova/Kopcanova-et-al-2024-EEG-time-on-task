%% ===============================================================================
%  ===============================================================================
%                           Behaviour plotted by the 5 difficulty levels
%                           and split by blocks
% =================================================================================
clear all; close all; clc;
lazypath = '\experiment_data\Behaviour\';
cd(lazypath);

Folders = dir(lazypath); %% NAME OF FOLDERS
Folders(1:2) = [];
for i = 1:length(Folders)
    FolderName = Folders(i).name;
    datafolder2 = [lazypath ,'\', FolderName]; %
    cd(datafolder2)% go inside the folder 
    subject = FolderName;
    Ind = FolderName; % here is just to find the name of the subject in case your folders are super long name e.i. find(FolderName == 'S'); %name of your folder 
    subjectname = subject; % in case you need to take some stuff out (1:Ind-4);
    
    cd([lazypath,'\',FolderName]);
    behData = dir('*_1.csv');
    beh = readtable(behData.name);


    block = [ones(180,1); ones(180,1)*2;ones(180,1)*3;ones(180,1)*4;ones(180,1)*5];
    accuracy = beh.key_resp_corr;
    contrast = beh.contrast;
    accuracy(isnan(accuracy))=[];
    contrast(isnan(contrast))=[];
    confidence = beh.key_resp_2_keys; confidence = confidence(~cellfun('isempty',confidence));
    respRT = beh.key_resp_rt; respRT(isnan(respRT)) = [];
    confrespRT = beh.key_resp_2_rt; confrespRT(isnan(confrespRT)) = [];
    for ii = 1:length(confidence)
        if confidence(ii,:) == "k"
            confidence2(ii,:) = 1;
        elseif confidence(ii,:) == "l"
            confidence2(ii,:) = 0;
        end
    end


    %behaviour by difficulty
    oneb1(i) = mean(accuracy(find(contrast == 0.075 & block == 1),:));
    oneb2(i) = mean(accuracy(find(contrast == 0.075 & block == 2),:));
    oneb3(i) = mean(accuracy(find(contrast == 0.075 & block == 3),:));
    oneb4(i) = mean(accuracy(find(contrast == 0.075 & block == 4),:));
    oneb5(i) = mean(accuracy(find(contrast == 0.075 & block == 5),:));

    twob1(i) = mean(accuracy(find(contrast == 0.15 & block == 1),:));
    twob2(i) = mean(accuracy(find(contrast == 0.15 & block == 2),:));
    twob3(i) = mean(accuracy(find(contrast == 0.15 & block == 3),:));
    twob4(i) = mean(accuracy(find(contrast == 0.15 & block == 4),:));
    twob5(i) = mean(accuracy(find(contrast == 0.15 & block == 5),:));

    threeb1(i) = mean(accuracy(find(contrast == 0.225 & block == 1),:));
    threeb2(i) = mean(accuracy(find(contrast == 0.225 & block == 2),:));
    threeb3(i) = mean(accuracy(find(contrast == 0.225 & block == 3),:));
    threeb4(i) = mean(accuracy(find(contrast == 0.225 & block == 4),:));
    threeb5(i) = mean(accuracy(find(contrast == 0.225 & block == 5),:));

    fourb1(i) = mean(accuracy(find(contrast == 0.3 & block == 1),:));
    fourb2(i) = mean(accuracy(find(contrast == 0.3 & block == 2),:));
    fourb3(i) = mean(accuracy(find(contrast == 0.3 & block == 3),:));
    fourb4(i) = mean(accuracy(find(contrast == 0.3 & block == 4),:));
    fourb5(i) = mean(accuracy(find(contrast == 0.3 & block == 5),:));

    fiveb1(i) = mean(accuracy(find(contrast == 0.375 & block == 1),:));
    fiveb2(i) = mean(accuracy(find(contrast == 0.375 & block == 2),:));
    fiveb3(i) = mean(accuracy(find(contrast == 0.375 & block == 3),:));
    fiveb4(i) = mean(accuracy(find(contrast == 0.375 & block == 4),:));
    fiveb5(i) = mean(accuracy(find(contrast == 0.375 & block == 5),:));
    

    %for the inset with collapsed difficulty across blocks
    b1(i) = mean(accuracy(find( block == 1),:));
    b2(i) = mean(accuracy(find( block == 2),:));
    b3(i) = mean(accuracy(find( block == 3),:));
    b4(i) = mean(accuracy(find( block == 4),:));
    b5(i) = mean(accuracy(find( block == 5),:));
    

end
    all_type1_b1 = [oneb1;twob1;threeb1;fourb1;fiveb1];
    all_type1_b2 = [oneb2;twob2;threeb2;fourb2;fiveb2];
    all_type1_b3 = [oneb3;twob3;threeb3;fourb3;fiveb3];
    all_type1_b4 = [oneb4;twob4;threeb4;fourb4;fiveb4];
    all_type1_b5 = [oneb5;twob5;threeb5;fourb5;fiveb5];
    
    mean_t1_b1 = mean(all_type1_b1,2);
    mean_t1_b2 = mean(all_type1_b2,2);
    mean_t1_b3 = mean(all_type1_b3,2);
    mean_t1_b4 = mean(all_type1_b4,2);
    mean_t1_b5 = mean(all_type1_b5,2);

    sem1_b1 = std(all_type1_b1(1,:))./sqrt(36);
    sem2_b1 = std(all_type1_b1(2,:))./sqrt(36);
    sem3_b1 = std(all_type1_b1(3,:))./sqrt(36);
    sem4_b1 = std(all_type1_b1(4,:))./sqrt(36);
    sem5_b1 = std(all_type1_b1(5,:))./sqrt(36);
    
    sem1_b2 = std(all_type1_b2(1,:))./sqrt(36);
    sem2_b2 = std(all_type1_b2(2,:))./sqrt(36);
    sem3_b2 = std(all_type1_b2(3,:))./sqrt(36);
    sem4_b2 = std(all_type1_b2(4,:))./sqrt(36);
    sem5_b2 = std(all_type1_b2(5,:))./sqrt(36);
    
    sem1_b3 = std(all_type1_b3(1,:))./sqrt(36);
    sem2_b3 = std(all_type1_b3(2,:))./sqrt(36);
    sem3_b3 = std(all_type1_b3(3,:))./sqrt(36);
    sem4_b3 = std(all_type1_b3(4,:))./sqrt(36);
    sem5_b3 = std(all_type1_b3(5,:))./sqrt(36);

    sem1_b4 = std(all_type1_b4(1,:))./sqrt(36);
    sem2_b4 = std(all_type1_b4(2,:))./sqrt(36);
    sem3_b4 = std(all_type1_b4(3,:))./sqrt(36);
    sem4_b4 = std(all_type1_b4(4,:))./sqrt(36);
    sem5_b4 = std(all_type1_b4(5,:))./sqrt(36);

    sem1_b5 = std(all_type1_b5(1,:))./sqrt(36);
    sem2_b5 = std(all_type1_b5(2,:))./sqrt(36);
    sem3_b5 = std(all_type1_b5(3,:))./sqrt(36);
    sem4_b5 = std(all_type1_b5(4,:))./sqrt(36);
    sem5_b5 = std(all_type1_b5(5,:))./sqrt(36);
  
% STATS - 5X5 WITHINPARTICIPANT ANOVA with block and difficulty as IV and
% performance as DV
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

t = table(subnums1',all_type1_b1(1,:)',all_type1_b1(2,:)',all_type1_b1(3,:)',all_type1_b1(4,:)',all_type1_b1(5,:)',...
    all_type1_b2(1,:)',all_type1_b2(2,:)',all_type1_b2(3,:)',all_type1_b2(4,:)',all_type1_b2(5,:)',...
    all_type1_b3(1,:)',all_type1_b3(2,:)',all_type1_b3(3,:)',all_type1_b3(4,:)',all_type1_b3(5,:)',...
    all_type1_b4(1,:)',all_type1_b4(2,:)',all_type1_b4(3,:)',all_type1_b4(4,:)',all_type1_b4(5,:)',....
    all_type1_b5(1,:)',all_type1_b5(2,:)',all_type1_b5(3,:)',all_type1_b5(4,:)',all_type1_b5(5,:)',...
    'VariableNames',...
    {'IDs','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12','v13','v14','v15',...
    'v16','v17','v18','v19','v20','v21','v22','v23','v24','v25'});

withinDesign = table([ones(1,5), ones(1,5)*2, ones(1,5)*3,ones(1,5)*4, ones(1,5)*5]',[1:5,1:5,1:5,1:5,1:5]','VariableNames',{'Block','Difficulty'});
withinDesign.Block = categorical(withinDesign.Block);
withinDesign.Difficulty = categorical(withinDesign.Difficulty);

rm = fitrm(t, 'v1-v25 ~ 1', 'WithinDesign', withinDesign);
%mauchly(rm)
AT = ranova(rm, 'WithinModel', 'Block*Difficulty');
%for GG corrected df:
tbl = epsilon(rm); eps = tbl.GreenhouseGeisser;
AT.df_gg = AT.DF .* tbl.GreenhouseGeisser;
%%   
ClrMap = colormap('cool');
Clr = [42, 84, 126, 168, 210];
x = 1:5;
figure;
plot(x,mean_t1_b1,'.-','Color',ClrMap(Clr(1),:),'LineWidth',1.5);hold on;
er = errorbar([1 2 3 4 5],mean_t1_b1,[sem1_b1 sem2_b1 sem3_b1 sem4_b1 sem5_b1],[sem1_b1 sem2_b1 sem3_b1 sem4_b1 sem5_b1]);
er.Color = ClrMap(Clr(1),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b2,'.-','Color',ClrMap(Clr(2),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b2,[sem1_b2 sem2_b2 sem3_b2 sem4_b2 sem5_b2],[sem1_b2 sem2_b2 sem3_b2 sem4_b2 sem5_b2]);
er.Color = ClrMap(Clr(2),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b3,'.-','Color',ClrMap(Clr(3),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b3,[sem1_b3 sem2_b3 sem3_b3 sem4_b3 sem5_b3],[sem1_b3 sem2_b3 sem3_b3 sem4_b3 sem5_b3]);
er.Color = ClrMap(Clr(3),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b4,'.-','Color',ClrMap(Clr(4),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b4,[sem1_b4 sem2_b4 sem3_b4 sem4_b4 sem5_b4],[sem1_b4 sem2_b4 sem3_b4 sem4_b4 sem5_b4]);
er.Color = ClrMap(Clr(4),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b5,'.-','Color',ClrMap(Clr(5),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b5,[sem1_b5 sem2_b5 sem3_b5 sem4_b5 sem5_b5],[sem1_b5 sem2_b5 sem3_b5 sem4_b5 sem5_b5]);
er.Color = ClrMap(Clr(5),:);                            
er.LineStyle = 'none';
title('Mean type-1 PF'); ylim([0.45 1]);
xlim([0.8 5.2]);
% legend('Block 1','','Block 2','','Block 3','','Block 4','','Block 5','','location','southeast');
% legend('box','off');
ylabel('Proportion correct'); xlabel('Stimulus contrast (low to high)');hold off;
xticks([1 2 3 4 5]); xticklabels([0.075 0.15 0.225 0.3 0.375]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12)
%inset
% axes('Position',[.55 .17 .3 .4]); box off;
% set(gca,'FontSize',14);
figure;
params = [b1';b2';b3';b4';b5'];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean proportion correct'); 
hold on;
xlabel('Block');
ylim([0.45 1]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',11)
%% =================================================================================
%                       accuracy by block only to follow up the sig ANOVA
%                       block effect
% ===================================================================================
clear all; close all; clc;
lazypath = '\experiment_data\Behaviour\';
cd(lazypath);

Folders = dir(lazypath); %% NAME OF FOLDERS
Folders(1:2) = [];
for i = 1:length(Folders)
    FolderName = Folders(i).name;
    datafolder2 = [lazypath ,'\', FolderName]; %
    cd(datafolder2)% go inside the folder 
    subject = FolderName;
    Ind = FolderName; % here is just to find the name of the subject in case your folders are super long name e.i. find(FolderName == 'S'); %name of your folder 
    subjectname = subject; % in case you need to take some stuff out (1:Ind-4);
    
    cd([lazypath,'\',FolderName]);
    behData = dir('*_1.csv');
    beh = readtable(behData.name);


    block = [ones(180,1); ones(180,1)*2;ones(180,1)*3;ones(180,1)*4;ones(180,1)*5];
    accuracy = beh.key_resp_corr;
    contrast = beh.contrast;
    accuracy(isnan(accuracy))=[];
    contrast(isnan(contrast))=[];
    confidence = beh.key_resp_2_keys; confidence = confidence(~cellfun('isempty',confidence));
    respRT = beh.key_resp_rt; respRT(isnan(respRT)) = [];
    confrespRT = beh.key_resp_2_rt; confrespRT(isnan(confrespRT)) = [];
    for ii = 1:length(confidence)
        if confidence(ii,:) == "k"
            confidence2(ii,:) = 1;
        elseif confidence(ii,:) == "l"
            confidence2(ii,:) = 0;
        end
    end


    %behaviour by difficulty
    b1(i) = mean(accuracy(find( block == 1),:));
    b2(i) = mean(accuracy(find( block == 2),:));
    b3(i) = mean(accuracy(find( block == 3),:));
    b4(i) = mean(accuracy(find( block == 4),:));
    b5(i) = mean(accuracy(find( block == 5),:));
end

ClrMap = colormap('cool');
Clr = [42, 84, 126, 168, 210];
figure;
params = [b1';b2';b3';b4';b5'];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean proportion correct'); hold on;
xlabel('Block');
% y = [1,2];
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; 

subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};
%stats (ANOVA?)
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',b1',b2',b3',b4',b5','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(2))
c = multcompare(rm,'Blocks')
tbl = epsilon(rm); eps = tbl.GreenhouseGeisser;
ranovatbl.df_gg = ranovatbl.DF .* tbl.GreenhouseGeisser
%% ===============================================================================================================
%                       CONFIDENCE
% ================================================================================================================
clear all; close all; clc;
lazypath = '\experiment_data\Behaviour\';
cd(lazypath);

Folders = dir(lazypath); %% NAME OF FOLDERS
Folders(1:2) = [];
for i = 1:length(Folders)
    FolderName = Folders(i).name;
    datafolder2 = [lazypath ,'\', FolderName]; %
    cd(datafolder2)% go inside the folder 
    subject = FolderName;
    Ind = FolderName; % here is just to find the name of the subject in case your folders are super long name e.i. find(FolderName == 'S'); %name of your folder 
    subjectname = subject; % in case you need to take some stuff out (1:Ind-4);
    
    cd([lazypath,'\',FolderName]);
    behData = dir('*_1.csv');
    beh = readtable(behData.name);


    block = [ones(180,1); ones(180,1)*2;ones(180,1)*3;ones(180,1)*4;ones(180,1)*5];
    accuracy = beh.key_resp_corr;
    contrast = beh.contrast;
    accuracy(isnan(accuracy))=[];
    contrast(isnan(contrast))=[];
    confidence = beh.key_resp_2_keys; confidence = confidence(~cellfun('isempty',confidence));
    respRT = beh.key_resp_rt; respRT(isnan(respRT)) = [];
    confrespRT = beh.key_resp_2_rt; confrespRT(isnan(confrespRT)) = [];
    for ii = 1:length(confidence)
        if confidence(ii,:) == "k"
            confidence2(ii,:) = 1;
        elseif confidence(ii,:) == "l"
            confidence2(ii,:) = 0;
        end
    end


    %behaviour by difficulty
    oneb1(i) = mean(confidence2(find(contrast == 0.075 & block == 1),:));
    oneb2(i) = mean(confidence2(find(contrast == 0.075 & block == 2),:));
    oneb3(i) = mean(confidence2(find(contrast == 0.075 & block == 3),:));
    oneb4(i) = mean(confidence2(find(contrast == 0.075 & block == 4),:));
    oneb5(i) = mean(confidence2(find(contrast == 0.075 & block == 5),:));

    twob1(i) = mean(confidence2(find(contrast == 0.15 & block == 1),:));
    twob2(i) = mean(confidence2(find(contrast == 0.15 & block == 2),:));
    twob3(i) = mean(confidence2(find(contrast == 0.15 & block == 3),:));
    twob4(i) = mean(confidence2(find(contrast == 0.15 & block == 4),:));
    twob5(i) = mean(confidence2(find(contrast == 0.15 & block == 5),:));

    threeb1(i) = mean(confidence2(find(contrast == 0.225 & block == 1),:));
    threeb2(i) = mean(confidence2(find(contrast == 0.225 & block == 2),:));
    threeb3(i) = mean(confidence2(find(contrast == 0.225 & block == 3),:));
    threeb4(i) = mean(confidence2(find(contrast == 0.225 & block == 4),:));
    threeb5(i) = mean(confidence2(find(contrast == 0.225 & block == 5),:));

    fourb1(i) = mean(confidence2(find(contrast == 0.3 & block == 1),:));
    fourb2(i) = mean(confidence2(find(contrast == 0.3 & block == 2),:));
    fourb3(i) = mean(confidence2(find(contrast == 0.3 & block == 3),:));
    fourb4(i) = mean(confidence2(find(contrast == 0.3 & block == 4),:));
    fourb5(i) = mean(confidence2(find(contrast == 0.3 & block == 5),:));

    fiveb1(i) = mean(confidence2(find(contrast == 0.375 & block == 1),:));
    fiveb2(i) = mean(confidence2(find(contrast == 0.375 & block == 2),:));
    fiveb3(i) = mean(confidence2(find(contrast == 0.375 & block == 3),:));
    fiveb4(i) = mean(confidence2(find(contrast == 0.375 & block == 4),:));
    fiveb5(i) = mean(confidence2(find(contrast == 0.375 & block == 5),:));

%
    b1(i) = mean(confidence2(find( block == 1),:));
    b2(i) = mean(confidence2(find( block == 2),:));
    b3(i) = mean(confidence2(find( block == 3),:));
    b4(i) = mean(confidence2(find( block == 4),:));
    b5(i) = mean(confidence2(find( block == 5),:));

    

end
    all_type1_b1 = [oneb1;twob1;threeb1;fourb1;fiveb1];
    all_type1_b2 = [oneb2;twob2;threeb2;fourb2;fiveb2];
    all_type1_b3 = [oneb3;twob3;threeb3;fourb3;fiveb3];
    all_type1_b4 = [oneb4;twob4;threeb4;fourb4;fiveb4];
    all_type1_b5 = [oneb5;twob5;threeb5;fourb5;fiveb5];
    
    mean_t1_b1 = mean(all_type1_b1,2);
    mean_t1_b2 = mean(all_type1_b2,2);
    mean_t1_b3 = mean(all_type1_b3,2);
    mean_t1_b4 = mean(all_type1_b4,2);
    mean_t1_b5 = mean(all_type1_b5,2);

    sem1_b1 = std(all_type1_b1(1,:))./sqrt(36);
    sem2_b1 = std(all_type1_b1(2,:))./sqrt(36);
    sem3_b1 = std(all_type1_b1(3,:))./sqrt(36);
    sem4_b1 = std(all_type1_b1(4,:))./sqrt(36);
    sem5_b1 = std(all_type1_b1(5,:))./sqrt(36);
    
    sem1_b2 = std(all_type1_b2(1,:))./sqrt(36);
    sem2_b2 = std(all_type1_b2(2,:))./sqrt(36);
    sem3_b2 = std(all_type1_b2(3,:))./sqrt(36);
    sem4_b2 = std(all_type1_b2(4,:))./sqrt(36);
    sem5_b2 = std(all_type1_b2(5,:))./sqrt(36);
    
    sem1_b3 = std(all_type1_b3(1,:))./sqrt(36);
    sem2_b3 = std(all_type1_b3(2,:))./sqrt(36);
    sem3_b3 = std(all_type1_b3(3,:))./sqrt(36);
    sem4_b3 = std(all_type1_b3(4,:))./sqrt(36);
    sem5_b3 = std(all_type1_b3(5,:))./sqrt(36);

    sem1_b4 = std(all_type1_b4(1,:))./sqrt(36);
    sem2_b4 = std(all_type1_b4(2,:))./sqrt(36);
    sem3_b4 = std(all_type1_b4(3,:))./sqrt(36);
    sem4_b4 = std(all_type1_b4(4,:))./sqrt(36);
    sem5_b4 = std(all_type1_b4(5,:))./sqrt(36);

    sem1_b5 = std(all_type1_b5(1,:))./sqrt(36);
    sem2_b5 = std(all_type1_b5(2,:))./sqrt(36);
    sem3_b5 = std(all_type1_b5(3,:))./sqrt(36);
    sem4_b5 = std(all_type1_b5(4,:))./sqrt(36);
    sem5_b5 = std(all_type1_b5(5,:))./sqrt(36);
  
% STATS - 5X5 WITHINPARTICIPANT ANOVA with block and difficulty as IV and
% performance as DV
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

t = table(subnums1',all_type1_b1(1,:)',all_type1_b1(2,:)',all_type1_b1(3,:)',all_type1_b1(4,:)',all_type1_b1(5,:)',...
    all_type1_b2(1,:)',all_type1_b2(2,:)',all_type1_b2(3,:)',all_type1_b2(4,:)',all_type1_b2(5,:)',...
    all_type1_b3(1,:)',all_type1_b3(2,:)',all_type1_b3(3,:)',all_type1_b3(4,:)',all_type1_b3(5,:)',...
    all_type1_b4(1,:)',all_type1_b4(2,:)',all_type1_b4(3,:)',all_type1_b4(4,:)',all_type1_b4(5,:)',....
    all_type1_b5(1,:)',all_type1_b5(2,:)',all_type1_b5(3,:)',all_type1_b5(4,:)',all_type1_b5(5,:)',...
    'VariableNames',...
    {'IDs','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12','v13','v14','v15',...
    'v16','v17','v18','v19','v20','v21','v22','v23','v24','v25'});

withinDesign = table([ones(1,5), ones(1,5)*2, ones(1,5)*3,ones(1,5)*4, ones(1,5)*5]',[1:5,1:5,1:5,1:5,1:5]','VariableNames',{'Block','Difficulty'});
withinDesign.Block = categorical(withinDesign.Block);
withinDesign.Difficulty = categorical(withinDesign.Difficulty);

rm = fitrm(t, 'v1-v25 ~ 1', 'WithinDesign', withinDesign);
mauchly(rm)
AT = ranova(rm, 'WithinModel', 'Block*Difficulty');  
%for GG corrected df:
tbl = epsilon(rm); eps = tbl.GreenhouseGeisser;
AT.df_gg = AT.DF .* tbl.GreenhouseGeisser
 %%  
ClrMap = colormap('cool');
Clr = [42, 84, 126, 168, 210];
x = 1:5;
figure;
plot(x,mean_t1_b1,'.-','Color',ClrMap(Clr(1),:),'LineWidth',1.5);hold on;
er = errorbar([1 2 3 4 5],mean_t1_b1,[sem1_b1 sem2_b1 sem3_b1 sem4_b1 sem5_b1],[sem1_b1 sem2_b1 sem3_b1 sem4_b1 sem5_b1]);
er.Color = ClrMap(Clr(1),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b2,'.-','Color',ClrMap(Clr(2),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b2,[sem1_b2 sem2_b2 sem3_b2 sem4_b2 sem5_b2],[sem1_b2 sem2_b2 sem3_b2 sem4_b2 sem5_b2]);
er.Color = ClrMap(Clr(2),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b3,'.-','Color',ClrMap(Clr(3),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b3,[sem1_b3 sem2_b3 sem3_b3 sem4_b3 sem5_b3],[sem1_b3 sem2_b3 sem3_b3 sem4_b3 sem5_b3]);
er.Color = ClrMap(Clr(3),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b4,'.-','Color',ClrMap(Clr(4),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b4,[sem1_b4 sem2_b4 sem3_b4 sem4_b4 sem5_b4],[sem1_b4 sem2_b4 sem3_b4 sem4_b4 sem5_b4]);
er.Color = ClrMap(Clr(4),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b5,'.-','Color',ClrMap(Clr(5),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b5,[sem1_b5 sem2_b5 sem3_b5 sem4_b5 sem5_b5],[sem1_b5 sem2_b5 sem3_b5 sem4_b5 sem5_b5]);
er.Color = ClrMap(Clr(5),:);                            
er.LineStyle = 'none';
title('Mean type-2 PF'); ylim([0 1]);
xlim([0.8 5.2]);
% legend('Block 1','','Block 2','','Block 3','','Block 4','','Block 5','','location','southeast');
% legend('box','off');
ylabel('Proportion confident'); xlabel('Evidence discriminability (low to high)');hold off;
xticks([1 2 3 4 5]); xticklabels([0.075 0.15 0.225 0.3 0.375]);

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12)

%inset
% axes('Position',[.55 .17 .3 .4]); box off;
% set(gca,'FontSize',14);
figure;
params = [b1';b2';b3';b4';b5'];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean proportion confident'); 
hold on;
xlabel('Block');
ylim([0 1]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',11)
%% ===============================================================================================================
%                       reaction time (response t1)
% ================================================================================================================
clear all; close all; clc;
lazypath = '\experiment_data\Behaviour\';
cd(lazypath);

Folders = dir(lazypath); %% NAME OF FOLDERS
Folders(1:2) = [];
for i = 1:length(Folders)
    FolderName = Folders(i).name;
    datafolder2 = [lazypath ,'\', FolderName]; %
    cd(datafolder2)% go inside the folder 
    subject = FolderName;
    Ind = FolderName; % here is just to find the name of the subject in case your folders are super long name e.i. find(FolderName == 'S'); %name of your folder 
    subjectname = subject; % in case you need to take some stuff out (1:Ind-4);
    
    cd([lazypath,'\',FolderName]);
    behData = dir('*_1.csv');
    beh = readtable(behData.name);


    block = [ones(180,1); ones(180,1)*2;ones(180,1)*3;ones(180,1)*4;ones(180,1)*5];
    accuracy = beh.key_resp_corr;
    contrast = beh.contrast;
    accuracy(isnan(accuracy))=[];
    contrast(isnan(contrast))=[];
    confidence = beh.key_resp_2_keys; confidence = confidence(~cellfun('isempty',confidence));
    respRT = beh.key_resp_rt; respRT(isnan(respRT)) = [];
    confrespRT = beh.key_resp_2_rt; confrespRT(isnan(confrespRT)) = [];
    for ii = 1:length(confidence)
        if confidence(ii,:) == "k"
            confidence2(ii,:) = 1;
        elseif confidence(ii,:) == "l"
            confidence2(ii,:) = 0;
        end
    end


    %behaviour by difficulty
    if i == 22
       rt = respRT(find(contrast == 0.075 & block == 1),:); %trials 3 and 6 overall are huge outliers so removing
       oneb1(i) = mean(rt([2,4:end],:),1);
       rtb1 = respRT(find( block == 1),:);
       b1(i) = mean(rtb1([1:2,4:5,7:end],:),1);
    else
        oneb1(i) = mean(respRT(find(contrast == 0.075 & block == 1),:));
        b1(i) = mean(respRT(find( block == 1),:));
    end
    oneb2(i) = mean(respRT(find(contrast == 0.075 & block == 2),:));
    oneb3(i) = mean(respRT(find(contrast == 0.075 & block == 3),:));
    oneb4(i) = mean(respRT(find(contrast == 0.075 & block == 4),:));
    oneb5(i) = mean(respRT(find(contrast == 0.075 & block == 5),:));

    twob1(i) = mean(respRT(find(contrast == 0.15 & block == 1),:));
    twob2(i) = mean(respRT(find(contrast == 0.15 & block == 2),:));
    twob3(i) = mean(respRT(find(contrast == 0.15 & block == 3),:));
    twob4(i) = mean(respRT(find(contrast == 0.15 & block == 4),:));
    twob5(i) = mean(respRT(find(contrast == 0.15 & block == 5),:));

    threeb1(i) = mean(respRT(find(contrast == 0.225 & block == 1),:));
    threeb2(i) = mean(respRT(find(contrast == 0.225 & block == 2),:));
    threeb3(i) = mean(respRT(find(contrast == 0.225 & block == 3),:));
    threeb4(i) = mean(respRT(find(contrast == 0.225 & block == 4),:));
    threeb5(i) = mean(respRT(find(contrast == 0.225 & block == 5),:));

    fourb1(i) = mean(respRT(find(contrast == 0.3 & block == 1),:));
    fourb2(i) = mean(respRT(find(contrast == 0.3 & block == 2),:));
    fourb3(i) = mean(respRT(find(contrast == 0.3 & block == 3),:));
    fourb4(i) = mean(respRT(find(contrast == 0.3 & block == 4),:));
    fourb5(i) = mean(respRT(find(contrast == 0.3 & block == 5),:));

    fiveb1(i) = mean(respRT(find(contrast == 0.375 & block == 1),:));
    fiveb2(i) = mean(respRT(find(contrast == 0.375 & block == 2),:));
    fiveb3(i) = mean(respRT(find(contrast == 0.375 & block == 3),:));
    fiveb4(i) = mean(respRT(find(contrast == 0.375 & block == 4),:));
    fiveb5(i) = mean(respRT(find(contrast == 0.375 & block == 5),:));

    
    b2(i) = mean(respRT(find( block == 2),:));
    b3(i) = mean(respRT(find( block == 3),:));
    b4(i) = mean(respRT(find( block == 4),:));
    b5(i) = mean(respRT(find( block == 5),:));

end
    all_type1_b1 = [oneb1;twob1;threeb1;fourb1;fiveb1];
    all_type1_b2 = [oneb2;twob2;threeb2;fourb2;fiveb2];
    all_type1_b3 = [oneb3;twob3;threeb3;fourb3;fiveb3];
    all_type1_b4 = [oneb4;twob4;threeb4;fourb4;fiveb4];
    all_type1_b5 = [oneb5;twob5;threeb5;fourb5;fiveb5];
    
    mean_t1_b1 = mean(all_type1_b1,2);
    mean_t1_b2 = mean(all_type1_b2,2);
    mean_t1_b3 = mean(all_type1_b3,2);
    mean_t1_b4 = mean(all_type1_b4,2);
    mean_t1_b5 = mean(all_type1_b5,2);

    sem1_b1 = std(all_type1_b1(1,:))./sqrt(36);
    sem2_b1 = std(all_type1_b1(2,:))./sqrt(36);
    sem3_b1 = std(all_type1_b1(3,:))./sqrt(36);
    sem4_b1 = std(all_type1_b1(4,:))./sqrt(36);
    sem5_b1 = std(all_type1_b1(5,:))./sqrt(36);
    
    sem1_b2 = std(all_type1_b2(1,:))./sqrt(36);
    sem2_b2 = std(all_type1_b2(2,:))./sqrt(36);
    sem3_b2 = std(all_type1_b2(3,:))./sqrt(36);
    sem4_b2 = std(all_type1_b2(4,:))./sqrt(36);
    sem5_b2 = std(all_type1_b2(5,:))./sqrt(36);
    
    sem1_b3 = std(all_type1_b3(1,:))./sqrt(36);
    sem2_b3 = std(all_type1_b3(2,:))./sqrt(36);
    sem3_b3 = std(all_type1_b3(3,:))./sqrt(36);
    sem4_b3 = std(all_type1_b3(4,:))./sqrt(36);
    sem5_b3 = std(all_type1_b3(5,:))./sqrt(36);

    sem1_b4 = std(all_type1_b4(1,:))./sqrt(36);
    sem2_b4 = std(all_type1_b4(2,:))./sqrt(36);
    sem3_b4 = std(all_type1_b4(3,:))./sqrt(36);
    sem4_b4 = std(all_type1_b4(4,:))./sqrt(36);
    sem5_b4 = std(all_type1_b4(5,:))./sqrt(36);

    sem1_b5 = std(all_type1_b5(1,:))./sqrt(36);
    sem2_b5 = std(all_type1_b5(2,:))./sqrt(36);
    sem3_b5 = std(all_type1_b5(3,:))./sqrt(36);
    sem4_b5 = std(all_type1_b5(4,:))./sqrt(36);
    sem5_b5 = std(all_type1_b5(5,:))./sqrt(36);
  
 % STATS - 5X5 WITHINPARTICIPANT ANOVA with block and difficulty as IV and
% performance as DV
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

t = table(subnums1',all_type1_b1(1,:)',all_type1_b1(2,:)',all_type1_b1(3,:)',all_type1_b1(4,:)',all_type1_b1(5,:)',...
    all_type1_b2(1,:)',all_type1_b2(2,:)',all_type1_b2(3,:)',all_type1_b2(4,:)',all_type1_b2(5,:)',...
    all_type1_b3(1,:)',all_type1_b3(2,:)',all_type1_b3(3,:)',all_type1_b3(4,:)',all_type1_b3(5,:)',...
    all_type1_b4(1,:)',all_type1_b4(2,:)',all_type1_b4(3,:)',all_type1_b4(4,:)',all_type1_b4(5,:)',....
    all_type1_b5(1,:)',all_type1_b5(2,:)',all_type1_b5(3,:)',all_type1_b5(4,:)',all_type1_b5(5,:)',...
    'VariableNames',...
    {'IDs','v1','v2','v3','v4','v5','v6','v7','v8','v9','v10','v11','v12','v13','v14','v15',...
    'v16','v17','v18','v19','v20','v21','v22','v23','v24','v25'});

withinDesign = table([ones(1,5), ones(1,5)*2, ones(1,5)*3,ones(1,5)*4, ones(1,5)*5]',[1:5,1:5,1:5,1:5,1:5]','VariableNames',{'Block','Difficulty'});
withinDesign.Block = categorical(withinDesign.Block);
withinDesign.Difficulty = categorical(withinDesign.Difficulty);

rm = fitrm(t, 'v1-v25 ~ 1', 'WithinDesign', withinDesign);
mauchly(rm)
AT = ranova(rm, 'WithinModel', 'Block*Difficulty');   
%for GG corrected df:
tbl = epsilon(rm); eps = tbl.GreenhouseGeisser;
AT.df_gg = AT.DF .* tbl.GreenhouseGeisser
%%
   
ClrMap = colormap('cool');
Clr = [42, 84, 126, 168, 210];
x = 1:5;
figure;
plot(x,mean_t1_b1,'.-','Color',ClrMap(Clr(1),:),'LineWidth',1.5);hold on;
er = errorbar([1 2 3 4 5],mean_t1_b1,[sem1_b1 sem2_b1 sem3_b1 sem4_b1 sem5_b1],[sem1_b1 sem2_b1 sem3_b1 sem4_b1 sem5_b1]);
er.Color = ClrMap(Clr(1),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b2,'.-','Color',ClrMap(Clr(2),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b2,[sem1_b2 sem2_b2 sem3_b2 sem4_b2 sem5_b2],[sem1_b2 sem2_b2 sem3_b2 sem4_b2 sem5_b2]);
er.Color = ClrMap(Clr(2),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b3,'.-','Color',ClrMap(Clr(3),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b3,[sem1_b3 sem2_b3 sem3_b3 sem4_b3 sem5_b3],[sem1_b3 sem2_b3 sem3_b3 sem4_b3 sem5_b3]);
er.Color = ClrMap(Clr(3),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b4,'.-','Color',ClrMap(Clr(4),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b4,[sem1_b4 sem2_b4 sem3_b4 sem4_b4 sem5_b4],[sem1_b4 sem2_b4 sem3_b4 sem4_b4 sem5_b4]);
er.Color = ClrMap(Clr(4),:);                            
er.LineStyle = 'none';
plot(x,mean_t1_b5,'.-','Color',ClrMap(Clr(5),:),'LineWidth',1.5);
er = errorbar([1 2 3 4 5],mean_t1_b5,[sem1_b5 sem2_b5 sem3_b5 sem4_b5 sem5_b5],[sem1_b5 sem2_b5 sem3_b5 sem4_b5 sem5_b5]);
er.Color = ClrMap(Clr(5),:);                            
er.LineStyle = 'none';
title('Mean response RT PF'); %ylim([0 1]);
xlim([0.8 5.2]);
legend('Block 1','','Block 2','','Block 3','','Block 4','','Block 5','','location','northeast');
legend('box','off');
ylabel('Mean response RT'); xlabel('Evidence discriminability (low to high)');hold off;
xticks([1 2 3 4 5]); xticklabels([0.075 0.15 0.225 0.3 0.375]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12)
%%
figure;
params = [b1';b2';b3';b4';b5'];
group = [ones(36,1); ones(36,1)*2; ones(36,1)*3; ones(36,1)*4; ones(36,1)*5];
params = [group,params];
b = boxplot(params(:,2),group,'Colors',ClrMap(Clr(1:5),:),'Symbol',''); xticklabels({'1','2','3','4','5'});
set(b,'LineWidth',1.1);
ylabel('Mean proportion RT'); 
hold on;
xlabel('Block');
% ylim([0 1]);
violinscatter(1,params(group==1,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(1),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(2,params(group==2,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(2),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(3,params(group==3,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(3),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(4,params(group==4,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(4),:),'MarkerFaceAlpha',0.3) ; 
violinscatter(5,params(group==5,2),0.1,'o','filled','MarkerFaceColor',ClrMap(Clr(5),:),'MarkerFaceAlpha',0.3) ; a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',11)

subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};
%stats (ANOVA?)
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',b1',b2',b3',b4',b5','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
eta = ranovatbl.SumSq(1)/(ranovatbl.SumSq(1)+ranovatbl.SumSq(2))
c = multcompare(rm,'Blocks')
%for GG corrected df:
tbl = epsilon(rm); eps = tbl.GreenhouseGeisser;
ranovatbl.df_gg = ranovatbl.DF .* tbl.GreenhouseGeisser
%% ===============================================================================================================
%                       reaction time (response t1) - breaking down the
%                       interaction
% ================================================================================================================
clear all; close all; clc;
lazypath = '\experiment_data\Behaviour\';
cd(lazypath);

Folders = dir(lazypath); %% NAME OF FOLDERS
Folders(1:2) = [];
for i = 1:length(Folders)
    FolderName = Folders(i).name;
    datafolder2 = [lazypath ,'\', FolderName]; %
    cd(datafolder2)% go inside the folder 
    subject = FolderName;
    Ind = FolderName; % here is just to find the name of the subject in case your folders are super long name e.i. find(FolderName == 'S'); %name of your folder 
    subjectname = subject; % in case you need to take some stuff out (1:Ind-4);
    
    cd([lazypath,'\',FolderName]);
    behData = dir('*_1.csv');
    beh = readtable(behData.name);


    block = [ones(180,1); ones(180,1)*2;ones(180,1)*3;ones(180,1)*4;ones(180,1)*5];
    accuracy = beh.key_resp_corr;
    contrast = beh.contrast;
    accuracy(isnan(accuracy))=[];
    contrast(isnan(contrast))=[];
    confidence = beh.key_resp_2_keys; confidence = confidence(~cellfun('isempty',confidence));
    respRT = beh.key_resp_rt; respRT(isnan(respRT)) = [];
    confrespRT = beh.key_resp_2_rt; confrespRT(isnan(confrespRT)) = [];
    for ii = 1:length(confidence)
        if confidence(ii,:) == "k"
            confidence2(ii,:) = 1;
        elseif confidence(ii,:) == "l"
            confidence2(ii,:) = 0;
        end
    end


    %behaviour by difficulty
    if i == 22
       rt = respRT(find(contrast == 0.075 & block == 1),:); %trials 3 and 6 overall are huge outliers so removing
       oneb1(i) = mean(rt([2,4:end],:),1);
    else
    oneb1(i) = mean(respRT(find(contrast == 0.075 & block == 1),:));
    end
    oneb2(i) = mean(respRT(find(contrast == 0.075 & block == 2),:));
    oneb3(i) = mean(respRT(find(contrast == 0.075 & block == 3),:));
    oneb4(i) = mean(respRT(find(contrast == 0.075 & block == 4),:));
    oneb5(i) = mean(respRT(find(contrast == 0.075 & block == 5),:));

    twob1(i) = mean(respRT(find(contrast == 0.15 & block == 1),:));
    twob2(i) = mean(respRT(find(contrast == 0.15 & block == 2),:));
    twob3(i) = mean(respRT(find(contrast == 0.15 & block == 3),:));
    twob4(i) = mean(respRT(find(contrast == 0.15 & block == 4),:));
    twob5(i) = mean(respRT(find(contrast == 0.15 & block == 5),:));

    threeb1(i) = mean(respRT(find(contrast == 0.225 & block == 1),:));
    threeb2(i) = mean(respRT(find(contrast == 0.225 & block == 2),:));
    threeb3(i) = mean(respRT(find(contrast == 0.225 & block == 3),:));
    threeb4(i) = mean(respRT(find(contrast == 0.225 & block == 4),:));
    threeb5(i) = mean(respRT(find(contrast == 0.225 & block == 5),:));

    fourb1(i) = mean(respRT(find(contrast == 0.3 & block == 1),:));
    fourb2(i) = mean(respRT(find(contrast == 0.3 & block == 2),:));
    fourb3(i) = mean(respRT(find(contrast == 0.3 & block == 3),:));
    fourb4(i) = mean(respRT(find(contrast == 0.3 & block == 4),:));
    fourb5(i) = mean(respRT(find(contrast == 0.3 & block == 5),:));

    fiveb1(i) = mean(respRT(find(contrast == 0.375 & block == 1),:));
    fiveb2(i) = mean(respRT(find(contrast == 0.375 & block == 2),:));
    fiveb3(i) = mean(respRT(find(contrast == 0.375 & block == 3),:));
    fiveb4(i) = mean(respRT(find(contrast == 0.375 & block == 4),:));
    fiveb5(i) = mean(respRT(find(contrast == 0.375 & block == 5),:));

    

end
%this one tests across block at each difficulty
    all_type1_b1 = [oneb1;oneb2;oneb3;oneb4;oneb5];
    all_type1_b2 = [twob1;twob2;twob3;twob4;twob5];
    all_type1_b3 = [threeb1;threeb2;threeb3;threeb4;threeb5];
    all_type1_b4 = [fourb1;fourb2;fourb3;fourb4;fourb5];
    all_type1_b5 = [fiveb1;fiveb2;fiveb3;fiveb4;fiveb5];

    %this one tests across difficulty at each block (is the sig one)
%     all_type1_b1 = [oneb1;twob1;threeb1;fourb1;fiveb1];
%     all_type1_b2 = [oneb2;twob2;threeb2;fourb2;fiveb2];
%     all_type1_b3 = [oneb3;twob3;threeb3;fourb3;fiveb3];
%     all_type1_b4 = [oneb4;twob4;threeb4;fourb4;fiveb4];
%     all_type1_b5 = [oneb5;twob5;threeb5;fourb5;fiveb5];

%     mean_t1_b1 = mean(all_type1_b1,2);
%     mean_t1_b2 = mean(all_type1_b2,2);
%     mean_t1_b3 = mean(all_type1_b3,2);
%     mean_t1_b4 = mean(all_type1_b4,2);
%     mean_t1_b5 = mean(all_type1_b5,2);


  
 % STATS - 5X5 WITHINPARTICIPANT ANOVA with block and difficulty as IV and
% performance as DV
subnums1 = {'155675','157061','192615','201786','205264','207302','215501','215990',...
    '222722','275945','336764','356885','388225','401502','403438','407302','414442',...
    '453240','472127','508708','517841','528695', '546350', '553536','556277', '571278','629010','632182',...
    '647517','652003','693339','711785','751325','762201','792475','804448'};

%stats (ANOVA?)
%no 1
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',all_type1_b1(1,:)',all_type1_b1(2,:)',all_type1_b1(3,:)',all_type1_b1(4,:)',all_type1_b1(5,:)','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
mauchly(rm)
c = multcompare(rm,'Blocks')

%no 2
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',all_type1_b2(1,:)',all_type1_b2(2,:)',all_type1_b2(3,:)',all_type1_b2(4,:)',all_type1_b2(5,:)','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
mauchly(rm)
c = multcompare(rm,'Blocks')

%no 3
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',all_type1_b3(1,:)',all_type1_b3(2,:)',all_type1_b3(3,:)',all_type1_b3(4,:)',all_type1_b3(5,:)','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
mauchly(rm)
c = multcompare(rm,'Blocks')

%no 4
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',all_type1_b4(1,:)',all_type1_b4(2,:)',all_type1_b4(3,:)',all_type1_b4(4,:)',all_type1_b4(5,:)','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
mauchly(rm)
c = multcompare(rm,'Blocks')

%no 5
%repeated measures ANOVA for CF/PW across blocks
t = table(subnums1',all_type1_b5(1,:)',all_type1_b5(2,:)',all_type1_b5(3,:)',all_type1_b5(4,:)',all_type1_b5(5,:)','VariableNames',...
    {'IDs','B1','B2','B3','B4','B5'});
Blocks = table(categorical(1:5)','VariableNames',{'Blocks'});
rm = fitrm(t, 'B1-B5 ~ 1','WithinDesign',Blocks);
ranovatbl = ranova(rm)
mauchly(rm)
c = multcompare(rm,'Blocks')