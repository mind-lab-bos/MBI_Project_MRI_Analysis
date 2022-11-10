%%%%Updated Marsbar Pipeline that includes ROI extraction of YA and OA Data%%%%

%defining and loading rois
cd '/home/mquinci/Wang_etal_ROIs'
roi_list = dir('*.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Extracting Young Adult Data%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define subjects
cd '/scratch/mquinci/mci_spm_musbid/pre';
sublistYA = struct('path',{'/scratch/mquinci/mci_spm_musbid/pre/sub-042-AZHA1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-043-SLEE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-044-LGAR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-045-CZAL1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-046-APOP1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-048-LWIE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-049-GWAI1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-050-SCAM1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-051-KCHU1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-052-KPLU1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-053-KWAL1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-054-ZCAR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-055-JLAW1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-056-ASUR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-057-IZIK1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-058-SAPE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-059-TBOY1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-060-MCHI1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-061-JCLA1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-062-SSOH1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-063-AWIL1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-064-STHO1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-065-MBIC1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-066-CKOB1/mus-bid_spm/liking'});
%sublistYA = struct('path',{'/scratch/mquinci/mci_spm_musbid/pre/sub-057-IZIK1/mus-bid_spm/liking'});

%define empty array
youngAdultdata_24 = cell(1, numel(sublistYA));

%loop through all participants and ROIs
for subs = 1:numel(sublistYA)
% Make marsbar design object
D = mardo(strcat(sublistYA(subs).path, '/', 'SPM.mat'));
% Set fmristat AR modelling
D = autocorr(D, 'fmristat', 2);

    for rois = 1:numel(roi_list)
        % Make marsbar ROI object
        R(rois, subs) = maroi(strcat(roi_list(rois).folder,'/',roi_list(rois).name));
        % Fetch data into marsbar data object
        Y  = get_marsy(R(rois), D, 'mean');
        %Summarize data
        youngAdultdata_24{subs}(:, rois) = summary_data(marsy(Y));
    end
end

%save data to a .mat file
cd '/home/mquinci/';
save('youngAdultdata_Wang_updated24.mat','youngAdultdata_updated');


%load onsets
cd '/home/mquinci/roiroi_scripts'
load('onsets.txt');

%isolate the 42 trs for each trial that correspond with the music listening
updatedDataYA = cell(1, numel(sublistYA));
for i = 1:size(youngAdultdata_24,2)
    
   for j =  1:length(onsets)

       updatedDataYA{i}((42*(j-1)+1):42*j, :) = youngAdultdata_24{i}(onsets(j):onsets(j) + 41, :);
       
   end
end

%average across each stim (producing 24 rows, one for each stim)
meanDataAllYA = cell(1,numel(sublistYA));
for i = 1:size(updatedDataYA, 2)
    
    for j= 1:length(onsets)
        
        meanDataAllYA{i}(j, :) = mean(updatedDataYA{i}((42*(j-1)+1):42*j, :));
        
    end
end


%load stim order
cd '/home/mquinci/ROI_Analyses/Liking/Young_Adults/stim_order'
StimListYA = dir('*.txt');
stim_orderYA  = cell(1, numel(StimListYA));
for stim = 1:numel(StimListYA)
  FileData     = load(StimListYA(stim).name);
  stim_orderYA{stim} = FileData;
  stim_orderYA{stim}(1,:)= [];
  stim_orderliking{stim} = stim_orderYA{stim}(:,2);
end
cd ../

liking_full_order = cell(1,numel(StimListYA));
hate_YA = cell(1,numel(sublistYA));
neutral_YA = cell(1, numel(sublistYA));
like_YA = cell(1, numel(sublistYA));
love_YA = cell(1, numel(sublistYA));
cd '/home/mquinci/roiroi_scripts/'
load('stimstarts.txt');
for sub = 1:size(updatedDataYA, 2)
    
    hate = 1;
    neutral = 1;
    like = 1;
    love = 1;
    
    for stim = 1:24
        
        if stim_orderliking{sub}(stim) == 1
            hate_YA{sub}((hate-1)*42+1:hate*42,:) = updatedDataYA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            hate=hate+1;
            
        elseif stim_orderliking{sub}(stim) == 2
            neutral_YA{sub}((neutral-1)*42+1:neutral*42,:) = updatedDataYA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            neutral=neutral+1;
            
        elseif stim_orderliking{sub}(stim) == 3
            like_YA{sub}((like-1)*42+1:like*42,:) = updatedDataYA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            like=like+1;
            
        elseif stim_orderliking{sub}(stim) == 4
            love_YA{sub}((love-1)*42+1:love*42,:) = updatedDataYA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            love=love+1;
            
        else 
            
        end
        
    end
end


for sub = 1:length(sublistYA)
 if numel(hate_YA{sub}) > 0   
    RhateYA(:,:,sub) = corr(hate_YA{sub});
   else
     RhateYA(:,:,sub) = NaN;
 end
   
 if numel(neutral_YA{sub}) > 0   
    RneutralYA(:,:,sub) = corr(neutral_YA{sub});
   else
     RneutralYA(:,:,sub) = NaN;
 end
 
 if numel(like_YA{sub}) > 0   
    RlikeYA(:,:,sub) = corr(like_YA{sub});
   else
     RlikeYA(:,:,sub) = NaN;
   end
   
 if numel(love_YA{sub}) > 0   
    RloveYA(:,:,sub) = corr(love_YA{sub});
   else
     RloveYA(:,:,sub) = NaN;
   end
 
 
end

%Define structs
Corr_Mats_YA.Hate.All = RhateYA(:,:,:);
Corr_Mats_YA.Neutral.All = RneutralYA(:,:,:);
Corr_Mats_YA.Like.All = RlikeYA(:,:,:);
Corr_Mats_YA.Love.All = RloveYA(:,:,:);


%xvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr', 'MPFC', 'Precuneus'};
%yvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr', 'MPFC', 'Precuneus'};


xvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr'};
yvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr'};


meanRhateYA = nanmean(RhateYA,3);
figure;
hateh = heatmap(xvalues, yvalues, meanRhateYA([1:36],[1:36]), 'ColorLimits',[0 1],'GridVisible','off');
colormap(parula)

meanRneutralYA = nanmean(RneutralYA,3);
figure;
neutralh = heatmap(xvalues, yvalues, meanRneutralYA([1:36],[1:36]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRlikeYA = nanmean(RlikeYA,3);
figure;
likeh = heatmap(xvalues, yvalues, meanRlikeYA([1:36],[1:36]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRloveYA= nanmean(RloveYA,3);
figure;
heatmap(xvalues, yvalues, meanRloveYA([1:36],[1:36]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

%Isolate the auditory-auditory section of the matrix (YA)
Corr_Mats_YA.Hate.Aud.Fullmat = Corr_Mats_YA.Hate.All([1:18],[1:18],:);
Corr_Mats_YA.Neutral.Aud.Fullmat = Corr_Mats_YA.Neutral.All([1:18],[1:18],:);
Corr_Mats_YA.Like.Aud.Fullmat = Corr_Mats_YA.Like.All([1:18],[1:18],:);
Corr_Mats_YA.Love.Aud.Fullmat = Corr_Mats_YA.Love.All([1:18],[1:18],:);

%Isolate the top half of the corr matrix for each liking rating
%(aud-aud)(YA)
for sub = 1:size(Corr_Mats_YA.Hate.All,3) 
    
    temp = Corr_Mats_YA.Hate.Aud.Fullmat(:,:,sub);
    Corr_Mats_YA.Hate.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_YA.Neutral.Aud.Fullmat(:,:,sub);
    Corr_Mats_YA.Neutral.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_YA.Like.Aud.Fullmat(:,:,sub);
    Corr_Mats_YA.Like.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_YA.Love.Aud.Fullmat(:,:,sub);
    Corr_Mats_YA.Love.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
end

Corr_Mats_YA.Hate.Aud.Mean = mean(Corr_Mats_YA.Hate.Aud.List);
Corr_Mats_YA.Neutral.Aud.Mean = mean(Corr_Mats_YA.Neutral.Aud.List);
Corr_Mats_YA.Like.Aud.Mean = mean(Corr_Mats_YA.Like.Aud.List);
Corr_Mats_YA.Love.Aud.Mean = mean(Corr_Mats_YA.Love.Aud.List);

%Isolate the Reward-Reward section of the matrix (YA)
Corr_Mats_YA.Hate.Rew.Fullmat = Corr_Mats_YA.Hate.All([19:36],[19:36],:);
Corr_Mats_YA.Neutral.Rew.Fullmat = Corr_Mats_YA.Neutral.All([19:36],[19:36],:);
Corr_Mats_YA.Like.Rew.Fullmat = Corr_Mats_YA.Like.All([19:36],[19:36],:);
Corr_Mats_YA.Love.Rew.Fullmat = Corr_Mats_YA.Love.All([19:36],[19:36],:);

%Isolate the top half of the corr matrix for each liking rating
%(Rew-Rew)(YA)
for sub = 1:size(Corr_Mats_YA.Hate.All,3) 
    
    temp = Corr_Mats_YA.Hate.Rew.Fullmat(:,:,sub);
    Corr_Mats_YA.Hate.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_YA.Neutral.Rew.Fullmat(:,:,sub);
    Corr_Mats_YA.Neutral.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_YA.Like.Rew.Fullmat(:,:,sub);
    Corr_Mats_YA.Like.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_YA.Love.Rew.Fullmat(:,:,sub);
    Corr_Mats_YA.Love.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
end

Corr_Mats_YA.Hate.Rew.Mean = mean(Corr_Mats_YA.Hate.Rew.List);
Corr_Mats_YA.Neutral.Rew.Mean = mean(Corr_Mats_YA.Neutral.Rew.List);
Corr_Mats_YA.Like.Rew.Mean = mean(Corr_Mats_YA.Like.Rew.List);
Corr_Mats_YA.Love.Rew.Mean = mean(Corr_Mats_YA.Love.Rew.List);


%Isolate the auditory-reward section of the matrix (YA)
Corr_Mats_YA.Hate.AudRew.Fullmat = Corr_Mats_YA.Hate.All([1:18],[19:36],:);
Corr_Mats_YA.Neutral.AudRew.Fullmat = Corr_Mats_YA.Neutral.All([1:18],[19:36],:);
Corr_Mats_YA.Like.AudRew.Fullmat = Corr_Mats_YA.Like.All([1:18],[19:36],:);
Corr_Mats_YA.Love.AudRew.Fullmat = Corr_Mats_YA.Love.All([1:18],[19:36],:);

%Get the mean for each participant across the matrix (aud-rew) (YA)
for sub = 1:size(Corr_Mats_YA.Hate.All,3) 
    
    temp = Corr_Mats_YA.Hate.AudRew.Fullmat(:,:,sub);
    Corr_Mats_YA.Hate.AudRew.List(:,sub) = temp(:);
    
    temp = Corr_Mats_YA.Neutral.AudRew.Fullmat(:,:,sub);
    Corr_Mats_YA.Neutral.AudRew.List(:,sub) = temp(:);
    
    temp = Corr_Mats_YA.Like.AudRew.Fullmat(:,:,sub);
    Corr_Mats_YA.Like.AudRew.List(:,sub) = temp(:);
    
    temp = Corr_Mats_YA.Love.AudRew.Fullmat(:,:,sub);
    Corr_Mats_YA.Love.AudRew.List(:,sub) = temp(:);
end
Corr_Mats_YA.Hate.AudRew.Mean = mean(Corr_Mats_YA.Hate.AudRew.List);
Corr_Mats_YA.Neutral.AudRew.Mean = mean(Corr_Mats_YA.Neutral.AudRew.List);
Corr_Mats_YA.Like.AudRew.Mean = mean(Corr_Mats_YA.Like.AudRew.List);
Corr_Mats_YA.Love.AudRew.Mean = mean(Corr_Mats_YA.Love.AudRew.List);

%Isolate mPFC-auditory connections (YA)
Corr_Mats_YA.Hate.AudmPFC.Fullmat = Corr_Mats_YA.Hate.All(37,1:18,:);
Corr_Mats_YA.Neutral.AudmPFC.Fullmat = Corr_Mats_YA.Neutral.All(37,1:18,:);
Corr_Mats_YA.Like.AudmPFC.Fullmat = Corr_Mats_YA.Like.All(37,1:18,:);
Corr_Mats_YA.Love.AudmPFC.Fullmat = Corr_Mats_YA.Love.All(37,1:18,:);

%Isolate Precun-auditory connections (YA)
Corr_Mats_YA.Hate.AudPrecun.Fullmat = Corr_Mats_YA.Hate.All(38,1:18,:);
Corr_Mats_YA.Neutral.AudPrecun.Fullmat = Corr_Mats_YA.Neutral.All(38,1:18,:);
Corr_Mats_YA.Like.AudPrecun.Fullmat = Corr_Mats_YA.Like.All(38,1:18,:);
Corr_Mats_YA.Love.AudPrecun.Fullmat = Corr_Mats_YA.Love.All(38,1:18,:);

%Isolate mPFC-Reward connections (YA)
Corr_Mats_YA.Hate.RewmPFC.Fullmat = Corr_Mats_YA.Hate.All(37,19:36,:);
Corr_Mats_YA.Neutral.RewmPFC.Fullmat = Corr_Mats_YA.Neutral.All(37,19:36,:);
Corr_Mats_YA.Like.RewmPFC.Fullmat = Corr_Mats_YA.Like.All(37,19:36,:);
Corr_Mats_YA.Love.RewmPFC.Fullmat = Corr_Mats_YA.Love.All(37,19:36,:);

%Isolate Precun-Reward connections (YA)
Corr_Mats_YA.Hate.RewPrecun.Fullmat = Corr_Mats_YA.Hate.All(38,19:36,:);
Corr_Mats_YA.Neutral.RewPrecun.Fullmat = Corr_Mats_YA.Neutral.All(38,19:36,:);
Corr_Mats_YA.Like.RewPrecun.Fullmat = Corr_Mats_YA.Like.All(38,19:36,:);
Corr_Mats_YA.Love.RewPrecun.Fullmat = Corr_Mats_YA.Love.All(38,19:36,:);

cd '/home/mquinci/';
save('Corr_Mats_YA.mat','Corr_Mats_YA');



















%combine the three files together and get the mean (one mean for hate
%aud-aud, one mean for hate rew-rew, and one mean for hate aud-rew. These
%would all go in separate plots. Dont forget about error bars

%Auditory-Auditory
for sub = 1:length(sublistYA) 
auditory_auditoryYA(sub,:) = [rdatahateaud_meanYA(:,sub), rdataneutralaud_meanYA(:,sub), rdatalikeaud_meanYA(:,sub), rdataloveaud_meanYA(:,sub)];
end

SE_audYA = nanstd(auditory_auditoryYA)/sqrt(length(auditory_auditoryYA));

auditory_auditorymeanYA = nanmean(auditory_auditoryYA);
X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
figure;  
bar(X,auditory_auditorymeanYA, 'FaceColor', [0.6350 0.0780 0.1840])
title('Young Adult Connectivity between Auditory-Auditory Regions')

hold on

er = errorbar(X,auditory_auditorymeanYA,SE_audYA);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%Reward-Reward
for sub = 1:length(sublistYA) 
reward_rewardYA(sub,:) = [rdatahaterew_meanYA(:,sub), rdataneutralrew_meanYA(:,sub), rdatalikerew_meanYA(:,sub), rdataloverew_meanYA(:,sub)];
end

SE_rewYA = nanstd(reward_rewardYA)/sqrt(length(reward_rewardYA));

reward_rewardmeanYA = nanmean(reward_rewardYA);
X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
figure;  
bar(X,reward_rewardmeanYA, 'FaceColor', [0.6350 0.0780 0.1840])
title('Young Adult Connectivity between Reward-Reward Regions')

hold on

er = errorbar(X,reward_rewardmeanYA,SE_rewYA);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%Auditory-Reward
for sub = 1:length(sublistYA) 
auditory_rewardYA(sub,:) = [rdatahateaudrew_meanYA(:,sub), rdataneutralaudrew_meanYA(:,sub), rdatalikeaudrew_meanYA(:,sub), rdataloveaudrew_meanYA(:,sub)];
end

SE_audrewYA = nanstd(auditory_rewardYA)/sqrt(length(auditory_rewardYA));

auditory_rewardmeanYA = nanmean(auditory_rewardYA);
X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
figure;  
bar(X,auditory_rewardmeanYA, 'FaceColor', [0.6350 0.0780 0.1840])
title('Young Adult Connectivity between Auditory-Reward Regions')

hold on

er = errorbar(X,auditory_rewardmeanYA,SE_audrewYA);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Extracting Older Adult Data%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define subjects
cd '/scratch/mquinci/mci_spm_musbid/pre';
sublistOA = struct('path',{'/scratch/mquinci/mci_spm_musbid/pre/sub-01-FKEE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-03-JSCH1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-05-LGEN1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-06-GHER1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-08-LDIR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-010-JPRE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-011-CWEI1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-012-MBLO1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-015-EROS1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-016-TCIE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-018-GCAS1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-019-DFAI1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-021-JKAY1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-022-DCAT1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-023-LOLE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-028-SCAE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-029-EHER1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-030-LJAC1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-032-ASHA1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-036-GSHA1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-037-AEHR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-039-ELAR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-040-BFOR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-047-TSTA1/mus-bid_spm/liking'});

%define empty array
olderAdultdata_24 = cell(1, numel(sublistOA));

%loop through all participants and ROIs
for subs = 1:numel(sublistOA)
% Make marsbar design object
D = mardo(strcat(sublistOA(subs).path, '/', 'SPM.mat'));
% Set fmristat AR modelling
D = autocorr(D, 'fmristat', 2);

    for rois = 1:numel(roi_list)
        % Make marsbar ROI object
        R(rois, subs) = maroi(strcat(roi_list(rois).folder,'/',roi_list(rois).name));
        % Fetch data into marsbar data object
        Y  = get_marsy(R(rois), D, 'mean');
        %Summarize data
        olderAdultdata_24{subs}(:, rois) = summary_data(marsy(Y));
    end
end

%save data to a .mat file
cd '/home/mquinci/';
save('olderAdultdata_updated24_10_6_22.mat','olderAdultdata_24');

%load onsets
cd '/scratch/mquinci/roiroi_scripts'
load('onsets.txt');

%isolate the 42 trs for each trial that correspond with the music listening
updatedDataOA = cell(1, numel(sublistOA));
for i = 1:size(olderAdultdata_24,2)
    
   for j =  1:length(onsets)

       updatedDataOA{i}((42*(j-1)+1):42*j, :) = olderAdultdata_24{i}(onsets(j):onsets(j) + 41, :);
       
   end
end

%average across each stim (producing 24 rows, one for each stim)
meanDataAllOA = cell(1,numel(sublistOA));
for i = 1:size(updatedDataOA, 2)
    
    for j= 1:length(onsets)
        
        meanDataAllOA{i}(j, :) = mean(updatedDataOA{i}((42*(j-1)+1):42*j, :));
        
    end
end

%load stim order
cd '/home/mquinci/ROI_Analyses/ROI_pipeline_OA/stim_order'
StimListOA = dir('*.txt');
stim_orderOA  = cell(1, numel(StimListOA));
for stim = 1:numel(StimListOA)
  FileData     = load(StimListOA(stim).name);
  stim_orderOA{stim} = FileData;
  stim_orderOA{stim}(1,:)= [];
  stim_orderlikingOA{stim} = stim_orderOA{stim}(:,2);
end
cd ../

%%%%%%***********************%%%%%%%
liking_full_order = cell(1,numel(StimListOA));
hate_OA = cell(1,numel(sublistOA));
neutral_OA = cell(1, numel(sublistOA));
like_OA = cell(1, numel(sublistOA));
love_OA = cell(1, numel(sublistOA));
cd '/scratch/mquinci/roiroi_scripts/'
load('stimstarts.txt');
for sub = 1:size(updatedDataOA, 2)
    
    hate = 1;
    neutral = 1;
    like = 1;
    love = 1;
    
    for stim = 1:24
        
        if stim_orderlikingOA{sub}(stim) == 1
            hate_OA{sub}((hate-1)*42+1:hate*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            hate=hate+1;
            
        elseif stim_orderlikingOA{sub}(stim) == 2
            neutral_OA{sub}((neutral-1)*42+1:neutral*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            neutral=neutral+1;
            
        elseif stim_orderlikingOA{sub}(stim) == 3
            like_OA{sub}((like-1)*42+1:like*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            like=like+1;
            
        elseif stim_orderlikingOA{sub}(stim) == 4
            love_OA{sub}((love-1)*42+1:love*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            love=love+1;
            
        else 
            
        end
        
    end
end


for sub = 1:length(sublistOA)
 if numel(hate_OA{sub}) > 0   
    RhateOA(:,:,sub) = corr(hate_OA{sub});
   else
     RhateOA(:,:,sub) = NaN;
 end
   
 if numel(neutral_OA{sub}) > 0   
    RneutralOA(:,:,sub) = corr(neutral_OA{sub});
   else
     RneutralOA(:,:,sub) = NaN;
 end
 
 if numel(like_OA{sub}) > 0   
    RlikeOA(:,:,sub) = corr(like_OA{sub});
   else
     RlikeOA(:,:,sub) = NaN;
   end
   
 if numel(love_OA{sub}) > 0   
    RloveOA(:,:,sub) = corr(love_OA{sub});
   else
     RloveOA(:,:,sub) = NaN;
   end
 
 
end

Corr_Mats_OA.Hate.All = RhateOA;
Corr_Mats_OA.Neutral.All = RneutralOA;
Corr_Mats_OA.Like.All = RlikeOA;
Corr_Mats_OA.Love.All = RloveOA;

%Isolate the auditory-auditory section of the matrix (OA)
Corr_Mats_OA.Hate.Aud.Fullmat = Corr_Mats_OA.Hate.All([1:18],[1:18],:);
Corr_Mats_OA.Neutral.Aud.Fullmat = Corr_Mats_OA.Neutral.All([1:18],[1:18],:);
Corr_Mats_OA.Like.Aud.Fullmat = Corr_Mats_OA.Like.All([1:18],[1:18],:);
Corr_Mats_OA.Love.Aud.Fullmat = Corr_Mats_OA.Love.All([1:18],[1:18],:);

%Isolate the top half of the corr matrix for each liking rating
%(aud-aud)(OA)
for sub = 1:size(Corr_Mats_OA.Hate.All,3) 
    
    temp = Corr_Mats_OA.Hate.Aud.Fullmat(:,:,sub);
    Corr_Mats_OA.Hate.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_OA.Neutral.Aud.Fullmat(:,:,sub);
    Corr_Mats_OA.Neutral.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_OA.Like.Aud.Fullmat(:,:,sub);
    Corr_Mats_OA.Like.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_OA.Love.Aud.Fullmat(:,:,sub);
    Corr_Mats_OA.Love.Aud.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
end

Corr_Mats_OA.Hate.Aud.Mean = mean(Corr_Mats_OA.Hate.Aud.List);
Corr_Mats_OA.Neutral.Aud.Mean = mean(Corr_Mats_OA.Neutral.Aud.List);
Corr_Mats_OA.Like.Aud.Mean = mean(Corr_Mats_OA.Like.Aud.List);
Corr_Mats_OA.Love.Aud.Mean = mean(Corr_Mats_OA.Love.Aud.List);


%Isolate the reward-reward section of the matrix (OA)
Corr_Mats_OA.Hate.Rew.Fullmat = Corr_Mats_OA.Hate.All([19:36],[19:36],:);
Corr_Mats_OA.Neutral.Rew.Fullmat = Corr_Mats_OA.Neutral.All([19:36],[19:36],:);
Corr_Mats_OA.Like.Rew.Fullmat = Corr_Mats_OA.Like.All([19:36],[19:36],:);
Corr_Mats_OA.Love.Rew.Fullmat = Corr_Mats_OA.Love.All([19:36],[19:36],:);

%Isolate the top half of the corr matrix for each liking rating
%(Rew-Rew)(OA)
for sub = 1:size(Corr_Mats_OA.Hate.All,3) 
    
    temp = Corr_Mats_OA.Hate.Rew.Fullmat(:,:,sub);
    Corr_Mats_OA.Hate.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_OA.Neutral.Rew.Fullmat(:,:,sub);
    Corr_Mats_OA.Neutral.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_OA.Like.Rew.Fullmat(:,:,sub);
    Corr_Mats_OA.Like.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
    temp = Corr_Mats_OA.Love.Rew.Fullmat(:,:,sub);
    Corr_Mats_OA.Love.Rew.List(:,sub) = temp(find(~tril(ones(size(temp)))));
    
end

Corr_Mats_OA.Hate.Rew.Mean = mean(Corr_Mats_OA.Hate.Rew.List);
Corr_Mats_OA.Neutral.Rew.Mean = mean(Corr_Mats_OA.Neutral.Rew.List);
Corr_Mats_OA.Like.Rew.Mean = mean(Corr_Mats_OA.Like.Rew.List);
Corr_Mats_OA.Love.Rew.Mean = mean(Corr_Mats_OA.Love.Rew.List);


%Isolate the auditory-reward section of the matrix (OA)
Corr_Mats_OA.Hate.AudRew.Fullmat = Corr_Mats_OA.Hate.All([1:18],[19:36],:);
Corr_Mats_OA.Neutral.AudRew.Fullmat = Corr_Mats_OA.Neutral.All([1:18],[19:36],:);
Corr_Mats_OA.Like.AudRew.Fullmat = Corr_Mats_OA.Like.All([1:18],[19:36],:);
Corr_Mats_OA.Love.AudRew.Fullmat = Corr_Mats_OA.Love.All([1:18],[19:36],:);

%Get the mean for each participant across the matrix (aud-rew) (OA)
for sub = 1:size(Corr_Mats_OA.Hate.All,3) 
    
    temp = Corr_Mats_OA.Hate.AudRew.Fullmat(:,:,sub);
    Corr_Mats_OA.Hate.AudRew.List(:,sub) = temp(:);
    
    temp = Corr_Mats_OA.Neutral.AudRew.Fullmat(:,:,sub);
    Corr_Mats_OA.Neutral.AudRew.List(:,sub) = temp(:);
    
    temp = Corr_Mats_OA.Like.AudRew.Fullmat(:,:,sub);
    Corr_Mats_OA.Like.AudRew.List(:,sub) = temp(:);
    
    temp = Corr_Mats_OA.Love.AudRew.Fullmat(:,:,sub);
    Corr_Mats_OA.Love.AudRew.List(:,sub) = temp(:);
end
Corr_Mats_OA.Hate.AudRew.Mean = mean(Corr_Mats_OA.Hate.AudRew.List);
Corr_Mats_OA.Neutral.AudRew.Mean = mean(Corr_Mats_OA.Neutral.AudRew.List);
Corr_Mats_OA.Like.AudRew.Mean = mean(Corr_Mats_OA.Like.AudRew.List);
Corr_Mats_OA.Love.AudRew.Mean = mean(Corr_Mats_OA.Love.AudRew.List);


%Isolate mPFC-auditory connections (OA)
Corr_Mats_OA.Hate.AudmPFC.Fullmat = Corr_Mats_OA.Hate.All(37,1:18,:);
Corr_Mats_OA.Neutral.AudmPFC.Fullmat = Corr_Mats_OA.Neutral.All(37,1:18,:);
Corr_Mats_OA.Like.AudmPFC.Fullmat = Corr_Mats_OA.Like.All(37,1:18,:);
Corr_Mats_OA.Love.AudmPFC.Fullmat = Corr_Mats_OA.Love.All(37,1:18,:);

%Isolate Precun-auditory connections (OA)
Corr_Mats_OA.Hate.AudPrecun.Fullmat = Corr_Mats_OA.Hate.All(38,1:18,:);
Corr_Mats_OA.Neutral.AudPrecun.Fullmat = Corr_Mats_OA.Neutral.All(38,1:18,:);
Corr_Mats_OA.Like.AudPrecun.Fullmat = Corr_Mats_OA.Like.All(38,1:18,:);
Corr_Mats_OA.Love.AudPrecun.Fullmat = Corr_Mats_OA.Love.All(38,1:18,:);

%Isolate mPFC-Reward connections (OA)
Corr_Mats_OA.Hate.RewmPFC.Fullmat = Corr_Mats_OA.Hate.All(37,19:36,:);
Corr_Mats_OA.Neutral.RewmPFC.Fullmat = Corr_Mats_OA.Neutral.All(37,19:36,:);
Corr_Mats_OA.Like.RewmPFC.Fullmat = Corr_Mats_OA.Like.All(37,19:36,:);
Corr_Mats_OA.Love.RewmPFC.Fullmat = Corr_Mats_OA.Love.All(37,19:36,:);

%Isolate Precun-Reward connections (OA)
Corr_Mats_OA.Hate.RewPrecun.Fullmat = Corr_Mats_OA.Hate.All(38,19:36,:);
Corr_Mats_OA.Neutral.RewPrecun.Fullmat = Corr_Mats_OA.Neutral.All(38,19:36,:);
Corr_Mats_OA.Like.RewPrecun.Fullmat = Corr_Mats_OA.Like.All(38,19:36,:);
Corr_Mats_OA.Love.RewPrecun.Fullmat = Corr_Mats_OA.Love.All(38,19:36,:);


cd '/home/mquinci/';
save('Corr_Mats_OA.mat','Corr_Mats_OA');


meanRhateOA = nanmean(RhateOA,3);
figure;
heatmap(xvalues, yvalues, meanRhateOA([1:36],[1:36]), 'ColorLimits',[0 1],'GridVisible','off');
colormap(parula)

meanRneutralOA = nanmean(RneutralOA,3);
figure;
heatmap(xvalues, yvalues, meanRneutralOA([1:36],[1:36]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRlikeYA = nanmean(RlikeOA,3);
figure;
heatmap(xvalues, yvalues, meanRlikeYA([1:36],[1:36]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRloveOA= nanmean(RloveOA,3);
figure;
heatmap(xvalues, yvalues, meanRloveOA([1:36],[1:36]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)
