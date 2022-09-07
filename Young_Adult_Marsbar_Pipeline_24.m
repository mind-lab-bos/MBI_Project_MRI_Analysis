%Young Adults

%defining and loading rois
cd '/scratch/mquinci/Wang_ROIs'
roi_list = dir('*.mat');

%define subjects
cd '/scratch/mquinci/mci_spm_musbid/pre';
sublistYA = struct('path',{'/scratch/mquinci/mci_spm_musbid/pre/sub-042-AZHA1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-043-SLEE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-044-LGAR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-045-CZAL1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-046-APOP1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-048-LWIE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-049-GWAI1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-050-SCAM1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-051-KCHU1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-052-KPLU1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-053-KWAL1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-054-ZCAR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-055-JLAW1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-056-ASUR1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-057-IZIK1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-058-SAPE1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-059-TBOY1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-060-MCHI1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-061-JCLA1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-062-SSOH1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-063-AWIL1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-064-STHO1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-065-MBIC1/mus-bid_spm/liking', '/scratch/mquinci/mci_spm_musbid/pre/sub-066-CKOB1/mus-bid_spm/liking'});

%define empty array
youngAdultdata_updated = cell(1, numel(sublistYA));

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
        youngAdultdata_updated{subs}(:, rois) = summary_data(marsy(Y));
    end
end

%reorder so HG are the first and second columns
for i = 1:size(youngAdultdata_updated,2)
       
       youngAdultdata_updated{i} = (youngAdultdata_updated{i}(:, [18, 19, 1:17, 20:38]));
    
end

%save data to a .mat file
cd '/scratch/mquinci/';
save('youngAdultdata_Wang_updated24.mat','youngAdultdata_updated');


%load onsets
cd '/home/mquinci/ROI_Analyses/Liking/Young_Adults'
load('onsets.txt');

%isolate the 42 trs for each trial that correspond with the music listening
updatedDataYA = cell(1, numel(sublistYA));
for i = 1:size(youngAdultdata_updated,2)
    
   for j =  1:length(onsets)

       updatedDataYA{i}((42*(j-1)+1):42*j, :) = youngAdultdata_updated{i}(onsets(j):onsets(j) + 41, :);
       
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

%orderedDataYA = cell(1, numel(sublistYA));
%for i = 1:size(meanDataAllYA,2)
       
       %orderedDataYA{i} = cat(2,meanDataAllYA{i},stim_orderliking{i});
       %orderedDataYA{i} = sortrows(orderedDataYA{i}(2:end, :), 39);
    
%end



liking_full_order = cell(1,numel(StimListYA));
hate_YA = cell(1,numel(sublistYA));
neutral_YA = cell(1, numel(sublistYA));
like_YA = cell(1, numel(sublistYA));
love_YA = cell(1, numel(sublistYA));
cd '/home/mquinci/ROI_Analyses/Liking/Young_Adults'
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


xvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr'};
yvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr'};


meanRhateYA = nanmean(RhateYA,3);
figure;
hateh = heatmap(xvalues, yvalues, meanRhateYA([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1],'GridVisible','off');
colormap(parula)

meanRneutralYA = nanmean(RneutralYA,3);
figure;
heatmap(xvalues, yvalues, meanRneutralYA([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRlikeYA = nanmean(RlikeYA,3);
figure;
heatmap(xvalues, yvalues, meanRlikeYA([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRloveYA= nanmean(RloveYA,3);
figure;
heatmap(xvalues, yvalues, meanRloveYA([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

%Isolate the auditory-auditory section of the matrix
rdatahateaudYA = RhateYA([1:2, 4:19],[1:2, 4:19],:);
rdataneutralaudYA = RneutralYA([1:2, 4:19],[1:2, 4:19],:);
rdatalikeaudYA = RlikeYA([1:2, 4:19],[1:2, 4:19],:);
rdataloveaudYA = RloveYA([1:2, 4:19],[1:2, 4:19],:);

%Isolate the top half of the corr matrix for each liking rating (aud-aud)
for sub = 1:length(sublistYA) 
    
    temp = rdatahateaudYA(:,:,sub);
    rdatahateaud_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatahateaud_meanYA = mean(rdatahateaud_allYA);
    
    temp = rdataneutralaudYA(:,:,sub);
    rdataneutralaud_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataneutralaud_meanYA = mean(rdataneutralaud_allYA);
    
    temp = rdatalikeaudYA(:,:,sub);
    rdatalikeaud_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatalikeaud_meanYA = mean(rdatalikeaud_allYA);
    
    temp = rdataloveaudYA(:,:,sub);
    rdataloveaud_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataloveaud_meanYA = mean(rdataloveaud_allYA);
    
end


%Isolate the reward-reward section of the matrix
rdatahaterewYA = RhateYA([21:38],[21:38],:);
rdataneutralrewYA = RneutralYA([21:38],[21:38],:);
rdatalikerewYA = RlikeYA([21:38],[21:38],:);
rdataloverewYA = RloveYA([21:38],[21:38],:);

%Isolate the top half of the corr matrix for each liking rating (rew-rew)
for sub = 1:length(sublistYA) 
    
    temp = rdatahaterewYA(:,:,sub);
    rdatahaterew_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatahaterew_meanYA = mean(rdatahaterew_allYA);
    
    temp = rdataneutralrewYA(:,:,sub);
    rdataneutralrew_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataneutralrew_meanYA = mean(rdataneutralrew_allYA);
    
    temp = rdatalikerewYA(:,:,sub);
    rdatalikerew_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatalikerew_meanYA = mean(rdatalikerew_allYA);
    
    temp = rdataloverewYA(:,:,sub);
    rdataloverew_allYA(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataloverew_meanYA = mean(rdataloverew_allYA);
    
end

%Isolate the auditory-reward section of the matrix
rdatahateaudrewYA = RhateYA([1:2, 4:19],[21:38],:);
rdataneutralaudrewYA = RneutralYA([1:2, 4:19],[21:38],:);
rdatalikeaudrewYA = RlikeYA([1:2, 4:19],[21:38],:);
rdataloveaudrewYA = RloveYA([1:2, 4:19],[21:38],:);

%Get the mean for each participant across the matrix (aud-rew)
for sub = 1:length(sublistYA) 
    
    rdatahateaudrew_meanYA(:,sub) = mean(rdatahateaudrewYA(:,:,sub), 'all');
    
    rdataneutralaudrew_meanYA(:,sub) = mean(rdataneutralaudrewYA(:,:,sub), 'all');
    
    rdatalikeaudrew_meanYA(:,sub) = mean(rdatalikeaudrewYA(:,:,sub), 'all');
    
    rdataloveaudrew_meanYA(:,sub) = mean(rdataloveaudrewYA(:,:,sub), 'all');
   
end

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
