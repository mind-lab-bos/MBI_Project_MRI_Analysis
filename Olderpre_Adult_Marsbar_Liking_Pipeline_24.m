%Older Adults pre

%defining and loading rois
cd '/Volumes/Promise_Pegasus/mci_spm_musbid/ROI_pipeline_OA/Wang_ROIs'
roi_list = dir('*.mat');

%define subjects
cd '/Volumes/Promise_Pegasus/mci_spm_musbid/pre';
sublistOA = struct('path',{'/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-01-FKEE1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-03-JSCH1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-05-LGEN1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-06-GHER1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-08-LDIR1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-010-JPRE1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-011-CWEI1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-012-MBLO1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-015-EROS1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-016-TCIE1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-018-GCAS1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-019-DFAI1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-021-JKAY1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-022-DCAT1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-023-LOLE1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-028-SCAE1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-029-EHER1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-030-LJAC1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-032-ASHA1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-036-GSHA1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-037-AEHR1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-039-ELAR1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-040-BFOR1/mus-bid_spm/liking', '/Volumes/Promise_Pegasus/mci_spm_musbid/pre/sub-047-TSTA1/mus-bid_spm/liking'});

%define empty array
olderAdultdata_updated = cell(1, numel(sublistOA));

marsbar('on')
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
        olderAdultdata_updated{subs}(:, rois) = summary_data(marsy(Y));
    end
end

%reorder so HG are the first and second columns
for i = 1:size(olderAdultdata_updated,2)
       
       olderAdultdata_updated{i} = (olderAdultdata_updated{i}(:, [18, 19, 1:17, 20:38]));
    
end

%save data to a .mat file
cd '/Volumes/Promise_Pegasus/mci_spm_musbid/ROI_pipeline_OA/';
save('olderAdultdata_Wang_updated.mat','olderAdultdata_updated');


%load onsets
cd '/scratch/mquinci/ROI_Analyses/Liking/Young_Adults'
load('onsets.txt');

%isolate the 42 trs for each trial that correspond with the music listening
updatedDataOA = cell(1, numel(sublistOA));
for i = 1:size(olderAdultdata_updated,2)
    
   for j =  1:length(onsets)

       updatedDataOA{i}((42*(j-1)+1):42*j, :) = olderAdultdata_updated{i}(onsets(j):onsets(j) + 41, :);
       
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
cd '/scratch/mquinci/ROI_Analyses/ROI_pipeline_OA/stim_order'
StimList = dir('*.txt');
stim_order  = cell(1, numel(StimList));
for stim = 1:numel(StimList)
  FileData     = load(StimList(stim).name);
  stim_order{stim} = FileData;
  stim_order{stim}(1,:)= [];
  stim_orderliking{stim} = stim_order{stim}(:,2);
end
cd ../


%load stim order
%orderedData = cell(1, numel(sublistOA));
%for i = 1:size(meanDataAllOA,2)
       
       %orderedData{i} = cat(2,meanDataAllOA{i},stim_orderliking{i});
       %orderedData{i} = sortrows(orderedData{i}(2:end, :), 39);
    
%end



liking_full_order = cell(1,numel(StimList));
hate_OA = cell(1,numel(sublistOA));
neutral_OA = cell(1, numel(sublistOA));
like_OA = cell(1, numel(sublistOA));
love_OA = cell(1, numel(sublistOA));
cd '/scratch/mquinci/ROI_Analyses/Liking/Young_Adults/'
load('stimstarts.txt');
for sub = 1:size(updatedDataOA, 2)
    
    hate = 1;
    neutral = 1;
    like = 1;
    love = 1;
    
    for stim = 1:24
        
        if stim_orderliking{sub}(stim) == 1
            hate_OA{sub}((hate-1)*42+1:hate*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            hate=hate+1;
            
        elseif stim_orderliking{sub}(stim) == 2
            neutral_OA{sub}((neutral-1)*42+1:neutral*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            neutral=neutral+1;
            
        elseif stim_orderliking{sub}(stim) == 3
            like_OA{sub}((like-1)*42+1:like*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            like=like+1;
            
        elseif stim_orderliking{sub}(stim) == 4
            love_OA{sub}((love-1)*42+1:love*42,:) = updatedDataOA{sub}((stimstarts(stim):stimstarts(stim+1)-1), :);
            love=love+1;
            
        else 
            
        end
        
    end
end


for sub = 1:length(sublistOA)
 if numel(hate_OA{sub}) > 0   
    Rhate(:,:,sub) = corr(hate_OA{sub});
   else
     Rhate(:,:,sub) = NaN;
 end
   
 if numel(neutral_OA{sub}) > 0   
    Rneutral(:,:,sub) = corr(neutral_OA{sub});
   else
     Rneutral(:,:,sub) = NaN;
 end
 
 if numel(like_OA{sub}) > 0   
    Rlike(:,:,sub) = corr(like_OA{sub});
   else
     Rlike(:,:,sub) = NaN;
   end
   
 if numel(love_OA{sub}) > 0   
    Rlove(:,:,sub) = corr(love_OA{sub});
   else
     Rlove(:,:,sub) = NaN;
   end
 
 
end


xvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr'};
yvalues = {'HGl', 'HGr','aSTGl','aSTGr','pSTGl', 'pSTGr','aMTGl','aMTGr','pMTGl', 'pMTGr', 'toMTGl', 'toMTGr', 'aITGl', 'aITGr', 'pITGl', 'pITGr', 'toITGl', 'toITGr', 'ICI', 'ICr', 'AC', 'PC', 'FOrbl', 'FOrbr', 'Caudatel', 'Caudater', 'Putamenl', 'Putamenr', 'Palliduml', 'Pallidumr', 'Hippocampusl', 'Hippocampusr', 'Amygdalal', 'Amygdalar', 'Accumbensl', 'Accumbensr'};


meanRhate = nanmean(Rhate,3);
figure;
hateh = heatmap(xvalues, yvalues, meanRhate([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRneutral = nanmean(Rneutral,3);
figure;
heatmap(xvalues, yvalues, meanRneutral([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRlike = nanmean(Rlike,3);
figure;
heatmap(xvalues, yvalues, meanRlike([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)

meanRlove = nanmean(Rlove,3);
figure;
heatmap(xvalues, yvalues, meanRlove([1:2, 4:19, 21:38],[1:2, 4:19, 21:38]), 'ColorLimits',[0 1], 'GridVisible','off');
colormap(parula)
colorbar('southoutside')

%Isolate the auditory-auditory section of the matrix
rdatahateaud = Rhate([1:2, 4:19],[1:2, 4:19],:);
rdataneutralaud = Rneutral([1:2, 4:19],[1:2, 4:19],:);
rdatalikeaud = Rlike([1:2, 4:19],[1:2, 4:19],:);
rdataloveaud = Rlove([1:2, 4:19],[1:2, 4:19],:);

%Isolate the top half of the corr matrix for each liking rating (aud-aud)
for sub = 1:length(sublistOA) 
    
    temp = rdatahateaud(:,:,sub);
    rdatahateaud_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatahateaud_mean = mean(rdatahateaud_all);
    
    temp = rdataneutralaud(:,:,sub);
    rdataneutralaud_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataneutralaud_mean = mean(rdataneutralaud_all);
    
    temp = rdatalikeaud(:,:,sub);
    rdatalikeaud_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatalikeaud_mean = mean(rdatalikeaud_all);
    
    temp = rdataloveaud(:,:,sub);
    rdataloveaud_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataloveaud_mean = mean(rdataloveaud_all);
    
end


%Isolate the reward-reward section of the matrix
rdatahaterew = Rhate([21:38],[21:38],:);
rdataneutralrew = Rneutral([21:38],[21:38],:);
rdatalikerew = Rlike([21:38],[21:38],:);
rdataloverew = Rlove([21:38],[21:38],:);

%Isolate the top half of the corr matrix for each liking rating (rew-rew)
for sub = 1:length(sublistOA) 
    
    temp = rdatahaterew(:,:,sub);
    rdatahaterew_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatahaterew_mean = mean(rdatahaterew_all);
    
    temp = rdataneutralrew(:,:,sub);
    rdataneutralrew_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataneutralrew_mean = mean(rdataneutralrew_all);
    
    temp = rdatalikerew(:,:,sub);
    rdatalikerew_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdatalikerew_mean = mean(rdatalikerew_all);
    
    temp = rdataloverew(:,:,sub);
    rdataloverew_all(:,sub) = temp(find(~tril(ones(size(temp)))));
    rdataloverew_mean = mean(rdataloverew_all);
    
end

%Isolate the auditory-reward section of the matrix
rdatahateaudrew = Rhate([1:2, 4:19],[21:38],:);
rdataneutralaudrew = Rneutral([1:2, 4:19],[21:38],:);
rdatalikeaudrew = Rlike([1:2, 4:19],[21:38],:);
rdataloveaudrew = Rlove([1:2, 4:19],[21:38],:);

%Get the mean for each participant across the matrix (aud-rew)
for sub = 1:length(sublistOA) 
    
    rdatahateaudrew_mean(:,sub) = mean(rdatahateaudrew(:,:,sub), 'all');
    
    rdataneutralaudrew_mean(:,sub) = mean(rdataneutralaudrew(:,:,sub), 'all');
    
    rdatalikeaudrew_mean(:,sub) = mean(rdatalikeaudrew(:,:,sub), 'all');
    
    rdataloveaudrew_mean(:,sub) = mean(rdataloveaudrew(:,:,sub), 'all');
   
end

%combine the three files together and get the mean (one mean for hate
%aud-aud, one mean for hate rew-rew, and one mean for hate aud-rew. These
%would all go in separate plots. Dont forget about error bars

%Auditory-Auditory
for sub = 1:length(sublistOA) 
auditory_auditory(sub,:) = [rdatahateaud_mean(:,sub), rdataneutralaud_mean(:,sub), rdatalikeaud_mean(:,sub), rdataloveaud_mean(:,sub)];
end

SE_aud = nanstd(auditory_auditory)/sqrt(length(auditory_auditory));

auditory_auditorymean = nanmean(auditory_auditory);
X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
figure;  
bar(X,auditory_auditorymean, 'FaceColor', [0.6350 0.0780 0.1840])
title('Older Adults Connectivity between Auditory-Auditory Regions')

hold on

er = errorbar(X,auditory_auditorymean,SE_aud);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%Reward-Reward
for sub = 1:length(sublistOA) 
reward_reward(sub,:) = [rdatahaterew_mean(:,sub), rdataneutralrew_mean(:,sub), rdatalikerew_mean(:,sub), rdataloverew_mean(:,sub)];
end

SE_rew = nanstd(reward_reward)/sqrt(length(reward_reward));

reward_rewardmean = nanmean(reward_reward);
X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
figure;  
bar(X,reward_rewardmean, 'FaceColor', [0.6350 0.0780 0.1840])
title('Older Adults Connectivity between Reward-Reward Regions')

hold on

er = errorbar(X,reward_rewardmean,SE_rew);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

%Auditory-Reward
for sub = 1:length(sublistOA) 
auditory_reward(sub,:) = [rdatahateaudrew_mean(:,sub), rdataneutralaudrew_mean(:,sub), rdatalikeaudrew_mean(:,sub), rdataloveaudrew_mean(:,sub)];
end

SE_audrew = nanstd(auditory_reward)/sqrt(length(auditory_reward));

auditory_rewardmean = nanmean(auditory_reward);
X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
figure;  
bar(X,auditory_rewardmean, 'FaceColor', [0.6350 0.0780 0.1840])
title('Older Adults Connectivity between Auditory-Reward Regions')

hold on

er = errorbar(X,auditory_rewardmean,SE_audrew);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

