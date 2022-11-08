%Auditory Connectivity
%define the data and error bars for both young and older adults
auditory_meanall=[auditory_auditorymeanYA; auditory_auditorymean]';
SE_audYA = nanstd(auditory_auditoryYA)/sqrt(length(auditory_auditoryYA));
SE_audOA = nanstd(auditory_auditory)/sqrt(length(auditory_auditory));
SE_audAll = [SE_audYA; SE_audOA]';

%plot the figure
figure
hBar = bar(auditory_meanall, 1); 
hBar(1).FaceColor = '#D95319' ;
hBar(2).FaceColor = '#0072BD' ;

%add y limit and title
ylim([.3 .6])
title(sprintf('Connectivity between Auditory Regions \n in Young and Older Adults'))

%add error bars
testxdataaud = hBar.XData
testxoffsetaud = hBar.XOffset
testadditionaud = testxdataaud + testxoffsetaud
testsubtractionaud = testxdataaud - testxoffsetaud
errorlocationsaud = [testadditionaud; testsubtractionaud]
hold on
errorbar(errorlocationsaud', auditory_meanall, SE_audAll, '.k')

%add x ticks 
set(gca, 'XTickLabel', {'Hate' 'Neutral' 'Like' 'Love'});

% legend display
set(hBar, {'DisplayName'}, {'Young Adults','Older Adults'}')
legend({'Young Adults','Older Adults'}) 

%Reward Connectivity
%define the data and error bars for both young and older adults
reward_meanall=[reward_rewardmeanYA; reward_rewardmean]';
SE_rewYA = nanstd(reward_rewardYA)/sqrt(length(reward_rewardYA));
SE_rewOA = nanstd(reward_reward)/sqrt(length(reward_reward));
SE_rewAll = [SE_rewYA; SE_rewOA]';

%plot the figure
figure
hBar = bar(reward_meanall, 1); 
hBar(1).FaceColor = '#D95319' ;
hBar(2).FaceColor = '#0072BD' ;

%add y limit and title
ylim([.3 .6])
title(sprintf('Connectivity between Reward Regions \n in Young and Older Adults'))

%add error bars
testxdatarew = hBar.XData
testxoffsetrew = hBar.XOffset
testadditionrew = testxdatarew + testxoffsetrew
testsubtractionrew = testxdatarew - testxoffsetrew
errorlocationsrew = [testadditionrew; testsubtractionrew]
hold on
errorbar(errorlocationsrew', reward_meanall, SE_rewAll, '.k')

%add x ticks 
set(gca, 'XTickLabel', {'Hate' 'Neutral' 'Like' 'Love'});

% legend display
set(hBar, {'DisplayName'}, {'Young Adults','Older Adults'}')
legend({'Young Adults','Older Adults'}) 


%Auditory Reward Connectivity
%define the data and error bars for both young and older adults
auditory_rewardmeanall=[auditory_rewardmeanYA; auditory_rewardmean]';
SE_audrewYA = nanstd(auditory_rewardYA)/sqrt(length(auditory_rewardYA));
SE_audrewOA = nanstd(auditory_reward)/sqrt(length(auditory_reward));
SE_audrewAll = [SE_audrewYA; SE_audrewOA]';

%plot the figure
figure
hBar = bar(auditory_rewardmeanall, 1);
hBar(1).FaceColor = '#D95319' ;
hBar(2).FaceColor = '#0072BD' ;

%add y limit and title
ylim([.3 .6])
title(sprintf('Connectivity between Auditory and Reward Regions \n in Young and Older Adults'))

%add error bars
testxdata = hBar.XData
testxoffset = hBar.XOffset
testaddition = testxdata + testxoffset
testsubtraction = testxdata - testxoffset
errorlocations = [testaddition; testsubtraction]
hold on
errorbar(errorlocations', auditory_rewardmeanall, SE_audrewAll, '.k')

%add x ticks 
set(gca, 'XTickLabel', {'Hate' 'Neutral' 'Like' 'Love'});

% legend display
set(hBar, {'DisplayName'}, {'Young Adults','Older Adults'}')
legend({'Young Adults','Older Adults'}) 

