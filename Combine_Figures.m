
%https://www.mathworks.com/matlabcentral/answers/102220-how-do-i-place-errorbars-on-my-grouped-bar-graph-using-function-errorbar-in-matlab#accepted_answer_111566

X = categorical({'Hate','Neutral','Like','Love'});
X = reordercats(X,{'Hate','Neutral','Like','Love'});
%T = table(X', auditory_rewardmeanYA', auditory_rewardmean');
auditory_rewardmeanAll = [auditory_rewardmeanYA; auditory_rewardmean]';

SE_audrewYA = nanstd(auditory_rewardYA)/sqrt(length(auditory_rewardYA));
SE_audrewOA = nanstd(auditory_reward)/sqrt(length(auditory_reward));
SE_audrewAll = [SE_audrewYA; SE_audrewOA]'

figure;  
h = bar(X,auditory_rewardmeanAll, 'grouped')
ylim([.3 .6])
title(sprintf('Connectivity between Auditory-Reward Regions \n in Young and Older Adults'))

% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'Young Adults','Older Adults'}')
% Legend will show names for each color
legend() 

%hold on

%er = errorbar(X, auditory_rewardmeanAll, SE_audrewAll);    
%er.Color = [0 0 0];                            
%er.LineStyle = 'none';  

%hold off

%Auditory Connectivity
auditory_auditorymeanAll = [auditory_auditorymeanYA; auditory_auditorymean]';

figure;  
h = bar(X,auditory_auditorymeanAll, 'grouped')
ylim([.3 .6])
title(sprintf('Connectivity between Auditory Regions \n in Young and Older Adults'))

% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'Young Adults','Older Adults'}')
% Legend will show names for each color
legend() 

%Reward Connectivity
reward_rewardmeanAll = [reward_rewardmeanYA; reward_rewardmean]';

figure;  
h = bar(X,reward_rewardmeanAll, 'grouped')
ylim([.3 .6])
title(sprintf('Connectivity between Reward Regions \n in Young and Older Adults'))

% set 3 display names for the 3 handles
set(h, {'DisplayName'}, {'Young Adults','Older Adults'}')
% Legend will show names for each color
legend()




y=[auditory_rewardmeanYA; auditory_rewardmean]'
StdErr = [SE_audrewYA; SE_audrewOA]'                                               % Create Standard Error Matrix
figure
hBar = bar(y, 0.8);                                                     % Return bar Handle
for k1 = 1:size(y,2)
    ctr(k1,:) = bsxfun(@plus, hBar(k1).XData, hBar(k1).XOffset');       % Note: XOffset Is An Undocumented Feature, This Selects The bar Centres
    ydt(k1,:) = hBar(k1).YData;                                         % Individual Bar Heights
end
hold on
errorbar(ctr, ydt, StdErr, '.r')                                        % Plot Error Bars
set(gca, 'XTickLabel', {'Rural' 'Exurban' 'Suburban' 'Urban'});





ytest=[0.095	0.037	0.132;	0.1	0.058	0.158;	0 0 0.02;	0 0 0]
StdErrtest = rand(3,4)*0.01   

