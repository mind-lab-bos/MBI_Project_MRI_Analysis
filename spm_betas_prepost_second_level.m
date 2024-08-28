function [outFiles] = spm_betas_prepost_second_level(batch,MCIstatus,intervention,subjectInfo,params)
    %{
        BMK - 16Jul2024
    plots of subject level pre v post betas for contrasts of interest

    %}
    

    % handle each condition in the model
    try
        nvars = size(params.glms,2); 
        fprintf('working on %d models...\n',nvars)
    catch
        fprintf('WARNING no variables/conditions found in paramter file!!!\n')
        return
    end

    % grab only subjects from batch(s) of interest
    % grab only subjects from batch(s) of interest
    batch{end+1} = '';
    MCIstatus{end+1} = '';
    intervention{end+1} = '';
    tmpBatch = cellfun(@(c)strcmp(c,{subjectInfo.batch}),batch,'UniformOutput',false);
    tmpStatus = cellfun(@(c)strcmp(c,{subjectInfo.MCI}),MCIstatus,'UniformOutput',false);
    tmpInter = cellfun(@(c)strcmp(c,{subjectInfo.intervention}),intervention,'UniformOutput',false);
    subidx = intersect(intersect(find(sum(vertcat(tmpBatch{:}))>0),find(sum(vertcat(tmpStatus{:}))>0)) ,find(sum(vertcat(tmpInter{:})>0)));
       if isempty(subidx)
        fprintf('WARNING no subject found for input combination!!!\n')
        return
    end
    isub=1;
    for iidx=1:length(subidx)
        sub2use(isub) = subjectInfo(subidx(iidx));
        isub = isub + 1;
    end

    fprintf('\n\n\n~~~~~Building second level beta plots with %d subjects~~~~~\n\n\n',length({sub2use.subjectID}))
    fprintf(batch{:})
    fprintf(intervention{:})
    
    cbatches = sprintf('%s_',batch{:}); cbatches=cbatches(1:end-1); cbatches=[cbatches sprintf('%s_',MCIstatus{:})]; cbatches=cbatches(1:end-1);
    cbatches=[cbatches sprintf('%s_',intervention{:})]; cbatches=cbatches(1:end-2);    
    % import all subs beta .csvs
    for isub=1:length(sub2use)
        for isess=1:length(sub2use(isub).subjectNUM)
            inpath = fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{isess})); 
            if exist(fullfile(inpath,sprintf('TconBetasOfInterest_%s_%s-sub-%s.csv',sub2use(isub).subjectID,sub2use(isub).sessNames{isess},sub2use(isub).subjectNUM{isess})),'file')~=2 %r
                fprintf(['WANRING: cannot find betas at ' savepath '!!!\n']);
                continue
            end
            if isub==1 && isess==1
                roiBetaTable = readtable(fullfile(inpath,sprintf('TconBetasOfInterest_%s_%s-sub-%s.csv',sub2use(isub).subjectID,sub2use(isub).sessNames{isess},sub2use(isub).subjectNUM{isess})));        
            else
                roiBetaTable = [roiBetaTable;readtable(fullfile(inpath,sprintf('TconBetasOfInterest_%s_%s-sub-%s.csv',sub2use(isub).subjectID,sub2use(isub).sessNames{isess},sub2use(isub).subjectNUM{isess})))];
            end
        end %for isess=1:length(sub2use(isub).subjectNUM)
    end %for isub=1:length(sub2use)
    
    writetable(roiBetaTable,fullfile(params.paths.spm2ndLevel,sprintf('%s_TconBetasOfInterest.csv',cbatches)));
    
    %grab the ROis in the csv
    colnames = roiBetaTable.Properties.VariableNames;
    % Find the intersection
    [roinames,ia] = setdiff(colnames,{'ID','NUM','session','contrast'}, 'stable');

    %Plotting the data
    uConts = unique(roiBetaTable.contrast);
    for icont=1:length(uConts)
        tiledlayout(length(roinames)/3,length(roinames)/2)
        for iroi=1:length(roinames)
            pre_msk = ~cellfun(@isempty,regexp('pre',roiBetaTable.session,'match')) ;
            post_msk = ~cellfun(@isempty,regexp('post',roiBetaTable.session,'match'));
            roiBetaTablePre = roiBetaTable(pre_msk,:);
            roiBetaTablePost = roiBetaTable(post_msk,:);
        
            roiBetaTablePreCconn = roiBetaTablePre(roiBetaTablePre.contrast== string(uConts(icont)),:);
            roiBetaTablePostCconn = roiBetaTablePost(roiBetaTablePost.contrast==string(uConts(icont)),:);
            
            both = [roiBetaTablePreCconn.(sprintf('%s',roinames{iroi})) roiBetaTablePostCconn.(sprintf('%s',roinames{iroi}))];
            nexttile
            title(extractAfter(replace(replace([extractBefore(roinames{iroi},'_nii')],'_',' '),digitsPattern(1),' '),' '))
            plotSpread(both,'xNames',{'pre','post'},'distributionColors',{'red' 'blue'}, ...
                'yLabel', 'Beta' )
            % Loop for line connections
            for k = 1 : size(roiBetaTablePreCconn.(sprintf('%s',roinames{iroi})), 1)
              plot([roiBetaTablePreCconn.(sprintf('%s',roinames{iroi}))(k), roiBetaTablePostCconn.(sprintf('%s',roinames{iroi}))(k)],...
                'ks-', 'MarkerSize', 5);
              hold on
            end
            hold off
        end
        drawnow;
        f = gcf;
        exportgraphics(f,fullfile(params.paths.spm2ndLevel,[cbatches '_preVpost_beta_' uConts{icont} '.png']),'Resolution',300)
        close all
    end

    outFiles = 0;
end
