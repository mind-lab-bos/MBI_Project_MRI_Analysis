function [outFiles3] = spm_preXpost_second_level(batch,MCIstatus,intervention,sessnames,subjectInfo,params)
    %{
        BMK - 05May2024
    Hardcoded group level Post and Pre contrasts

    %}
    
    % set up batch job
    spm('defaults','fmri');
    spm_jobman('initcfg');
    designstats = struct;

    % handle each condition in the model
    try
        nvars = size(params.glms,2); 
        fprintf('working on %d models...\n',nvars)
    catch
        fprintf('WARNING no variables/conditions found in paramter file!!!\n')
        return
    end
    
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

    fprintf('\n\n\n~~~~~Building second level pre X post model(s) with %d subjects~~~~~\n\n\n',length({sub2use.subjectID}))
    fprintf(batch{:})
    fprintf(intervention{:})
    
    cbatches = sprintf('%s_',batch{:}); cbatches=cbatches(1:end-1); cbatches=[cbatches sprintf('%s_',MCIstatus{:})]; cbatches=cbatches(1:end-1);
    cbatches=[cbatches sprintf('%s_',intervention{:})]; cbatches=cbatches(1:end-2);

    modelCounter = 1; %keeps track of second level models
    %sessnames = unique(vertcat(subjectInfo.sessNames),'stable');
    for isess=1:length(sessnames) %loop over each session
        % create a table of all subjects included in this contrast
        varTypes = ["string" "string" "string" "string" "string" "string" repelem(["string"],length(params.allContrasts))];
        varNames = ["ID" "NUM" "batch" "status" "intervention" "session" params.allContrasts]; 
        outTable = table('Size',[0 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
        for ivar=1:nvars %loop over each variable/model 
    
            % add the contrasts to the model
            for icont=1:length(params.glms(ivar).Tconts) %loop through eav t contrast for this model 
                % copy factorial_design params in one go
                designstats.matlabbatch{modelCounter}.spm.stats.factorial_design = params.fmri.factorial_design;
            
                % save dir (batch-interventionXmodel)
                savepath = fullfile(params.paths.spm2ndLevel,sessnames{isess},params.glms(ivar).name,params.glms(ivar).Tconts{icont},cbatches); 
                mkdir(savepath);
                if exist(fullfile(savepath,'SPM.mat'),'file')==2 %remove existing SPM.mat if present to avoid GUI prompt
                    delete(fullfile(savepath,'SPM.mat'));
                end
            
                subPairCounter = 0; % for each subject, figure out the con images to load
                fprintf('\n\nWorking on contrast: %s\n\',params.glms(ivar).Tconts{icont})
                for isub=1:length({sub2use.subjectID}) %loop through each subject and grab SPM file with contrast info
                   
                    try %current contrast
                        outTable(isub,1:6) = {sub2use(isub).subjectID,sub2use(isub).subjectNUM{isess},sub2use(isub).batch,sub2use(isub).MCI,sub2use(isub).intervention,sub2use(isub).sessNames{isess}}; %update table that will hold all are onset and response info
                        clear cSPM1 curcont preConName
                        cSPM1 = load(fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{isess}),params.glms(ivar).name,'SPM.mat'));
                        curcont = find(ismember({cSPM1.SPM.xCon.name},params.glms(ivar).Tconts{icont})==1); %find our contrast of interest 
                        preConName = cSPM1.SPM.xCon(curcont).Vcon.fname;
                        cstrng = sprintf('['); %gather contrasts for table 
                        for itmp=1:length(cSPM1.SPM.xCon(curcont).c) %   THERE CHANCE THAT THIS INDEX IS MESSED UP WHEN MISSING A CONT...
                            cstrng = [cstrng sprintf('%d ',cSPM1.SPM.xCon(curcont).c(itmp))];
                        end
                        cstrng = [cstrng(1:end-1) ']'];
                        outTable(isub,find(ismember(outTable.Properties.VariableNames,params.glms(ivar).Tconts{icont})==1)) = {cstrng}; %update table
                    catch
                        fprintf('WARNING could not find contrasts for subject: %s; session: %s; sub-%s;!!!\n',sub2use(isub).subjectID,sessnames{isess},'NA')
                        outTable(isub,1:6) = {sub2use(isub).subjectID,'NA',sub2use(isub).batch,sub2use(isub).MCI,sub2use(isub).intervention,sessnames{isess}}; %update table that will hold all are onset and response info
                        outTable(isub,find(ismember(outTable.Properties.VariableNames,params.glms(ivar).Tconts{icont})==1)) = {'[]'}; %update table
                        continue
                    end
                    subPairCounter = subPairCounter + 1;
                    designstats.matlabbatch{modelCounter}.spm.stats.factorial_design.des.t1.scans{subPairCounter,1} = [fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{isess}),params.glms(ivar).name,preConName) ',1'];
                end %for isub=1:length({sub2use.subjectID})
                if subPairCounter > 2 %make sure we have enough subjects
                    designstats.matlabbatch{modelCounter}.spm.stats.factorial_design.dir = {savepath};
                    modelCounter = modelCounter + 1;
                else
                    fprintf('WARNING not enough subjects (%d) for contrast!!!\n',subPairCounter)
                    continue
                end
            end %for icont=1:length(params.glms(ivar).Tconts)
        end %for ivar=1:nvars
        %save out table with contrasts
        writetable(outTable,fullfile(params.paths.spm2ndLevel,sessnames{isess},sprintf('contrasts_%s.csv',sessnames{isess})));
    end %for isess=1:length(sub2use(isub).sessNames)

    
   
    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles = spm_jobman('run',designstats.matlabbatch);
    
    modelestimation = struct;
    % estimate each of the models
    for imod=1:length(outFiles)
        % copy params and model SPM.mat
        modelestimation.matlabbatch{imod}.spm.stats.fmri_est = params.fmri.fmri_est;
        modelestimation.matlabbatch{imod}.spm.stats.fmri_est.spmmat = outFiles{imod}.spmmat;
    end %imod=1:length(outFiles)

    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles2 = spm_jobman('run',modelestimation.matlabbatch);

    contrast = struct;
    % contrasts for each 2nd level model 
    for imod=1:length(outFiles2)
        [cpath,~,~] = fileparts(outFiles2{imod}.spmmat); %grab the contrast name
        capthtmp = extractBefore(cpath(end-length([char(batch) char(intervention)])-2:-1:1),"/");
        
        % copy params, model SPM.mat, and pre vs post contrast
        contrast.matlabbatch{imod}.spm.stats.con.delete = params.fmri.stats.con.delete;
        contrast.matlabbatch{imod}.spm.stats.con.spmmat = outFiles2{imod}.spmmat;
        contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.weights = 1;
        contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.weights = 1;
        if contains(cpath,'/pre/')
            contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.name = sprintf('F_Pre__%s',capthtmp(end:-1:1));
            contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.name = sprintf('T_Pre__%s',capthtmp(end:-1:1));
        elseif contains(cpath,'/post/')
            contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.name = sprintf('F_Post__%s',capthtmp(end:-1:1));
            contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.name = sprintf('T_Post__%s',capthtmp(end:-1:1));
        end
    end %imod=1:length(outFiles2)

    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles3 = spm_jobman('run',contrast.matlabbatch);

end
