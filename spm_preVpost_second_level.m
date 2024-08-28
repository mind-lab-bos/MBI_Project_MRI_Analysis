function [outFiles3] = spm_preVpost_second_level(batch,intervention,subjectInfo,params)
    %{
        BMK - 05May2024
    Hardcoded group level Post>Pre contrasts

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
    subidx = intersect(find(strcmp({subjectInfo.batch},batch)==1),find(strcmp({subjectInfo.intervention},intervention)==1));
    if isempty(subidx)
        fprintf('WARNING no subject found for %s %s!!!\n',batch, intervention)
        return
    end
    isub=1;
    for iidx=1:length(subidx)
        sub2use(isub) = subjectInfo(subidx(iidx));
        isub = isub + 1;
    end

    fprintf('\n\n\n~~~~~Building second level pre V post model(s) with %d %s %s subjects~~~~~\n\n\n',length({sub2use.subjectID}),batch,intervention)
    
    % create a table of all subjects included in this contrast
    varTypes = ["string" "string" "string" "string" "string" "string" repelem(["string"],length(params.allContrasts))];
    varNames = ["ID" "NUM" "batch" "status" "intervention" "session" params.allContrasts]; 
    outTable = table('Size',[0 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);

    modelCounter = 1; %keeps track of second level models
    for ivar=1:nvars %loop over each variable/model 
        
        % add the contrasts to the model
        for icont=1:length(params.glms(ivar).Tconts) %loop through eav t contrast for this model 
            % copy factorial_design params in one go
            designstats.matlabbatch{modelCounter}.spm.stats.factorial_design = params.fmri.factorial_design;
            designstats.matlabbatch{modelCounter}.spm.stats.factorial_design.des.pt.gmsca =  params.fmri.factorial_design_paired.des.pt.gmsca;
            designstats.matlabbatch{modelCounter}.spm.stats.factorial_design.des.pt.ancova = params.fmri.factorial_design_paired.des.pt.ancova;
        
            % save dir (one for each subXmodel
            savepath = fullfile(params.paths.spm2ndLevel,'postVpre',params.glms(ivar).name,params.glms(ivar).Tconts{icont},[batch '_' intervention]); 
            mkdir(savepath);
            if exist(fullfile(savepath,'SPM.mat'),'file')==2 %remove existing SPM.mat if present to avoid GUI prompt
                delete(fullfile(savepath,'SPM.mat'));
            end
        
            subPairCounter = 0; % for each subject, figure out the con images to load
            fprintf('\n\nWorking on contrast: %s\n\',params.glms(ivar).Tconts{icont})
            for isub=1:length({sub2use.subjectID}) %loop through each subject and grab SPM file with contrast info
                outTable((2*isub)-1,1) = {sub2use(isub).subjectID}; %update table that will hold all are onset and response info
                try %pre contrast
                    clear cSPM1 curcont postConName cstrng
                    cSPM1 = load(fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{1}),params.glms(ivar).name,'SPM.mat'));
                    curcont = find(ismember({cSPM1.SPM.xCon.name},params.glms(ivar).Tconts{icont})==1); %find our contrast of interest 
                    preConName = cSPM1.SPM.xCon(curcont).Vcon.fname;
                    cstrng = sprintf('['); %gather contrasts for table
                    for itmp=1:length(cSPM1.SPM.xCon(curcont).c)
                        cstrng = [cstrng sprintf('%d ',cSPM1.SPM.xCon(curcont).c(itmp))];
                    end
                    cstrng = [cstrng(1:end-1) ']'];
                    outTable((2*isub)-1,find(ismember(outTable.Properties.VariableNames,params.glms(ivar).Tconts{icont})==1)) = {cstrng}; %update table
                    outTable((2*isub)-1,2:6) = {sub2use(isub).subjectNUM{1},sub2use(isub).batch,sub2use(isub).MCI,sub2use(isub).intervention,sub2use(isub).sessNames{1}};
                   
                catch
                    fprintf('WARNING could not find contrasts for subject: %s; session: %s; sub-%s;!!!\n',sub2use(isub).subjectID,'pre','NA')
                    outTable((2*isub)-1,find(ismember(outTable.Properties.VariableNames,params.glms(ivar).Tconts{icont})==1)) = {'[]'}; %update table
                    outTable((2*isub)-1,2:6) = {'NA',sub2use(isub).batch,sub2use(isub).MCI,sub2use(isub).intervention,'pre'};
                    continue
                end
                outTable((2*isub),1) = {sub2use(isub).subjectID}; %update table that will hold all are onset and response info
                try %post contrast
                    clear cSPM2 curcont postConName cstrng
                    cSPM2 = load(fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{2}),params.glms(ivar).name,'SPM.mat'));
                    curcont = find(ismember({cSPM2.SPM.xCon.name},params.glms(ivar).Tconts{icont})==1); %find our contrast of interest 
                    postConName = cSPM2.SPM.xCon(curcont).Vcon.fname;
                    cstrng = sprintf('['); %gather contrasts for table
                    for itmp=1:length(cSPM2.SPM.xCon(curcont).c)
                        cstrng = [cstrng sprintf('%d ',cSPM2.SPM.xCon(curcont).c(itmp))];
                    end
                    cstrng = [cstrng(1:end-1) ']'];
                    outTable((2*isub),find(ismember(outTable.Properties.VariableNames,params.glms(ivar).Tconts{icont})==1)) = {cstrng}; %update table
                    outTable((2*isub),2:6) = {sub2use(isub).subjectNUM{2},sub2use(isub).batch,sub2use(isub).MCI,sub2use(isub).intervention,sub2use(isub).sessNames{2}};
                   
                catch
                    fprintf('WARNING could not find contrasts for subject: %s; session: %s; sub-%s;!!!\n',sub2use(isub).subjectID,'post','NA')
                    outTable((2*isub),find(ismember(outTable.Properties.VariableNames,params.glms(ivar).Tconts{icont})==1)) = {'[]'}; %update table
                    outTable((2*isub),2:6) = {'NA',sub2use(isub).batch,sub2use(isub).MCI,sub2use(isub).intervention,'pre'};
                    continue
                end
                subPairCounter = subPairCounter + 1;
                designstats.matlabbatch{modelCounter}.spm.stats.factorial_design.des.pt.pair(subPairCounter).scans = {
                     [fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{1}),params.glms(ivar).name,preConName) ',1']
                     [fullfile(params.paths.spm1stLevel,sprintf('sub-%s',sub2use(isub).subjectNUM{2}),params.glms(ivar).name,postConName) ',1']};
                
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
    writetable(outTable,fullfile(params.paths.spm2ndLevel,'postVpre',sprintf('contrasts_postVpre.csv')));
   
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
        
        contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.name = sprintf('F_Post_V_Pre__%s',capthtmp(end:-1:1));
        contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
        contrast.matlabbatch{imod}.spm.stats.con.consess{1}.fcon.weights = [-1 1];
        contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.name = sprintf('T_Post_V_Pre__%s',capthtmp(end:-1:1));
        contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        contrast.matlabbatch{imod}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
        
    end %imod=1:length(outFiles2)

    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles3 = spm_jobman('run',contrast.matlabbatch);

end
