function [outFiles] = spm_build_first_level(isub,isess,params)
    %{
        BMK - 01May2024

    %}
    
    fprintf('\n\n\n~~~~~Building model(s) for subject: %s; session: %d; sub-%s;~~~~~\n\n\n',isub.subjectID,isess,isub.subjectNUM{isess})
    
    % could loop through pre and post as two sessions in one model, but for
    % now keep seperate (always use index 1 here) 
    csess = 1; %

    % set up batch job
    spm('defaults','fmri');
    spm_jobman('initcfg');
    designstats = struct;

    % open log file and gather condition, stimulus, and response info
    try
        fid = fopen(fullfile(params.paths.logs,sprintf('%s_%s',isub.batch,isub.sessNames{isess}),isub.logfiles{isess}));
    catch
        fprintf('WANRING no log file found!!!\n')
        return
    end
    %init table to hold stim and response info to use later on
    varTypes = ["double","double","double","double","double","double","double","double","double","double","string"]; %repelem(["double"],32)];
    varNames = ["order","scanOnset","timeOnset","scanQALIKOnset","timeQALIKOnset","scanQAFAMOnset","timeQAFAMOnset","songName","Liking","Familiarity","SongType"];%params.chanNames];
    outTable = table('Size',[24 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
    dumbCounter = 1;
    tline = fgetl(fid);%ignore 000?
    while ischar(tline)
        %disp(tline)
        tline=fgetl(fid);
        if ~ischar(tline)
            break
        end
        tlinesplit = strsplit(tline,', '); %split to get songid
        tlinesplit2 = strsplit(tlinesplit{2},' ');%split to get ratings
        rtlinesplit3 = regexprep(tlinesplit2{2},'[;]','');
        %figure out songType based on song name (first col)
        tmp = find(params.stims.song2cond.SS==str2num(tlinesplit{1}));
        songType = '';
        if ~isempty(tmp)
            %SS song
            songType = 'ss_listening';
        end 
        if isempty(songType)
            tmp2 = find(params.stims.song2cond.FW==str2num(tlinesplit{1}));
            if ~isempty(tmp2)
                %FW song
                songType = 'fw_listening';
            end
        end
        if isempty(songType)
            tmp3 = find(params.stims.song2cond.BP==str2num(tlinesplit{1}));
            if ~isempty(tmp3)
                %BP song
                songType = 'bp_listening';
            end
        end
        %figure out scanOnsets figure out timeOnsets
        if dumbCounter==1
            %first stim first TR onset
            scanOnset = 1;
            timeOnset = 0;
            scanQALIKOnset = params.fmri.musdur + 1;
            timeQALIKOnset = params.fmri.musdur*params.fmri.fmri_spec.timing.RT;
            scanQAFAMOnset = params.fmri.musdur+params.fmri.LIKratingsdur + 1;
            timeQAFAMOnset = params.fmri.musdur*params.fmri.fmri_spec.timing.RT;
        else
            scanOnset = ((dumbCounter-1)*(params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur))+1;%+1;
            timeOnset = (((dumbCounter-1)*(params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur))*params.fmri.fmri_spec.timing.RT);%+params.fmri.fmri_spec.timing.RT;
            scanQALIKOnset = ((dumbCounter-1)*(params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur))+params.fmri.musdur + 1;
            timeQALIKOnset = (((dumbCounter-1)*(params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur))+params.fmri.musdur)*params.fmri.fmri_spec.timing.RT;
            scanQAFAMOnset = ((dumbCounter-1)*(params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur))+params.fmri.musdur+params.fmri.LIKratingsdur+1;
            timeQAFAMOnset = (((dumbCounter-1)*(params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur))+params.fmri.musdur+params.fmri.LIKratingsdur)*params.fmri.fmri_spec.timing.RT;
      
        end
        %update table that will hold all are onset and response info
        outTable(dumbCounter,:) = {dumbCounter,scanOnset,timeOnset,scanQALIKOnset,timeQALIKOnset,scanQAFAMOnset,timeQAFAMOnset,str2num(tlinesplit{1}),str2num(tlinesplit2{1}),str2num(rtlinesplit3),songType}; 
        dumbCounter = dumbCounter + 1;
    end %while ischar(tline)
    fclose(fid);

    % grab the 4d nifti
    try
        niiFile = fullfile(params.paths.data,sprintf('sub-%s',isub.subjectNUM{isess}),isub.niifiles{isess});
        func4D_spm = spm_vol(niiFile); %load in the 4D .nii
        niiSize = size(func4D_spm);
        %func4D_spm_split = spm_file_split(func4D_spm); %get the names of each TR .nii
    catch
        fprintf('WARNING no 4D .nii found!!!\n')
        return
    end

    % handle each condition in the model
    try
        nvars = size(params.glms,2); 
        fprintf('working on %d models...\n',nvars)
    catch
        fprintf('WARNING no variables/conditions found in paramter file!!!\n')
        return
    end
    for ivar=1:nvars %loop over each variable/model
        % copy fmri_spec params in one go
        designstats.matlabbatch{ivar}.spm.stats.fmri_spec = params.fmri.fmri_spec;
    
        % save dir (one for each subXmodel
        savepath = fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),params.glms(ivar).name); 
        mkdir(savepath);
        if exist(fullfile(savepath,'SPM.mat'),'file')==2 %remove existing SPM.mat if present to avoid GUI prompt
            delete(fullfile(savepath,'SPM.mat'));
        end
        designstats.matlabbatch{ivar}.spm.stats.fmri_spec.dir = {savepath};
    
        % add the scans to the model
        designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).scans = {};
        for iimg=1:niiSize(1)
            designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).scans{iimg,1} = [func4D_spm(iimg).fname ',' num2str(iimg) ]; %[func4D_spm_split(iimg).fname ',1'];
        end
        
        % add conditions for each variable/model
        conds2include = {}; %if we are missing a level (icond) need to remove that from the model.
        conds2resp2include = {};
        condsFound = {};
        for icond=1:length(params.glms(ivar).conds) %loop over each condition in variable/model
            switch params.glms(ivar).name
                case 'Liking'
                    matches = regexp(condsFound,['^' extractBefore(params.glms(ivar).conds{icond},'_')],'match');
                    if sum(~cellfun(@isempty,matches))>0 && sum(strcmp(params.glms(ivar).conds{icond},{'hate_liking_rating','neutral_liking_rating','like_liking_rating','love_liking_rating','hate_fam_rating','neutral_fam_rating','like_fam_rating','love_fam_rating'}))>0
                        conds2include{end+1} = params.glms(ivar).conds{icond};
                    elseif sum(strcmp(params.glms(ivar).conds{icond},{'hate_liking_rating','neutral_liking_rating','like_liking_rating','love_liking_rating','hate_fam_rating','neutral_fam_rating','like_fam_rating','love_fam_rating'}))==0
                        if length(outTable.scanOnset(outTable.Liking==params.glms(ivar).conds2resp{icond}))>0  
                            conds2include{end+1} = params.glms(ivar).conds{icond};
                            condsFound{end+1} = extractBefore(params.glms(ivar).conds{icond},'_');
                            conds2resp2include{end+1} = params.glms(ivar).conds2resp{icond};
                        end
                    end
                case 'Familiarity'
                    matches = regexp(condsFound,['^' extractBefore(params.glms(ivar).conds{icond},'_')],'match');
                    if sum(~cellfun(@isempty,matches))>0 && sum(strcmp(params.glms(ivar).conds{icond},{'vunf_liking_rating','unf_liking_rating','fam_liking_rating','vfam_liking_rating','vunf_fam_rating','unf_fam_rating','fam_fam_rating','vfam_fam_rating'}))>0
                        conds2include{end+1} = params.glms(ivar).conds{icond};
                    elseif sum(strcmp(params.glms(ivar).conds{icond},{'vunf_liking_rating','unf_liking_rating','fam_liking_rating','vfam_liking_rating','vunf_fam_rating','unf_fam_rating','fam_fam_rating','vfam_fam_rating'}))==0
                        if length(outTable.scanOnset(outTable.Familiarity==params.glms(ivar).conds2resp{icond}))>0  
                            conds2include{end+1} = params.glms(ivar).conds{icond};
                            condsFound{end+1} = extractBefore(params.glms(ivar).conds{icond},'_');
                            conds2resp2include{end+1} = params.glms(ivar).conds2resp{icond};
                        end
                    end
                case 'SongType'
                    if sum(strcmp(params.glms(ivar).conds{icond},{'ss_liking_rating','fw_liking_rating','bp_liking_rating','ss_fam_rating','fw_fam_rating','bp_fam_rating'}))>0
                        conds2include{end+1} = params.glms(ivar).conds{icond};
                    else
                        if length(outTable.scanOnset(outTable.SongType==params.glms(ivar).conds{icond}))>0
                            conds2include{end+1} = params.glms(ivar).conds{icond};
                        end
                    end
            end %switch params.glms(ivar).name
        end %for icond=1:length(params.glms(ivar).conds)

        for icond=1:length(conds2include) %loop over each condition in variable/model
            designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).name = conds2include{icond}; %'self-select';
            designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).tmod = params.fmri.tmod;
            designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).pmod = params.fmri.pmod;
            designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).orth = params.fmri.orth;
            if strmatch(conds2include{icond},'q&a','exact')
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur+params.fmri.FAMratingsdur;
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset;
            elseif strmatch(conds2include{icond},'ss_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.SongType=='ss_listening');
            elseif strmatch(conds2include{icond},'fw_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.SongType=='fw_listening');
            elseif strmatch(conds2include{icond},'bp_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.SongType=='bp_listening');
            elseif strmatch(conds2include{icond},'ss_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.SongType=='ss_listening');
            elseif strmatch(conds2include{icond},'fw_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.SongType=='fw_listening');
            elseif strmatch(conds2include{icond},'bp_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.SongType=='bp_listening');
            
            elseif strmatch(conds2include{icond},'hate_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Liking==1);
            elseif strmatch(conds2include{icond},'neutral_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Liking==2);
            elseif strmatch(conds2include{icond},'like_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Liking==3);
            elseif strmatch(conds2include{icond},'love_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Liking==4);
            elseif strmatch(conds2include{icond},'hate_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Liking==1);
            elseif strmatch(conds2include{icond},'neutral_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Liking==2);
            elseif strmatch(conds2include{icond},'like_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Liking==3);
            elseif strmatch(conds2include{icond},'love_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Liking==4);

            elseif strmatch(conds2include{icond},'vunf_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Familiarity==1);
            elseif strmatch(conds2include{icond},'unf_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Familiarity==2);
            elseif strmatch(conds2include{icond},'fam_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Familiarity==3);
            elseif strmatch(conds2include{icond},'vfam_liking_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.LIKratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQALIKOnset(outTable.Familiarity==4);
            elseif strmatch(conds2include{icond},'vunf_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Familiarity==1);
            elseif strmatch(conds2include{icond},'unf_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Familiarity==2);
            elseif strmatch(conds2include{icond},'fam_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Familiarity==3);
            elseif strmatch(conds2include{icond},'vfam_fam_rating','exact')
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.FAMratingsdur;
                 designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanQAFAMOnset(outTable.Familiarity==4);

            else
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).duration = params.fmri.musdur;
                switch params.glms(ivar).name %handle onsets diff for each variable/model
                    case 'Liking'
                        designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanOnset(outTable.Liking==conds2resp2include{icond});
                    case 'Familiarity'
                        designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanOnset(outTable.Familiarity==conds2resp2include{icond});
                    case 'SongType'
                        designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(icond).onset = outTable.scanOnset(outTable.SongType==conds2include{icond});
                 
                end %switch params.glms(ivar).name
            end %if strcmp(conds2include{icond},'q&a')
        end %icond=1:length(params.glms(ivar).conds) 
        
        % now handle case where we have missing response (e.g., 0 for liking or familiarity)
        if ~strcmp(params.glms(ivar).name,'SongType') && length(outTable.scanOnset(outTable.(sprintf('%s',params.glms(ivar).name))==0))>0
            for imiss=1:length(outTable.scanOnset(outTable.(sprintf('%s',params.glms(ivar).name))==0))
                condidx = size(designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond,2);
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(condidx+1).name = sprintf('NA%d',imiss); %'';
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(condidx+1).duration = params.fmri.musdur+params.fmri.LIKratingsdur+params.fmri.FAMratingsdur;
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(condidx+1).tmod = params.fmri.tmod;
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(condidx+1).pmod = params.fmri.pmod;
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(condidx+1).orth = params.fmri.orth;
                %find imiss (nth) index in table so we can add the onset
                tmpidx = find(outTable.(sprintf('%s',params.glms(ivar).name))==0,imiss); 
                designstats.matlabbatch{ivar}.spm.stats.fmri_spec.sess(csess).cond(condidx+1).onset = [outTable.scanOnset(tmpidx(end))];
            end %for imiss=1:length(outTable.scanOnset(outTable.Liking==0))
        end % 
    end %for ivar=1:nvars
   
    %save out table with subXsess responses and onsets
    writetable(outTable,fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),sprintf('taskData_%s_%s-sub-%s.csv',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess})));

    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles = spm_jobman('run',designstats.matlabbatch);
    
    % save out DMs as pdfs
    try
        for iDM=1:length(outFiles)
            [cpath,~,~] = fileparts(outFiles{iDM}.spmmat); %grab the contrast name
            capthtmp = extractBefore(cpath(end:-1:1),"/");
            load([outFiles{iDM}.spmmat{:}]); %load mat
            spm_DesRep('DesMtx',SPM.xX) %grab design matrix
            saveas(gcf, fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),capthtmp(end:-1:1),sprintf('DM_%s_%s-sub-%s.pdf',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess})),'pdf')
            close all
            clear SPM
        end
    catch
        fprintf('WANRING  could not find design matrix!!!\n')
    end

end
