function [outFiles] = spm_contrasts_first_level(isub,isess,inFilesEst,params)
%{
    BMK - 01May2024

adjusts contrasts for missing levels if possible. Note, ignores linear
liking/familiarity if less than 3 levels and ignores quadratic if missing
any values (less than 4 levels)

%}

    % set up batch job
    spm('defaults','fmri');
    spm_jobman('initcfg');
    contrast = struct;
    
    % handle each condition in the model
    try
        nvars = size(inFilesEst,2); 
        fprintf('\n\n~~~working on contrasts for %d models...~~~\n\n',nvars)
    catch
        fprintf('WARNING no 1st level models found in paramter file!!!\n')
        return
    end

    % load up the fmri log to see if all levels exist for a contrast
    ctab = readtable(fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),sprintf('taskData_%s_%s-sub-%s.csv',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess})));
    
    for ivar=1:nvars %loop through each model
        if exist(inFilesEst{ivar}.spmmat{:},'file')~=2 
            fprintf('WANRING  could not find SPM.mat!!!\n')
            return
        end
        currglm = params.glms(ivar); %transfer over contrast/model info temp;
        ctrials = unique(ctab.(sprintf('%s',currglm.name))); %grab all resp/levels for this var
        if ~strcmp(currglm.name,'SongType')
            ctrials(find(ctrials==0)) = []; %remove 0 if missing trial present
        end
        % handle missing response levels
        matches = regexp(currglm.conds,'_rating','match');
        if length(ctrials) < sum(cellfun(@isempty,matches))
            if length(ctrials)<=1 || ~params.adjust4missingContrasts
                %not enough trials for our contrasts
                fprintf('\nNOT ENOUGH LEVELS/TRIALS FOR SUBJECT: %s; SESSION: %s-%s MODEL: %s\n\n',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess},currglm.name)
                continue
            end
            if length(ctrials)>1
                %have more than 1 level, figure out the right contrast
                fprintf('\nmissing level(s) for SUBJECT: %s; SESSION: %s-%s MODEL: %s\n\nadjusting...\n',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess},currglm.name)
                % fix the contrasts (update F and t conts seperately incase they don't match)
                for icont=1:length(currglm.Fconts)
                    a=1;
                end %icont=1:length(currglm.Fconts)
                for icont=1:length(currglm.Tconts)
                    switch currglm.Tconts{icont} %t conts  
                        case {'T_Lall_listening_>_all_rating','T_Fall_listening_>_all_rating'}
                            if (~ismember(1,ctrials) || ~ismember(2,ctrials)) && (ismember(1,ctrials) || ismember(2,ctrials)) && (ismember(3,ctrials) && ismember(4,ctrials)) 
                                %if 1 or 2 is missing, but not both, and we have 3 and 4
                                currglm.Tweights{icont} = [2 2 2 -1 -1 -1 -1 -1 -1 0];  %[2 2 2 2 -1 -1 -1 -1 -1 -1 -1 -1]
                            elseif (~ismember(3,ctrials) || ~ismember(4,ctrials)) && (ismember(3,ctrials) || ismember(4,ctrials)) && (ismember(1,ctrials) && ismember(2,ctrials)) 
                                %if 3 or 4 is missing, but not both, and we have 1 and 2
                                currglm.Tweights{icont} = [2 2 2 -1 -1 -1 -1 -1 -1 0];
                            elseif (length(unique(ctrials))==2)
                                currglm.Tweights{icont} = [2 2 -1 -1 -1 -1 0];
                            elseif (length(unique(ctrials))==1)
                                currglm.Tweights{icont} = [2 -1 -1 0];
                            else %remove if no levels (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                            end 
                        case {'T_linear_liking_listening','T_linear_fam_listening'}
                            if (~ismember(1,ctrials) || ~ismember(2,ctrials)) && (ismember(1,ctrials) || ismember(2,ctrials)) && (ismember(3,ctrials) && ismember(4,ctrials)) %if 1 or 2 is missing, but not both, and we have 3 and 4
                                currglm.Tweights{icont} = [-4 1 3 0 0 0 0 0 0 0]; 
                            elseif (~ismember(3,ctrials) || ~ismember(4,ctrials)) && (ismember(3,ctrials) || ismember(4,ctrials)) && (ismember(1,ctrials) && ismember(2,ctrials)) %if 3 or 4 is missing, but not both, and we have 1 and 2
                                currglm.Tweights{icont} = [-3 -1 4 0 0 0 0 0 0 0];
                            else %remove if only two levels (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                            end 
                        
                        case {'T_linear_liking_liking_rating','T_linear_fam_liking_rating'}
                            if (~ismember(1,ctrials) || ~ismember(2,ctrials)) && (ismember(1,ctrials) || ismember(2,ctrials)) && (ismember(3,ctrials) && ismember(4,ctrials)) %if 1 or 2 is missing, but not both, and we have 3 and 4
                                currglm.Tweights{icont} = [0 0 0 -4 1 3 0 0 0 0]; 
                            elseif (~ismember(3,ctrials) || ~ismember(4,ctrials)) && (ismember(3,ctrials) || ismember(4,ctrials)) && (ismember(1,ctrials) && ismember(2,ctrials)) %if 3 or 4 is missing, but not both, and we have 1 and 2
                                currglm.Tweights{icont} = [0 0 0 -3 -1 4 0 0 0 0];
                            else %remove if only two levels (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                            end 
                        case {'T_linear_liking_fam_rating','T_linear_fam_fam_rating'}
                            if (~ismember(1,ctrials) || ~ismember(2,ctrials)) && (ismember(1,ctrials) || ismember(2,ctrials)) && (ismember(3,ctrials) && ismember(4,ctrials)) %if 1 or 2 is missing, but not both, and we have 3 and 4
                                currglm.Tweights{icont} = [0 0 0 0 0 0 -4 1 3 0]; 
                            elseif (~ismember(3,ctrials) || ~ismember(4,ctrials)) && (ismember(3,ctrials) || ismember(4,ctrials)) && (ismember(1,ctrials) && ismember(2,ctrials)) %if 3 or 4 is missing, but not both, and we have 1 and 2
                                currglm.Tweights{icont} = [0 0 0 0 0 0 -3 -1 4 0];
                            else %remove if only two levels (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                            end 
                        case {'T_linear_liking_all_rating','T_linear_fam_all_rating'}
                            if (~ismember(1,ctrials) || ~ismember(2,ctrials)) && (ismember(1,ctrials) || ismember(2,ctrials)) && (ismember(3,ctrials) && ismember(4,ctrials)) %if 1 or 2 is missing, but not both, and we have 3 and 4
                                currglm.Tweights{icont} = [0 0 0 -4 1 3 -4 1 3 0]; 
                            elseif (~ismember(3,ctrials) || ~ismember(4,ctrials)) && (ismember(3,ctrials) || ismember(4,ctrials)) && (ismember(1,ctrials) && ismember(2,ctrials)) %if 3 or 4 is missing, but not both, and we have 1 and 2
                                currglm.Tweights{icont} = [0 0 0 -3 -1 4 -3 -1 4 0];
                            else %remove if only two levels (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                            end 

                        case {'T_hate_listening_>_hate_rating','T_vunf_listening_>_vunf_rating'}
                             if ismember(1,ctrials) && sum(ismember(4,ctrials) + ismember(2,ctrials) + ismember(3,ctrials)==2) %we have hate and 2 more
                                currglm.Tweights{icont} = [2 0 0 -1 0 0 -1 0 0 0]; %[2 0 0 0 -1 0 0 0 -1 0 0 0]
                             elseif ismember(1,ctrials) && sum(ismember(4,ctrials) + ismember(2,ctrials) + ismember(3,ctrials)==1) %we hate love and 1 more
                                 currglm.Tweights{icont} = [2 0 -1 0 -1 0 0];
                             elseif ismember(1,ctrials) && sum(ismember(4,ctrials) + ismember(2,ctrials) + ismember(3,ctrials)==0)
                                 currglm.Tweights{icont} = [2 -1 -1 0];
                             else %remove if only two levels from same category (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                             end
                        case {'T_neutral_listening_>_neutral_rating','T_unf_listening_>_unf_rating'}
                             if ismember(2,ctrials) && ismember(1,ctrials)  && sum(ismember(4,ctrials) + ismember(3,ctrials)==1) %we have neutral and hate and 1 more
                                currglm.Tweights{icont} = [0 2 0 0 -1 0 0 -1 0 0]; 
                             elseif ismember(2,ctrials) && ~ismember(1,ctrials)  && sum(ismember(4,ctrials) + ismember(3,ctrials)==2) %we have neutral and not hate but 2 more
                                 currglm.Tweights{icont} = [2 0 0 -1 0 0 -1 0 0 0];
                             elseif ismember(2,ctrials) && ismember(1,ctrials)  && sum(ismember(4,ctrials) + ismember(3,ctrials)==0) %we have neutral and hate only
                                 currglm.Tweights{icont} = [0 2 0 -1 0 -1 0];
                             elseif ismember(2,ctrials) && ~ismember(1,ctrials)  && sum(ismember(4,ctrials) + ismember(3,ctrials)==1) %we have neutral and not hate but 1 more
                                 currglm.Tweights{icont} = [2 0 -1 0 -1 0 0];
                             elseif ismember(2,ctrials) && ~ismember(1,ctrials)  && sum(ismember(4,ctrials) + ismember(3,ctrials)==0) %we have neutral only
                                 currglm.Tweights{icont} = [2 -1 -1 0];
                             else %remove if only two levels from same category (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                             end
                        case {'T_like_listening_>_like_rating','T_fam_listening_>_fam_rating'}
                             if ismember(1,ctrials) && ismember(2,ctrials) && ismember(3,ctrials) && ~ismember(4,ctrials) %we have like and hate and 1 between
                                currglm.Tweights{icont} = [0 0 2 0 0 -1 0 0 -1 0]; 
                             elseif ~ismember(1,ctrials) && ismember(2,ctrials) && ismember(3,ctrials) && ismember(4,ctrials) %we have like and not hate but 2 more
                                 currglm.Tweights{icont} = [0 2 0 0 -1 0 0 -1 0 0];
                             elseif ismember(1,ctrials) && ~ismember(2,ctrials) && ismember(3,ctrials) && ismember(4,ctrials) %we have hate like and love
                                 currglm.Tweights{icont} = [0 2 0 0 -1 0 0 -1 0 0];
                             elseif (ismember(1,ctrials) || ismember(2,ctrials)) && ismember(3,ctrials) && ~ismember(4,ctrials) %we have like and hate only or neutral
                                 currglm.Tweights{icont} = [0 2 0 -1 0 -1 0];
                             elseif ~ismember(1,ctrials) && ~ismember(2,ctrials) && ismember(3,ctrials) && ismember(4,ctrials) %we have like and not hate or (neutral)
                                 currglm.Tweights{icont} = [2 0 -1 0 -1 0 0];
                             elseif ~ismember(1,ctrials) && ~ismember(2,ctrials) && ismember(3,ctrials) && ~ismember(4,ctrials) %we have neutral only
                                 currglm.Tweights{icont} = [2 -1 -1 0];
                             else %remove if only two levels from same category (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                             end
                         case {'T_love_listening_>_love_rating','T_vfam_listening_>_vfam_rating'}
                             if ismember(4,ctrials) && sum(ismember(1,ctrials) + ismember(2,ctrials) + ismember(3,ctrials)==2) %we have love and 2 more
                                currglm.Tweights{icont} = [0 0 2 0 0 -1 0 0 -1 0]; %
                             elseif ismember(4,ctrials) && sum(ismember(1,ctrials) + ismember(2,ctrials) + ismember(3,ctrials)==1) %we have  and 1 more
                                 currglm.Tweights{icont} = [0 2 0 -1 0 -1 0];
                             elseif ismember(4,ctrials) && sum(ismember(1,ctrials) + ismember(2,ctrials) + ismember(3,ctrials)==0)
                                 currglm.Tweights{icont} = [2 -1 -1 0];
                             else %remove if only two levels from same category (set to zero)
                                currglm.Tweights{icont} = [];
                                currglm.Tconts{icont} = [];
                                currglm.Tsessrep{icont} = [];
                             end 

                        
                    end %switch currglm.Tconts{icont}
                end %icont=1:length(currglm.Tconts)
            end %length(ctrials)==3
        end %if length(ctrials) < size(currglm.conds,2)
        
        if length(ctrials) > size(currglm.conds,2)
            fprintf('\n\nWARNING found more response levels than expected for sub:\n%s\nIgnoring trial!!!\n',sprintf('%s_%s-sub-%s',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess}))
        end

        % add some params
        contrast.matlabbatch{ivar}.spm.stats.con.spmmat = {[inFilesEst{ivar}.spmmat{:}]};
        contrast.matlabbatch{ivar}.spm.stats.con.delete = params.fmri.stats.con.delete;
        
        % add the contrasts
        nconts = 0;
        for icont=1:length(currglm.Fconts) %loop through each F contrasts
            if ~isempty(currglm.Fconts{icont})
                nconts = nconts + 1;
                contrast.matlabbatch{ivar}.spm.stats.con.consess{nconts}.fcon.name = currglm.Fconts{icont};
                contrast.matlabbatch{ivar}.spm.stats.con.consess{nconts}.fcon.weights = currglm.Fweights{icont};
                contrast.matlabbatch{ivar}.spm.stats.con.consess{nconts}.fcon.sessrep = currglm.Fsessrep{icont};
            end
        end %icont=1:length(currglm.Fconts)
        for icont=1:length(currglm.Tconts) %loop through each T contrasts
            if ~isempty(currglm.Tconts{icont})
                nconts = nconts + 1;
                contrast.matlabbatch{ivar}.spm.stats.con.consess{nconts}.tcon.name = currglm.Tconts{icont};
                contrast.matlabbatch{ivar}.spm.stats.con.consess{nconts}.tcon.weights = currglm.Tweights{icont};
                contrast.matlabbatch{ivar}.spm.stats.con.consess{nconts}.tcon.sessrep = currglm.Tsessrep{icont};
            end
        end %icont=1:length(currglm.Tconts)
    end %for ivar=1:nvars
    
    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles = spm_jobman('run',contrast.matlabbatch);

    % create a table of all contrasts run for this subject 
    varTypes = ["string" "string" "string" repelem(["string"],length(params.allContrasts))];
    varNames = ["ID" "NUM" "session" params.allContrasts]; 
    outTable = table('Size',[1 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
    try
        for iDM=1:length(outFiles)
            load([outFiles{iDM}.spmmat{:}]); %load mat
            outTable(1,1:3) = {isub.subjectID,isub.subjectNUM{isess},isub.sessNames{isess},}; %update table that will hold all are onset and response info
            for icont=1:size(SPM.xCon,2) %add all contrast codes
                cstrng = sprintf('['); 
                for itmp=1:length(SPM.xCon(icont).c)
                    cstrng = [cstrng sprintf('%d ',SPM.xCon(icont).c(itmp))];
                end
                cstrng = [cstrng(1:end-1) ']'];
                outTable(1,find(ismember(outTable.Properties.VariableNames,SPM.xCon(icont).name)==1)) = {cstrng}; 
            end
            clear SPM
        end %iDM=1:length(outFiles)  
    catch
        fprintf('WANRING  could not find contrast(s)!!!\n')
    end
    %save out table with contrasts
    writetable(outTable,fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),sprintf('contrasts_%s_%s-sub-%s.csv',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess})));
end