function runPipeline_SPMfirstLevel(slurmIN)
    %{
    BMK - 01May2024
    
    INPUTS:
    *params_gMBI.m - make changes here (paths, subject info, etc.)
    *uses mus-bid log files
    *uses preprocessed 4d mus-bid niftis
    
    OUTPUTS:
    *taskData*.csv - table with event order,onsets, and responses (in sub dirs)
    *contrasts*.csv - table with all contrasts evaluated (in sub dirs)
    *SPM.mat - constructed 1st level models (in sub/model dirs)
    *DM*.pdf - design matrix (in sub/model dirs)
    
    %}
    warning('off')
    
    params_gMBI; %load params
    
    % handle slurm inputs if needed
    try
        sub2use = subjectInfo(slurmIN);
    catch
        sub2use = subjectInfo; %run all subs in params if none specified by SLURM
    end
    
    fprintf('\n\n\n\n~~~~~~~~~~running SPM First Level pipeline for %d subjects~~~~~~~~~~\n\n\n\n',length({sub2use.subjectID}))
    sub2use %print info to terminal/logfile
    fprintf('\n\n\n')

    % Loop through subjects listed in params_gMBI.m
    for isub=1:length({sub2use.subjectID})
        for isess=1:length(sub2use(isub).sessNames) 
    
            % build first level models 
            %outFilesDM = spm_build_first_level(sub2use(isub),isess,params);
            
            % estimate first level models
            %outFilesEst = spm_estimate_first_level(outFilesDM,params);
            
            % evaluate contrasts for first level models
            %outFilesCnt = spm_contrasts_first_level(sub2use(isub),isess,outFilesEst,params);
            
            %fprintf('\n\n~~~~~finished SPM first level for subject: %s; session: %d - sub-%s;~~~~~\n\n',sub2use(isub).subjectID,isess,sub2use(isub).subjectNUM{isess})
    
            % extract betas for conditions and ROIs of interest 
            COIs = {'T_linear_liking_listening','T_linear_fam_listening','T_ss_listening_>_fw_listening','T_ss_listening_>_os_listening'};
            outTableBs = spm_grabBetas_first_level(COIs,sub2use(isub),isess,params);
            
        end %for isess=1:length(subjectInfo(isub).sessNames)
        
    end %for isub=1:length({subjectInfo.subjectID})
end
