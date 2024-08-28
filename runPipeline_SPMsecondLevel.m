function runPipeline_SPMsecondLevel(slurmIN)
    %{
    BMK - 05May2024
    
    INPUTS:
    *params_gMBI.m - make changes here (paths, subject info, etc.)
    *first level spm models (from runPipline_SPMfirstLevel.m)
    
    OUTPUTS:
    
    
    %}
    warning('off')
    
    params_gMBI; %load params
    level2conts = {'preXpost_mci','pre_healthy','preVpost','healthyVmciPre','healthyVmciPost','betas_pre&post_mci'}; %2 level contrast options  
    % 

    % handle slurm inputs if needed
    try
        conts2use = level2conts(slurmIN);
    catch
        conts2use = level2conts; %run all subs in params if none specified by SLURM
    end
    
    fprintf('\n\n\n\n~~~~~~~~~~running SPM Second Level pipeline for %d contrasts~~~~~~~~~~\n\n\n\n',length(conts2use))
    subjectInfo %print info to terminal/logfile
    conts2use
    fprintf('\n\n\n') 

    startdir = pwd;

    for icont=1:length(conts2use)
        switch conts2use{icont}
            case 'preXpost_mci'
                % build, estimate, contrasts for second level pre and post models (1 sample t for each session)
                outFilesCnt1 = spm_preXpost_second_level({'gmbi'},{'MCI'},{'gMBI'},{'pre','post'},subjectInfo,params); %batch, MCIstatus, interventionType, sessnames
                % visual contrast results
                spm_visualize_contrasts(outFilesCnt1,params);
                cd(startdir)
            
            case 'pre_healthy'
                % build, estimate, contrasts for second level pre and post models (1 sample t for each session)
                outFilesCnt1 = spm_preXpost_second_level({'gmbi','healthy','mbi'},{'healthy'},{'gMBI','MBI','NA'},{'pre'},subjectInfo,params); %batch, MCIstatus, interventionType, sessnames
                % visual contrast results
                spm_visualize_contrasts(outFilesCnt1,params);
                cd(startdir)
                
            case 'preVpost'
                % build, estimate, contrasts for second level pre vs post models (paired t)
                outFilesCnt2 = spm_preVpost_second_level('gmbi','gMBI',subjectInfo,params); %batch, interventionType
                % visual contrast results
                spm_visualize_contrasts(outFilesCnt2,params);
                cd(startdir)
            
            case 'healthyVmciPre'
                % build, estimate, contrasts for second level healthy vs mci models (2 sample t for pre)
                outFilesCnt3 = spm_healthyVmci_second_level({'gmbi','mbi'},{'healthy','MCI'},{1,1},subjectInfo,params); %batches, MCIstatus, sessions (healthy->1,MCI->1/2)
                % visual contrast results  
                spm_visualize_contrasts(outFilesCnt3,params);
                cd(startdir)
            
            case 'healthyVmciPost'
                % build, estimate, contrasts for second level healthy vs mci models (2 sample t, pre for healthy post for mci)
                outFilesCnt4 = spm_healthyVmci_second_level({'gmbi','mbi'},{'healthy','MCI'},{1,2},subjectInfo,params); %batches, MCIstatus, sessions (healthy->1,MCI->1/2)
                % visual contrast results  
                spm_visualize_contrasts(outFilesCnt4,params);
                cd(startdir)
             
            case 'betas_pre&post_mci'
                % collect and plot subject-level betas from Contrasts of Interest 
                outFilesCnt5 = spm_betas_prepost_second_level({'gmbi'},{'healthy','MCI'},{'gMBI'},subjectInfo,params); %batch, MCIstatus, interventionType
            
            case 'betas_pre&post_healthy'
                % collect and plot subject-level betas from Contrasts of Interest 
                outFilesCnt6 = spm_betas_prepost_second_level({'gmbi','mbi'},{'healthy'},{'gMBI','MBI','NA'},subjectInfo,params); %batch, MCIstatus, interventionType


        end %switch
    
        % 
    end %for imod=1:length(2ndLevelModels)
    
    fprintf('\n\n~~~~~finished SPM second level models for %d contrasts~~~~~\n\n',length(conts2use))
    
end
