function [outFiles] = spm_estimate_first_level(inFilesDM,params)
%{
    BMK - 01May2024

%}

    % set up batch job
    spm('defaults','fmri');
    spm_jobman('initcfg');
    modelestimation = struct;
    
    % handle each condition in the model
    try
        nvars = size(inFilesDM,2); 
        fprintf('\n\n~~~working on estimating %d models...~~~\n\n',nvars)
    catch
        fprintf('WARNING no 1st level models found in paramter file!!!\n')
        return
    end
    
    for ivar=1:nvars %loop through each model
        if exist(inFilesDM{ivar}.spmmat{:},'file')~=2 
            fprintf('WANRING  could not find SPM.mat!!!\n')
            return
        end
        modelestimation.matlabbatch{ivar}.spm.stats.fmri_est.spmmat = {[inFilesDM{ivar}.spmmat{:}]};
        modelestimation.matlabbatch{ivar}.spm.stats.fmri_est.write_residuals = params.fmri.fmri_est.write_residuals;
        modelestimation.matlabbatch{ivar}.spm.stats.fmri_est.method.Classical = params.fmri.fmri_est.method.Classical;
    end %for ivar=1:nvars
    
    % run the job (1 sub, 1 session, params.glms(*) number of models)
    outFiles = spm_jobman('run',modelestimation.matlabbatch);

end