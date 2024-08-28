function [outImg] = spm_visualize_contrasts(inConts,params)
%{
BMK - 06May2024
based on code from https://github.com/bramzandbelt/slice_display
%}
    
   fprintf('\n\n\n~~~~~Visualizing %d contrasts from %d model(s);~~~~~\n\n\n',length(inConts),length({params.glms.name}))
    
   % set up batch job for thresholding in SPM
    spm('defaults','fmri');
    spm_jobman('initcfg');
    
    % Get custom colormaps
    load(fullfile(params.paths.display,'colormaps.mat'));
    

    for imod=1:length(inConts) %loop throguh each model %DONT WORK: not enough subs? 2 (quad liking_missing 9/11?); 3 (T_allInto_V_allMeh missing none); 7 (T_quad_familiarity missing 6);
        cSPM = load(inConts{imod}.spmmat{:});
        [cpath,cname,cext] = fileparts(inConts{imod}.spmmat{:});
        
        for icont=1:size(cSPM.SPM.xCon,2) %loop through each contrast
            % run some stat thresholding before ploting
            thresh = struct;
            thresh.matlabbatch{1}.spm.stats.results.spmmat = inConts{imod}.spmmat;
            thresh.matlabbatch{1}.spm.stats.results.conspec = params.fmri.conspec;
            thresh.matlabbatch{1}.spm.stats.results.conspec.titlestr = cSPM.SPM.xCon(icont).name;
            thresh.matlabbatch{1}.spm.stats.results.units = params.fmri.results.units;
            thresh.matlabbatch{1}.spm.stats.results.print = params.fmri.results.print;
            thresh.matlabbatch{1}.spm.stats.results.write.none = params.fmri.results.write.none; 
            thresh.matlabbatch{1}.spm.stats.results.export{1}.binary.basename = [cSPM.SPM.xCon(icont).name,'_thrs'];
            try
                % run the 2nd level threshold
                outFiles = spm_jobman('run',thresh.matlabbatch);
                close all
                tmpt = niftiread(fullfile(cpath,sprintf('spmF_%04d_%s_thrs.nii',1,cSPM.SPM.xCon(icont).name)));
                if sum(sum(sum(tmpt)))>0
                
                    % make/save out pdf
                    layers = sd_config_layers('init',{'truecolor','blob'}); %Initialize empty layers
                    settings = sd_config_settings('init');
                    %define layers           
                    layers(1).color.file = fullfile(spm('Dir'),'canonical','single_subj_T1.nii'); %Layer 1: anatomical image
                    layers(1).color.map = gray(256);
                    layers(2).color.file = fullfile(cpath,cSPM.SPM.xCon(icont).Vcon.fname); %Layer 2: thresholded t-map
                    layers(2).color.map = CyBuBkRdYl;
                    layers(2).color.label = [cSPM.SPM.xCon(icont).STAT '-value'];
                    layers(2).mask.file  = fullfile(cpath,sprintf('spmF_%04d_%s_thrs.nii',1,cSPM.SPM.xCon(icont).name));
                    %specify settings
                    settings.slice.orientation = 'sagittal'; %axial; coronal
                    settings.slice.disp_slices = -70:4:70;
                    settings.fig_specs.n.slice_column = 6;
                    settings.fig_specs.title = sprintf('%s; height: %1.3f; extent: %d; nScans: %d',strrep(cSPM.SPM.xCon(icont).name,'_',' '),params.fmri.conspec.thresh,params.fmri.conspec.extent,cSPM.SPM.nscan); %remove _ replace with space
                    %display/save
                    sd_display(layers,settings);
                    tmpsag = fullfile(cpath,sprintf('%s_sag_thrs.pdf',cSPM.SPM.xCon(icont).name));
                    saveas(gcf, tmpsag,'pdf'); %save current figure
                    close all
        
                    % copy previous and save axial
                    settings.slice.orientation = 'axial'; 
                    settings.slice.disp_slices = -56:4:84;
                    %display/save
                    sd_display(layers,settings);
                    tmphor = fullfile(cpath,sprintf('%s_hor_thrs.pdf',cSPM.SPM.xCon(icont).name));
                    saveas(gcf, tmphor,'pdf'); %save current figure
                    close all
        
                    % copy previous and save coronal
                    settings.slice.orientation = 'coronal'; 
                    settings.slice.disp_slices = -94:4:70;
                    %display/save
                    sd_display(layers,settings);
                    tmpcor = fullfile(cpath,sprintf('%s_cor_thrs.pdf',cSPM.SPM.xCon(icont).name));
                    saveas(gcf, tmpcor,'pdf'); %save current figure
                    close all
                else
                    fprintf("NO VOXELS FOUND FOR %s with height thresh: %1.3f; extent thresh: %d\n",cSPM.SPM.xCon(icont).name,params.fmri.conspec.thresh,params.fmri.conspec.extent)
                end %if sum(sum(sum(tmpt)))>0

                % append_pdfs and delete ind. copies
                %append_pdfs(fullfile(cpath,sprintf('%s_thrs.pdf',cSPM.SPM.xCon(icont).name)),[tmpsag,tmphor,tmpcor])   
            catch
                fprintf("ERROR, COULD NOT DISPLAY CONTRAST!!!\n")
            end

        end %for icont=1:size(cSPM.SPM.xCon,2)
    end %for imod=1:length(inConts)

end