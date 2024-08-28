function [outTable] = spm_grabBetas_first_level(COIs,isub,isess,params)
    %{
        BMK - 08Jul2024
    script loads up an estimated GLM and extraxt
    %}
    
    fprintf('\n\n\n~~~~~Grabbing Betas for subject: %s; session: %d; sub-%s;~~~~~\n\n\n',isub.subjectID,isess,isub.subjectNUM{isess})
    
    % handle each condition in the model
    try
        nvars = size(params.glms,2); 
        fprintf('working on %d models...\n',nvars)
    catch
        fprintf('WARNING no variables/conditions found in paramter file!!!\n')
        return
    end

    % get subset of ROIs we want to examine.
    roi_filtlist = struct();
    tmpcounter = 1;
    for iROI=1:length({params.ROIs.roi_files.name})
        cmatches = regexp({params.ROIs.roi_files(iROI).name},params.ROIs.rois2include,'match');
        if sum(~cellfun(@isempty,cmatches))==1
            if tmpcounter==1
                roi_filtlist = params.ROIs.roi_files(iROI);
            else
                roi_filtlist(tmpcounter) = params.ROIs.roi_files(iROI);
            end
            tmpcounter = tmpcounter + 1;
        end
    end

    % create a table of all contrasts run for this subject 
    varTypes = ["string" "string" "string" "string" repelem(["double"],length({roi_filtlist.name}))];
    varNames = ["ID" "NUM" "session" "contrast" {roi_filtlist.name}]; 
    outTable = table('Size',[0 length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
    for ivar=1:nvars %loop over each variable/model
        % save dir (one for each subXmodel
        savepath = fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),params.glms(ivar).name); 
        if exist(fullfile(savepath,'SPM.mat'),'file')~=2 %remove existing SPM.mat if present to avoid GUI prompt
            fprintf(['WANRING: attempting to grab betas from non existent model ' savepath '!!!\n']);
            continue
        end
        % find first level SPM GLM (one for each subXmodel)
        SPMin = dir(fullfile(savepath,'SPM.mat'));

        if ~isempty(SPMin) && length({SPMin.name})==1
            for icont=1:length(params.glms(ivar).Tconts) %loop over each condition in variable/model

                %if it's a contrast we are interested in 
                if length(strmatch(params.glms(ivar).Tconts{icont},COIs,'exact'))>0
                    fprintf('Working on contrast: %s\n',params.glms(ivar).Tconts{icont})

                    %find the correct beta file
                    cSPM = load(fullfile(SPMin.folder,SPMin.name));
                    cSPM.SPM.xCon(1).name
                    
                    matches = regexp({cSPM.SPM.xCon.name},params.glms(ivar).Tconts{icont},'match');
                    conIdx = find(~cellfun(@isempty,matches)==1);
                    if ~isempty(conIdx) && length(conIdx)==1
                        %update the out table
                        ctableidx = size(outTable,1)+1;
                        outTable(ctableidx,1:3) = {isub.subjectID,isub.subjectNUM{isess},isub.sessNames{isess},};
                        outTable(ctableidx,find(ismember(outTable.Properties.VariableNames,'contrast')==1)) = {params.glms(ivar).Tconts{icont}};
                        for iROI=1:length(roi_filtlist)
                            %extract betas from con file for contrast and ROI
                            cBeta = spm_summarise(fullfile(SPMin.folder,cSPM.SPM.xCon(conIdx).Vcon.fname),...
                                fullfile(roi_filtlist(iROI).folder,roi_filtlist(iROI).name),@mean);
                            %save out betas
                            outTable(ctableidx,find(ismember(outTable.Properties.VariableNames,roi_filtlist(iROI).name)==1)) = {cBeta}; 
                        end
                    elseif isempty(conIdx)
                        % need to update table with null result if contrast missing
                        %update the out table
                        ctableidx = size(outTable,1)+1;
                        outTable(ctableidx,1:3) = {isub.subjectID,isub.subjectNUM{isess},isub.sessNames{isess},};
                        outTable(ctableidx,find(ismember(outTable.Properties.VariableNames,'contrast')==1)) = {params.glms(ivar).Tconts{icont}};
                       
                    else
                        fprintf('more than 1 match')
                    end
                else
                    fprintf('Skipping contrast: %s\n',params.glms(ivar).Tconts{icont})
                end

            end %for icont=1:length(params.glms(ivar).Tconts)
            
        else
            fprintf('WARNING no first level SPM.mat found for subject: %s; session: %d; sub-%s; skipping!!!\n',isub.subjectID,isess,isub.subjectNUM{isess})
            continue
        end % if ~isempty(SPMin) && length({SPMin.name})==1

    end %for ivar=1:nvars
    
    %save out table with betas
    writetable(outTable,fullfile(params.paths.spm1stLevel,sprintf('sub-%s',isub.subjectNUM{isess}),sprintf('TconBetasOfInterest_%s_%s-sub-%s.csv',isub.subjectID,isub.sessNames{isess},isub.subjectNUM{isess})));


end
