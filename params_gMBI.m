%params_gMBI.m
%{
BMK - 01May2024

Parameters for gamma MBI fMRI mus-bid task preprocessing+pipeline

%}

%% path parameters (CHANGE THESE)
user = 'xi.wu'; %name of home folder on Discovery
homedir = sprintf('/scratch/%s/SPM_SMPC/',user);
params.paths.data = '/work/mindlab/Projects/GammaMBI/MRI/Preprocessed_Data_Conn/'; %mus-bid .nii files must be here!
params.paths.logs = fullfile(homedir,'data/mus-bid_task_logs'); %mus-bid task log files must be here!
params.paths.derivatives = fullfile(homedir,'data/derivatives');
params.paths.spm1stLevel = fullfile(params.paths.derivatives,'SPM/1stLevel');
params.paths.spm2ndLevel = fullfile(params.paths.derivatives,'SPM/2ndLevel');
params.paths.connProject = fullfile(params.paths.derivatives,'CONN/');

params.paths.homedir = homedir;

% add SPM to the path
params.paths.spmdir = fullfile('/work/mindlab/Programs/spm12/');
addpath(params.paths.spmdir);

% add CONN to the path
params.paths.conn = fullfile('/work/mindlab/Programs/conn/');
addpath(params.paths.conn);

% add slice_display to the path
params.paths.display = fullfile(params.paths.homedir,'code','slice_display-master');
addpath(params.paths.display);

% add Panel to the path
params.paths.panel = fullfile(params.paths.homedir,'code','panel');
addpath(params.paths.panel);

% add append_pdfs to the path
params.paths.pdfs = fullfile(params.paths.homedir,'code','append_pdfs');
addpath(params.paths.pdfs);

% add plotSpread to the path
params.paths.pspread = fullfile(params.paths.homedir,'code','plotSpread','plotSpread');
addpath(params.paths.pspread);

%% fMRI params
%~~ First Level ~~%
params.fmri.musdur = 42; %duration of music in scans (42) see timing_units above
params.fmri.LIKratingsdur = 9; %duration of Q&A in scans (18) 
params.fmri.FAMratingsdur = 9;

%see https://www.andysbrainblog.com/andysbrainblog/tag/SPM.mat "What's in an SPM.mat Files?
params.fmri.tmod = 0;
params.fmri.pmod = struct('name', {}, 'param', {}, 'poly', {});
params.fmri.orth = 1;
% fmri_spec
params.fmri.fmri_spec.timing.units = 'scans'; %make sure this matches musdur and ratingsdur below
params.fmri.fmri_spec.timing.RT = 0.475;
params.fmri.fmri_spec.timing.fmri_t = 16;
params.fmri.fmri_spec.timing.fmri_t0 = 8;
params.fmri.fmri_spec.sess.multi = {''};
params.fmri.fmri_spec.sess.regress = struct('name', {}, 'val', {});
params.fmri.fmri_spec.sess.multi_reg = {''};
params.fmri.fmri_spec.sess.hpf = Inf;
params.fmri.fmri_spec.fact = struct('name', {}, 'levels', {});
params.fmri.fmri_spec.bases.hrf.derivs = [0 0];
params.fmri.fmri_spec.volt = 1;
params.fmri.fmri_spec.global = 'None';
params.fmri.fmri_spec.mthresh = 0.8;
params.fmri.fmri_spec.mask = {''};
params.fmri.fmri_spec.cvi = 'AR(1)';
% fmri_est
params.fmri.fmri_est.write_residuals = 0;
params.fmri.fmri_est.method.Classical = 1;
% stats
params.fmri.stats.con.delete = 0;

%~~ Second Level~~ %
% factorial_design
params.fmri.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
params.fmri.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
params.fmri.factorial_design.masking.tm.tm_none = 1;
params.fmri.factorial_design.masking.im = 1;
params.fmri.factorial_design.masking.em = {''};
params.fmri.factorial_design.globalc.g_omit = 1;
params.fmri.factorial_design.globalm.gmsca.gmsca_no = 1;
params.fmri.factorial_design.globalm.glonorm = 1;

params.fmri.factorial_design_paired.des.pt.gmsca = 0;
params.fmri.factorial_design_paired.des.pt.ancova = 0;

params.fmri.factorial_design_2samp.des.t2.dept = 0;
params.fmri.factorial_design_2samp.des.t2.variance = 1;
params.fmri.factorial_design_2samp.des.t2.gmsca = 0;
params.fmri.factorial_design_2samp.des.t2.ancova = 0;

% conspec
params.fmri.conspec.contrasts = 1; %use all existing contrasts?
params.fmri.conspec.threshdesc = 'none'; %family-wise error
params.fmri.conspec.thresh = 0.001; %p-value threshold
params.fmri.conspec.extent = 0; %vox extent
params.fmri.conspec.mask.none = 1; %??
params.fmri.results.units = 1;
params.fmri.results.print = 'pdf';
params.fmri.results.write.none = 1; %write.none WHATS THIS DO?

%~~ CONN ~~ %
params.conn.Setup.analyses = 1; % 1 for ROI-to-ROI analyses, 2 for seed-to-voxel analyses, and 3 for voxel-to-voxel analyses
params.conn.Analysis.type = 1; % 1 for ROI-to-ROI analyses, 2 for seed-to-voxel analyses, and 3 for voxel-to-voxel analyses
%Integer identifying the connectivity measure used. Current valid values are: 
%1 Correlation (bivariate) [default]; 2 Correlation (semipartial); 3 Regression (bivariate); 4 Regression (multivariate)
params.conn.Analysis.measure = 1; 
%Integer identifying the type of intra-block weighting performed: 1 none; 2 hrf [default]; 3 hanning
%By default the scans within each block (or session if no task conditions are defined) 
%is weighted using a hrf-shaped weighting function. This effectively weights down the 
%contribution of the first few scans on the estimation of connectivity measures. 
%Alternatively you can define no weighting (1) or a hanning weighting function (3) that 
%weights more heavily only in the middle scans within each block.
params.conn.Analysis.weight = 2; 
%Specifies ROI sources for connectivity analyses (connectivity measures will be estimated 
%between these sources and every voxel -for seed-to-voxel analyses-, and between these 
%sources and every ROI -for ROI-to-ROI analyses-). If this field is left unset the toolbox 
%will by default use as sources all the ROIs that have not been defined as confounds.
%params.conn.Analysis.sources = {'MPFC','PCC'};
%This field contains a cell array ranging over each source ROI, where the 
%BATCH.Analysis.sources.names{nsource} element contains a char array with the name of the 
%source nsource. Valid source names can be any ROI name (defined in BATCH.Setup.rois.names). 
%Note that each of the sources defined here can be multivariate (defined by multiple time-series, see below)
params.conn.Analysis.sources.names = {'MPFC','PCC'}; %(specifies two ROI sources for the current experiment)
%This field contains a cell array ranging over each source ROI, where the 
%BATCH.Analysis.sources.dimensions{nsource} element contains an integer characterizing the 
%number of components (time-series) characterizing the source nsource. Valid dimension numbers 
%range between 1 and the number of dimensions extracted for the original source (specified in 
%BATCH.Setup.rois.dimensions; typically just one characterizing the average BOLD signal within this ROI) 
params.conn.Analysis.sources.dimensions = {1,1};
%This field contains a cell array ranging over each source ROI, where the 
%BATCH.Analysis.sources.deriv{nsource} element contains an integer describing the order of 
%additional time-derivatives to be used as additional source time-series for source nsource 
%(typically 0 to indicate the raw BOLD signal, and no time-derivatives used)
params.conn.Analysis.sources.deriv = {0,0}; 
params.conn.Analysis.done = 1; % 1/0 variable specifying whether to run the first-level analysis steps
%'Yes'/'No' variable specifying whether to overwrite existing files if they already 
%exist when performing the first-level analysis steps. This field is only relevant 
%when BATCH.Analysis.done is set to 1. The default value is 'Yes'. Set overwrite to 
%'No' if you are adding new subjects (or ROIs) to an already processed experiment, 
%and want to estimate the first- level connectivity measures only for the new subjects (or ROIs).
params.conn.Analysis.overwrite = 'Yes';

%% analysis params
%specify song id to condition map 
params.stims.song2cond.SS = [1:6]; %self selected condition
params.stims.song2cond.FW = [7:16]; %familiar western condition
params.stims.song2cond.BP = [17:24];

%each glms(*) makes a model with the specified conditions and contrasts
params.adjust4missingContrasts = true; %if true, adjust for contrasts when possible in spm_contrasts_first_level.m

%ROIs %used in spm_grabBetas_first_level.m
params.ROIs.roi_files = dir('/scratch/xi.wu/Wang_etal_ROIs/SMPC24/nii/*.nii');
params.ROIs.rois2include = {'1bilateral_STG.nii','2bilateral_HPC.nii','3bilateral_ACC.nii','4bilateral_SMA.nii','5bilateral_MFG.nii','6mPFC.nii'};

%%%
params.glms(1).name = 'Liking';
params.glms(1).conds = {'hate_listening','neutral_listening','like_listening','love_listening','hate_liking_rating','neutral_liking_rating','like_liking_rating','love_liking_rating','hate_fam_rating','neutral_fam_rating','like_fam_rating','love_fam_rating'}; %makes one model
params.glms(1).conds2resp = {1,2,3,4}; %make sure this is same order as conds
%add F-contrasts (order needs to match between Fconts, Fweights, and Fsessrep)
params.glms(1).Fconts = {};
params.glms(1).Fweights = {}; %order of weights needs to match glms(1).conds
params.glms(1).Fsessrep = {};
%add t-contrasts (order needs to match between Fconts, Fweights, and Fsessrep)
params.glms(1).Tconts = {'T_linear_liking_listening',  'T_linear_liking_liking_rating', 'T_linear_liking_fam_rating', 'T_linear_liking_all_rating', 'T_hate_listening_>_hate_rating','T_neutral_listening_>_neutral_rating','T_like_listening_>_like_rating','T_love_listening_>_love_rating','T_Lall_listening_>_all_rating'};
params.glms(1).Tweights = {[-3 -1 1 3 0 0 0 0 0 0 0 0],[0 0 0 0 -3 -1 1 3 0 0 0 0],    [0 0 0 0 0 0 0 0 -3 -1 1 3],  [0 0 0 0 -3 -1 1 3 -3 -1 1 3], [2 0 0 0 -1 0 0 0 -1 0 0 0],     [0 2 0 0 0 -1 0 0 0 -1 0 0],           [0 0 2 0 0 0 -1 0 0 0 -1 0],     [0 0 0 2 0 0 0 -1 0 0 0 -1],    [2 2 2 2 -1 -1 -1 -1 -1 -1 -1 -1]}; %order of weights needs to match glms(1).conds
params.glms(1).Tsessrep = {'none','none','none','none','none','none','none','none','none'};

%%% 
params.glms(2).name = 'Familiarity';
params.glms(2).conds = {'vunf_listening','unf_listening','fam_listening','vfam_listening','vunf_liking_rating','unf_liking_rating','fam_liking_rating','vfam_liking_rating','vunf_fam_rating','unf_fam_rating','fam_fam_rating','vfam_fam_rating'};
params.glms(2).conds2resp = {1,2,3,4}; %make sure this is same order as conds
%add F-contrasts (order needs to match between Fconts, Fweights, and Fsessrep)
params.glms(2).Fconts = {};
params.glms(2).Fweights = {}; %order of weights needs to match glms(1).conds
params.glms(2).Fsessrep = {};
%add t-contrasts (order needs to match between Fconts, Fweights, and Fsessrep)
params.glms(2).Tconts = {'T_linear_fam_listening',     'T_linear_fam_liking_rating','T_linear_fam_fam_rating',  'T_linear_fam_all_rating',    'T_vunf_listening_>_vunf_rating','T_unf_listening_>_unf_rating','T_fam_listening_>_fam_rating','T_vfam_listening_>_vfam_rating','T_Fall_listening_>_all_rating'};
params.glms(2).Tweights = {[-3 -1 1 3 0 0 0 0 0 0 0 0],[0 0 0 0 -3 -1 1 3 0 0 0 0], [0 0 0 0 0 0 0 0 -3 -1 1 3],[0 0 0 0 -3 -1 1 3 -3 -1 1 3],[2 0 0 0 -1 0 0 0 -1 0 0 0],     [0 2 0 0 0 -1 0 0 0 -1 0 0],    [0 0 2 0 0 0 -1 0 0 0 -1 0],  [0 0 0 2 0 0 0 -1 0 0 0 -1],     [2 2 2 2 -1 -1 -1 -1 -1 -1 -1 -1]}; %order of weights needs to match glms(1).conds
params.glms(2).Tsessrep = {'none','none','none','none','none','none','none','none','none'};

%%%
params.glms(3).name = 'SongType';
params.glms(3).conds = {'ss_listening','fw_listening','bp_listening','ss_liking_rating','fw_liking_rating','bp_liking_rating','ss_fam_rating','fw_fam_rating','bp_fam_rating'};%{'self-select','familiar-western','updated-bp','selfSelect_q&a','familiarWestern_q&a','updated-bp_q&a'};
%add F-contrasts (order needs to match between Fconts, Fweights, and Fsessrep)
params.glms(3).Fconts = {};
%{'F_selfSelect_V_familiarWestern','F_selfSelect_V_BP','F_familiarWestern_V_BP','F_allFamiliar_V_allNew','F_selfSelect_V_otherSelected','F_Smusic_V_q&a','F_selfSelect_V_selfSelect_q&a','F_familiarWestern_V_familiarWestern_q&a','F_BP_V_BP_q&a','F_selfSelect_q&a_V_BP_q&a','F_selfSelect_q&a_V_familiarWestern_q&a','F_familiarWestern_q&a_V_BP_q&a'};
params.glms(3).Fweights = {};%{[1 -1 0 0 0 0],                [1 0 -1 0 0 0],     [0 1 -1 0 0 0],         [1 1 -2 0 0 0],          [2 -1 -1 0 0 0],               [1 1 1 -1 -1 -1], [1 0 0 -1 0 0],                   [0 1 0 0 -1 0],                           [0 0 1 0 0 -1],  [0 0 0 1 0 -1],             [0 0 0 1 -1 0],                          [0 0 0 0 1 -1]}; %order of weights needs to match glms(1).conds
params.glms(3).Fsessrep = {};
%add t-contrasts (order needs to match between Fconts, Fweights, and Fsessrep)
params.glms(3).Tconts = {'T_ss_listening_>_os_listening','T_ss_listening_>_all_rating','T_fw_listening_>_all_rating','T_bp_listening_>_all_rating','T_ss_listening_>_ss_rating','T_fw_listening_>_fw_rating','T_bp_listening_>_bp_rating','T_ss_listening_>_fw_listening','T_ss_listening_>_bp_listening','T_fw_listening_>_bp_listening','T_ss_rating_>_bp_rating','T_ss_rating_>_fw_rating','T_fw_rating_>_bp_rating','T_Sall_listening_>_all_rating'};params.glms(3).Tweights = {[2 0 0 -1 0 0 -1 0 0],     [0 2 0 0 -1 0 0 -1 0],       [0 0 2 0 0 -1 0 0 -1],      [1 -1 0 0 0 0 0 0 0],           [1 0 -1 0 0 0 0 0 0],           [0 1 -1 0 0 0 0 0 0],           [0 0 0 1 0 -1 1 0 -1],    [0 0 0 1 -1 0 1 -1 0],    [0 0 0 0 1 -1 0 1 -1],    [2 2 2 -1 -1 -1 -1 -1 -1]}; %order of weights needs to match glms(1).conds
params.glms(3).Tweights = {[2 -1 -1 0 0 0 0 0 0],        [6 0 0 -1 -1 -1 -1 -1 -1],[0 6 0 -1 -1 -1 -1 -1 -1],[0 0 6 -1 -1 -1 -1 -1 -1],  [2 0 0 -1 0 0 -1 0 0],     [0 2 0 0 -1 0 0 -1 0],       [0 0 2 0 0 -1 0 0 -1],      [1 -1 0 0 0 0 0 0 0],           [1 0 -1 0 0 0 0 0 0],           [0 1 -1 0 0 0 0 0 0],           [0 0 0 1 0 -1 1 0 -1],    [0 0 0 1 -1 0 1 -1 0],    [0 0 0 0 1 -1 0 1 -1],    [2 2 2 -1 -1 -1 -1 -1 -1]}; %order of weights needs to match glms(1).conds
params.glms(3).Tsessrep = {'none','none','none','none','none','none','none','none','none','none','none','none','none','none'};

params.allContrasts = {}; %list of all contrasts below (needs to match! used to make output table)
for iglm=1:size(params.glms,2)
    params.allContrasts = {params.allContrasts{:},params.glms(iglm).Fconts{:},params.glms(iglm).Tconts{:}};
end



%% Subjects (ADD NEW SUBS HERE)

%{
~~~TEMPLATE TO PASTE~~~
isub=isub+1;
subjectInfo(isub).batch = ''; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = '';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'MCI'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre';'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'';''}; %in params.paths.data
subjectInfo(isub).logfiles = {'_mus-bid-biddata.txt';'_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/wausub-';'func/swausub-'}; %in params.paths.data

%}

isub=0; %init subject counter

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = 'FKEE';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'001'}; %in params.paths.data
subjectInfo(isub).logfiles = {'190716FKEE1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-01_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = 'JSCH';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'003'}; %in params.paths.data
subjectInfo(isub).logfiles = {'190717JSCH1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-03_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = 'GHER';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'006'}; %in params.paths.data
subjectInfo(isub).logfiles = {'190904GHER1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-06_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = 'JPRE';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'010'}; %in params.paths.data
subjectInfo(isub).logfiles = {'191113JPRE1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-010_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = 'MBLO';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'012'}; %in params.paths.data
subjectInfo(isub).logfiles = {'191204MBLO1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-012_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %(gmbi/mbi/healthy/) should match log folder names; also name of second level anals folders
subjectInfo(isub).subjectID = 'JKAY';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %%MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'021'}; %in params.paths.data
subjectInfo(isub).logfiles = {'200204JKAY1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-021_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'SCAE';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'028'}; %in params.paths.data
subjectInfo(isub).logfiles = {'201210SCAE1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-028_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'LJAC';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'030'}; %in params.paths.data
subjectInfo(isub).logfiles = {'210202LJAC1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-030_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'ASHA';
subjectInfo(isub).intervention = 'MBI';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'032'}; %in params.paths.data
subjectInfo(isub).logfiles = {'210223ASHA1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-032_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'mbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'BFOR';
subjectInfo(isub).intervention = 'MBI'; % lights optional
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'040'}; %in params.paths.data
subjectInfo(isub).logfiles = {'210719BFOR1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-040_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'DMUN';
subjectInfo(isub).intervention = 'gMBI'; 
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'087'; '090'}; %in params.paths.data
subjectInfo(isub).logfiles = {'220811DMUN1_mus-bid-biddata.txt'; '220811DMUN3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-087_task-musbid_bold.nii'; 'func/swausub-090_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'GHER';
subjectInfo(isub).intervention = 'gMBI'; 
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'088'; '093'}; %in params.paths.data
subjectInfo(isub).logfiles = {'220914GHER1_mus-bid-biddata.txt'; '220914GHER3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-088_task-musbid_bold.nii'; 'func/swausub-093_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'JDEN';
subjectInfo(isub).intervention = 'gMBI'; 
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'089'; '092'}; %in params.paths.data
subjectInfo(isub).logfiles = {'220919JDEN1_mus-bid-biddata.txt'; '220913JDEN3_mus-bid-biddata.txt'}; % happened before GHER3 tho
subjectInfo(isub).niifiles = {'func/swausub-089_task-musbid_run-02_bold.nii'; 'func/swausub-092_task-musbid_run-03_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'EKAS';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'091'; '094'}; %in params.paths.data
subjectInfo(isub).logfiles = {'221103EKAS1_mus-bid-biddata.txt'; '221103EKAS3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-091_task-musbid_bold.nii'; 'func/swausub-094_task-musbid_bold.nii'}; %in params.paths.data
%{
isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'PMUR';
subjectInfo(isub).intervention = 'gMBI'; % did gMBI, but dropped out before posttest
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'102'}; %in params.paths.data
subjectInfo(isub).logfiles = {'230612PMUR1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-102_task-musbid_bold.nii'}; %in params.paths.data
%}
isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'FGOD';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'103'; '109'}; %in params.paths.data
subjectInfo(isub).logfiles = {'230803FGOD1_mus-bid-biddata.txt'; '230803FGOD3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-103_task-musbid_run-02_bold.nii'; 'func/swausub-109_task-musbid_run-05_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'RFIS';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'104'; '111'}; %in params.paths.data
subjectInfo(isub).logfiles = {'230911RFIS1_mus-bid-biddata.txt'; '230911RFIS3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-104_task-musbid_bold.nii'; 'func/swausub-111_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1; 
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'JPRI';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'105'; '112'}; %in params.paths.data
subjectInfo(isub).logfiles = {'230913JPRI1_mus-bid-biddata.txt'; '230913JPRI3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-105_task-musbid_run-03_bold.nii'; 'func/swausub-112_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'PMUL';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'106'; '113'}; %in params.paths.data
subjectInfo(isub).logfiles = {'230925PMUL1_mus-bid-biddata.txt'; '230925PMUL3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-106_task-musbid_bold.nii'; 'func/swausub-113_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'LJOS';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'107';'114'}; %in params.paths.data
subjectInfo(isub).logfiles = {'231002LJOS1_mus-bid-biddata.txt'; '231002LJOS3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-107_task-musbid_bold.nii'; 'func/swausub-114_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'JPEN';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'MCI'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'108'; '115'}; %in params.paths.data
subjectInfo(isub).logfiles = {'231005JPEN1_mus-bid-biddata.txt'; '231005JPEN3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-108_task-musbid_bold.nii'; 'func/swausub-115_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'gmbi'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'JCHA';
subjectInfo(isub).intervention = 'gMBI';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'; 'post'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'110'; '118'}; %in params.paths.data
subjectInfo(isub).logfiles = {'231013JCHA1_mus-bid-biddata.txt'; '231013JCHA3_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-110_task-musbid_bold.nii'; 'func/swausub-118_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'healthy'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'RFOL';
subjectInfo(isub).intervention = 'NA';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'122'}; %in params.paths.data
subjectInfo(isub).logfiles = {'240327RFOL1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-122_task-musbid_bold.nii'}; %in params.paths.data

isub=isub+1;
subjectInfo(isub).batch = 'healthy'; %should match log folder names (gmbi/ or healthy/...)
subjectInfo(isub).subjectID = 'TWHI';
subjectInfo(isub).intervention = 'NA';
subjectInfo(isub).MCI = 'healthy'; %MCI status: MCI, healthy, NA
subjectInfo(isub).sessNames = {'pre'}; %should match log folder names
subjectInfo(isub).subjectNUM = {'124'}; %in params.paths.data
subjectInfo(isub).logfiles = {'240405TWHI1_mus-bid-biddata.txt'}; %in params.paths.logs
subjectInfo(isub).niifiles = {'func/swausub-124_task-musbid_bold.nii'}; %in params.paths.data

fprintf(sprintf('Hello %s! Good luck with your work today!\n',user))
