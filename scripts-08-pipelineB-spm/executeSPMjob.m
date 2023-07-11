function [] = executeSPMjob(taskName,spaces,ssKernel,hpfValue,spmFolder,fmriPrepFolder,bidsFolder,subjectID)
%EXECUTESPMJOB Perform batch in SPM12
%   Spatial smoothing, GLM
%
% Inputs:
%   taskName
%   spaces - a list of spaces to be analysed
%   ssKernel - spatial smoothing kernel (in mm)
%   hpfValue - cut-off value for high-pass filtering (in seconds)
%   spmFolder
%   fmriPrepFolder
%   bidsFolder
%   subjectID
%
% Author: Alexandre Sayal
% CIBIT, University of Coimbra
% February 2022
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

fprintf('[executeSPMjob] Started for run %s.\n',taskName)

%% Convert .tsv to SPM multiple conditions file
tsvFile = fullfile(bidsFolder,subjectID,'ses-01','func',[subjectID '_ses-01_' taskName '_events.tsv']);
tsvTable = importTSVProtocol(tsvFile);

names = cellstr(unique(tsvTable.trial_type));

% Remove baseline condition from matrix - in this case, it is condition
% called 'Static'
%names(strcmp(names,'Static')) = [];

onsets = cell(size(names));
durations = cell(size(names));

for cc = 1:length(names)
    
    onsets{cc} = tsvTable.onset(tsvTable.trial_type==names{cc})';
    durations{cc} = tsvTable.duration(tsvTable.trial_type==names{cc})';
    
end

save(fullfile(spmFolder,['protocol_' taskName '.mat']),...
    'names','onsets','durations')

%% Fetch run info

% Get number of volumes, slices, and TR
funcImage = fullfile(fmriPrepFolder,subjectID,'ses-01','func',...
    [subjectID '_ses-01_' taskName '_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz']);

[~,tr]=system(sprintf("echo $(fslinfo %s | grep -m 1 pixdim4 | awk '{print $2}')",funcImage));
tr = str2double(tr);

% taskProperName = strsplit(taskName,'_');
% taskProperName = taskProperName{1};

spm('defaults', 'FMRI');

%% Iterate
for sp = 1:length(spaces)
    
    % Start important folders
    mkdir(fullfile(spmFolder,['3d_' taskName '_' spaces{sp}]))
    mkdir(fullfile(spmFolder,['model_' taskName '_' spaces{sp}]))
    
    clear matlabbatch
    
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files = {funcImage};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir = {spmFolder};
    matlabbatch{1}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep = true;
    
    matlabbatch{2}.spm.util.split.vol(1) = cfg_dep('Gunzip Files: Gunzipped Files', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{':'}));
    matlabbatch{2}.spm.util.split.outdir = {fullfile(spmFolder,['3d_' taskName '_' spaces{sp}])};
    
    matlabbatch{3}.spm.spatial.smooth.data(1) = cfg_dep('4D to 3D File Conversion: Series of 3D Volumes', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','splitfiles'));
    matlabbatch{3}.spm.spatial.smooth.fwhm = [ssKernel ssKernel ssKernel];
    matlabbatch{3}.spm.spatial.smooth.dtype = 0;
    matlabbatch{3}.spm.spatial.smooth.im = 0;
    matlabbatch{3}.spm.spatial.smooth.prefix = 'ss';
    
    matlabbatch{4}.spm.stats.fmri_spec.dir = {fullfile(spmFolder,['model_' taskName '_' spaces{sp}])};
    matlabbatch{4}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{4}.spm.stats.fmri_spec.timing.RT = tr;
    matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t = 16;
    matlabbatch{4}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{4}.spm.stats.fmri_spec.sess.scans(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{4}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi = {fullfile(spmFolder,['protocol_' taskName '.mat'])};
    matlabbatch{4}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
    matlabbatch{4}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(spmFolder,[subjectID '_ses-01_' taskName '_desc-PNMmodel.mat'])};
    matlabbatch{4}.spm.stats.fmri_spec.sess.hpf = hpfValue;
    matlabbatch{4}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
    matlabbatch{4}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{4}.spm.stats.fmri_spec.volt = 1;
    matlabbatch{4}.spm.stats.fmri_spec.global = 'None';
    matlabbatch{4}.spm.stats.fmri_spec.mthresh = 0.8;
    matlabbatch{4}.spm.stats.fmri_spec.mask = {''};
    matlabbatch{4}.spm.stats.fmri_spec.cvi = 'FAST';
    
    matlabbatch{5}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{5}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{5}.spm.stats.fmri_est.method.Classical = 1;
    
    % This is of course task-specific
    matlabbatch{6}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{6}.spm.stats.con.consess{1}.tcon.name = 'Coherent > Static';
    matlabbatch{6}.spm.stats.con.consess{1}.tcon.weights = [0 0 0 1 0 0 0 0 0 0 0 0 -1];
    matlabbatch{6}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{6}.spm.stats.con.consess{2}.tcon.name = 'Incoherent > Static';
    matlabbatch{6}.spm.stats.con.consess{2}.tcon.weights = [0 0 0 0 0 0 0 0 1 0 0 0 -1];
    matlabbatch{6}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    matlabbatch{6}.spm.stats.con.consess{3}.tcon.name = 'NonAdapt > Static';
    matlabbatch{6}.spm.stats.con.consess{3}.tcon.weights = [0 0 0 0 0 0 0 0 0 0 1 0 -1];
    matlabbatch{6}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
    matlabbatch{6}.spm.stats.con.consess{4}.tcon.name = 'Coherent+Incoherent+NonAdapt > Static';
    matlabbatch{6}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1 0 0 0 0 1 0 1 0 -3];
    matlabbatch{6}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
    matlabbatch{6}.spm.stats.con.consess{5}.fcon.name = 'Effects of interest';
    matlabbatch{6}.spm.stats.con.consess{5}.fcon.weights = eye(13);
    matlabbatch{6}.spm.stats.con.consess{5}.fcon.sessrep = 'none';
    matlabbatch{6}.spm.stats.con.consess{6}.tcon.name = 'Coh_aCoh > Coherent';
    matlabbatch{6}.spm.stats.con.consess{6}.tcon.weights = [1 0 0 -1 0 0 0 0 0 0 0 0 0];
    matlabbatch{6}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
    matlabbatch{6}.spm.stats.con.consess{7}.tcon.name = 'InCoh_aCoh > Coherent';
    matlabbatch{6}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 -1 0 1 0 0 0 0 0 0 0];
    matlabbatch{6}.spm.stats.con.consess{7}.tcon.sessrep = 'none';  
    matlabbatch{6}.spm.stats.con.consess{8}.tcon.name = 'Coh_aInCoh > Incoherent';
    matlabbatch{6}.spm.stats.con.consess{8}.tcon.weights = [0 1 0 0 0 0 0 0 -1 0 0 0 0];
    matlabbatch{6}.spm.stats.con.consess{8}.tcon.sessrep = 'none'; 
    matlabbatch{6}.spm.stats.con.consess{9}.tcon.name = 'InCoh_aInCoh > Incoherent';
    matlabbatch{6}.spm.stats.con.consess{9}.tcon.weights = [0 0 0 0 0 0 1 0 -1 0 0 0 0];
    matlabbatch{6}.spm.stats.con.consess{9}.tcon.sessrep = 'none';     
    matlabbatch{6}.spm.stats.con.consess{10}.tcon.name = 'Coh_aNA > NonAdapt';
    matlabbatch{6}.spm.stats.con.consess{10}.tcon.weights = [0 0 1 0 0 0 0 0 0 0 -1 0 0];
    matlabbatch{6}.spm.stats.con.consess{10}.tcon.sessrep = 'none';      
    matlabbatch{6}.spm.stats.con.consess{11}.tcon.name = 'InCoh_aNA > NonAdapt';
    matlabbatch{6}.spm.stats.con.consess{11}.tcon.weights = [0 0 0 0 0 0 0 1 0 0 -1 0 0];
    matlabbatch{6}.spm.stats.con.consess{11}.tcon.sessrep = 'none';       
    matlabbatch{6}.spm.stats.con.delete = 0;
    
    matlabbatch{7}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{6}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{7}.spm.stats.results.conspec.titlestr = '';
    matlabbatch{7}.spm.stats.results.conspec.contrasts = Inf;
    matlabbatch{7}.spm.stats.results.conspec.threshdesc = 'none';
    matlabbatch{7}.spm.stats.results.conspec.thresh = 0.001;
    matlabbatch{7}.spm.stats.results.conspec.extent = 50;
    matlabbatch{7}.spm.stats.results.conspec.conjunction = 1;
    matlabbatch{7}.spm.stats.results.conspec.mask.none = 1;
    matlabbatch{7}.spm.stats.results.units = 1;
    matlabbatch{7}.spm.stats.results.export = cell(1,0);
    
    matlabbatch{8}.spm.util.cat.vols(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
    matlabbatch{8}.spm.util.cat.name = fullfile(spmFolder,[subjectID '_ses-01_' taskName '_space-' spaces{sp} '_desc-SS_bold.nii']);
    matlabbatch{8}.spm.util.cat.dtype = 4;
    matlabbatch{8}.spm.util.cat.RT = tr;
    
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.files(1) = cfg_dep('3D to 4D File Conversion: Concatenated 4D Volume', substruct('.','val', '{}',{8}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','mergedfile'));
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.outdir = {''};
    matlabbatch{9}.cfg_basicio.file_dir.file_ops.cfg_gzip_files.keep = false;
    
    matlabbatch{10}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.files = {fullfile(fmriPrepFolder,subjectID,'ses-01','anat',[subjectID '_ses-01_run-1_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz'])};
    matlabbatch{10}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.outdir = {spmFolder};
    matlabbatch{10}.cfg_basicio.file_dir.file_ops.cfg_gunzip_files.keep = true;
        
    % RUN
    spm_jobman('run', matlabbatch);
    
    % Delete temporary files
    delete(fullfile(spmFolder,[subjectID '_ses-01_' taskName '_space-' spaces{sp} '_desc-preproc_bold.nii'])) % unzipped fmriprep output file    
    delete(fullfile(spmFolder,['3d_' taskName '_' spaces{sp}],'sub-*.nii')) % 3d volumes unsmoothed
    
end

fprintf('[executeSPMjob] Finished for run %s.\n',taskName)

end
