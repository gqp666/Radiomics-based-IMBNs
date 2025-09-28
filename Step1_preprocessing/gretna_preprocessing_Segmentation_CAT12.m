function gretna_preprocessing_Segmentation_CAT12(Data_path, File_filter, Para)

%==========================================================================
% This function is used to perform tissue segmentation of images (typically
% structural MRI images) with the CAT12 toolbox.
%
%
% Syntax: function gretna_preprocessing_Segmentation_CAT12(Data_path, File_filter, Para)
%
% Inputs:
%         Data_path:
%                   The directory & filename of a .txt file that contains
%                   the directory of those files to be processed (can be
%                   obtained by gretna_gen_data_path.m).
%       File_filter:
%                   The prefix of those files to be processed.
%   Para (optional):
%           Para.Segment.Nthreads:
%                   The number of threads for parallel calculation.
%           Para.Segment.Regularisation:
%                   'European':   Affine regularisation using European
%                                 template.
%                   'East Asian': Affine regularisation using East Asian
%                                 template.
%           Para.Segment.Registration:
%                   'Shooting': Spatial shooting registration
%                               (Ashburner, 2008).
%                   'Dartel':   Spatial dartel registration
%                               (Ashburner, 2011).
%           Para.Segment.Type:
%                   'yes': VBM with surface and thickness estimation.
%                   'no':  VBM without surface and thickness estimation.
%           Para.Segment.TPM_path:
%                   The TPM used to initial spm segment. It is ok to use
%                   SPM default for very old/young brains. Nevertheless,
%                   for children data, it is recommended to use
%                   customized TPMs created with the Template-O-Matic
%                   toolbox.
%
% Ningkai WAng,IBRR, SCNU, Guangzhou, 2020/09/14, ningkai.wang.1993@gmail.com
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/10/24, jinhui.wang.1982@gmail.com
%==========================================================================

%% spm_input
if nargin == 2
    Para.Segment.Nthreads       = spm_input('Number of Threads',                 1,'e',[],1);
    Para.Segment.Regularisation = spm_input('Affine Regularisation',             2, 'European|East Asian');
    Para.Segment.Registration   = spm_input('Spatial Registration',              3, 'Shooting|Dartel');
    Para.Segment.Type           = spm_input('Whether Surface and CT Estimation', 4,'yes|no');
    Para.Segment.TPM            = spm_input('TPM',                               5, 'SPM Default|Study Specific');
    
    if strcmp(Para.Segment.TPM, 'SPM Default')
        if exist('spm','file') == 2
            spm_dir = which('spm');
            [pathstr, ~, ~] = fileparts(spm_dir);
            Para.Segment.TPM_path = fullfile(pathstr, 'tpm', 'TPM.nii');
        else
            error('Cannot find SPM toolbox in Matlab search path of your computer!!')
        end
    else
        Para.Segment.TPM_path = spm_input('Enter TPM Path', 6, 's');
    end
end
close

%% update batch parameters
load gretna_Segmentation_CAT12.mat
batch_segment = matlabbatch;

batch_segment{1}.spm.tools.cat.estwrite.nproc = Para.Segment.Nthreads;

if exist('cat12','file') == 2
    cat_dir = which('cat12');
    [pathstr, ~, ~] = fileparts(cat_dir);
    
    % affine regularisation
    switch lower(Para.Segment.Regularisation)
    case 'european'
        batch_segment{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
    case 'east asian'
        batch_segment{1}.spm.tools.cat.estwrite.opts.affreg = 'eastern';
    end

    % spatial registration
    switch lower(Para.Segment.Registration)
    case 'shooting'
        batch_segment{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm{1} = ...
            fullfile(pathstr, 'templates_volumes','Template_0_IXI555_MNI152_GS.nii');
    case 'dartel'
        batch_segment{1}.spm.tools.cat.estwrite.extopts.registration.dartel.darteltpm{1} = ...
            fullfile(pathstr, 'templates_volumes', 'Template_1_IXI555_MNI152.nii');
    end
else
    error('Cannot find CAT toolbox in Matlab search path of your computer!!')
end

% surface estimation
if strcmpi(Para.Segment.Type, 'yes')
    batch_segment{1}.spm.tools.cat.estwrite.output.surface = 1;
else
    batch_segment{1}.spm.tools.cat.estwrite.output.surface = 0;
end

% TPM
batch_segment{1}.spm.tools.cat.estwrite.opts.tpm{1} = Para.Segment.TPM_path;

%% Update batch data
fid      = fopen(Data_path);
Dir_data = textscan(fid, '%s');
fclose(fid);

Num_subs = size(Dir_data{1},1);
Sour_all = cell(Num_subs,1);

for isub = 1:Num_subs
    cd([Dir_data{1}{isub}])
    Sour_ind = spm_select('ExtList', pwd, ['^' File_filter '.*' filesep '.nii$'],inf);

    if isempty(Sour_ind)
        Sour_ind = spm_select('ExtList', pwd, ['^' File_filter '.*' filesep '.img$'],inf);
    end

    Sour_all{isub,1} = [Dir_data{1}{isub} filesep Sour_ind];
end

batch_segment{1}.spm.tools.cat.estwrite.data = Sour_all;

% Run batch
spm_jobman('run',batch_segment);

return