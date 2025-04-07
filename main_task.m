%% [METAC] exploration data set - task data


%% remove default toolboxes that could be conflicting

% List all paths in search path and break them into a cell array
dirs = regexp(path,['[^;]*'],'match');

% Return all search path directories containing this string
strsToFind = {'tapas', 'hgf-toolbox', 'comb_obs_models', 'mc_obs', 'VBA-toolbox'}; 

% loop over strings & entries
for i = 1:numel(strsToFind)
    strToFind = char(strsToFind(i));
    % Index to cell entries containing the desired string
    whichCellEntry = find(cellfun(@(dirs) contains(dirs, strToFind), dirs) == 1);
    % remove from path
    for j = 1:length(whichCellEntry)
        rmpath(char(dirs(whichCellEntry(j))))
    end
end

%% add hgf toolbox
addpath(genpath('hgf-toolbox'));
addpath('comb_obs_models');
addpath('mc_obs');
cd('VBA-toolbox')
VBA_setup();
cd ..

%% load task data
% dat = load_task_discovery_set();
load(fullfile('data', 'discovery_set_tmp.mat'));

%% preprocess task data
dat.pdat = metac_preproc(dat);

%% visualise raw data
metac_plot_raw_data(dat);


%% MC: autoreg obs models
dat = metac_mc_autoreg_obs(dat);



%% INT + MC: 3-level binary HGF MAB + autoreg obs model
dat = metac_int_hgf_binary_mc_autoreg_obs_mab(dat);


%% MC: cont eHGF on raw PE (scale of PE ???)
dat = metac_raw_PE_mc_ehgf(dat);

%% INT: 3-level binary HGF MAB
dat = metac_int_hgf_binary_mab(dat);

%% INT: 3-level binary eHGF MAB (ERROR..)
dat = metac_int_ehgf_binary_mab(dat);

%% INT: classical 3-level binary eHGF (NOT APPROPRIATE!!!)
dat = metac_int_ehgf_binary(dat);











