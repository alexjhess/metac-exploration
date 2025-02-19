%% [METAC] exploration data set - SBF Analysis


%% remove default toolboxes that could be conflicting

% List all paths in search path and break them into a cell array
dirs = regexp(path,['[^;]*'],'match');

% Return all search path directories containing this string
strsToFind = {'tapas', 'hgf-toolbox', 'comb_obs_models', 'VBA-toolbox'}; 

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
cd('VBA-toolbox')
VBA_setup();
cd ..

%% load task data
% dat = load_discovery_set();
load(fullfile('data', 'discovery_set_tmp.mat'));

%% fit behavioural models

% calc ABF from GBF

% dat = metac_fit_behaviour(dat);
load(fullfile('data', 'discovery_set_fits_tmp.mat'));






%% load ppids

id = readtable(fullfile('data', 'metac_ppids_exploration_set.txt'));
ppidstr = cell(size(id,1),1);
for i = 1: size(id,1)
    ppidstr{i} = ['METAC_', num2str(id{i,1})];
end


%% load questionnaire data

% datadir = fullfile('P:METAC_Iglesias', 'Data');
datadir = fullfile('data');
f = dir(fullfile(datadir, 'METAC1*.csv'));

%% load ppids of discovery set
df = readtable(fullfile(f.folder, f.name));
idx = find(contains(df.ppid, ppidstr));
df_discovery = df(idx,:);

%% extract FAS scores and features
y_fas = df_discovery.fas_exp_total_score(find(~isnan(df_discovery.fas_exp_total_score)));
X = [df_discovery.demo_age(find(~isnan(df_discovery.demo_age))),...
    df_discovery.demo_gender(find(~isnan(df_discovery.demo_gender)))...
    ];

%% extract other features
% avg overall control, tolerance, ...



%% create table for winning mod 1
m = dat.ffx.idx;
pars = dat.param_mat(dat.mod(m).obs_idx(1:end-1),:)';
tab = table(y_fas, X(:,1), X(:,2), ...
    pars(:,1), pars(:,2), pars(:,3), pars(:,4));
tab.Properties.VariableNames = {'FAS', 'age', 'gender',...
    'gamma', 'shift', 'scale', 'w'};

%% normalize values ?
tab_norm = normalize(tab);


%% write csv file with data for Bayesian ANCOVA
save_path = fullfile('data', ['tmp_data.csv']);
writetable(tab, save_path);

save_path2 = fullfile('data', ['tmp_norm_data.csv']);
writetable(tab_norm, save_path2);

