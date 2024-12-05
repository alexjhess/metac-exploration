%% [METAC] exploration data set - task data


%% remove default toolboxes that could be conflicting

% List all paths in search path and break them into a cell array
dirs = regexp(path,['[^;]*'],'match');

% Return all search path directories containing this string
strsToFind = {'tapas'}; 

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
addpath('comb_obs_models')


%% load data

basedir = fullfile('R:METAC', 'behavior', 'raw');
fname = fullfile('behavior', 'task', '*.csv');

id = readtable(fullfile('data', 'metac_ppids_exploration_task_v3_1.txt'));

for n = 1%:size(id,1)
    ppid = ['TNU_METAC_', num2str(id{n,1})];
    d = dir(fullfile(basedir, ppid, fname));
    for i = 1:size(d,1)
        if contains(d(i).name, 'experiment')
            % specs stored in rows 1-2
            dat = readtable(fullfile(d(i).folder, d(i).name), 'NumHeaderLines', 2)
        end
    end
end


%% get responses

u_bin = dat.jSuccess;
y_pred = dat.prediction;
u_pe = u_bin - y_pred;
y_mc = dat.control;


%% plot raw data

figure;
plot(y_mc)
hold on;
plot(1-abs(u_pe))
hold off;

corr(y_mc,1-abs(u_pe))


%% build int HGF

sim = tapas_simModel(u_bin,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'tapas_beta_obs',...
    log(128),...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(sim)


%% fit in HGF

est = tapas_fitModel(y_pred,...
    u_bin,...
    tapas_ehgf_binary_config,...
    tapas_beta_obs_config,...
    tapas_quasinewton_optim_config ...
    );
tapas_ehgf_binary_plotTraj(est)


%% ____________
% build metacognitive HGF

sim = tapas_simModel((1.-abs(u_pe)),...
    'tapas_ehgf',...
    [0 1 1 1 0 0 1 -2 -2 10],...
    'tapas_beta_obs',...
    log(128),...
    12345 ...
    );
tapas_ehgf_plotTraj(sim)


%% fit metacognitive HGF

est = tapas_fitModel(10.*y_mc,...
    10.*(1.-abs(u_pe)),...
    tapas_ehgf_config,...
    tapas_gaussian_obs_config,...
    tapas_quasinewton_optim_config ...
    );
tapas_ehgf_plotTraj(est)


%% ____________
% build autoregr mc obs model

sim = tapas_simModel(u_bin,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_autoreg_pe_inobs',... 'mc_null_inobs',...
    [0.05 1 -2 -1.5 -0.04 0.005],... [0.05 -2 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(sim)
hold on
plot(sim.y(:,2))

%% null model
sim = tapas_simModel(u_bin,...
    'tapas_ehgf_binary',...
    [NaN 0 1 NaN 1 1 NaN 0 0 1 1 NaN -2 2],...
    'mc_null_inobs',...
    [0.05 0.5 0.005],...
    12345 ...
    );
tapas_ehgf_binary_plotTraj(sim)
hold on
plot(sim.y(:,2))


%% fit metacognitive HGF

est = tapas_fitModel(y_mc,...
    u_bin,...
    tapas_ehgf_binary_config,...
    mc_autoreg_pe_inobs_config,...
    tapas_quasinewton_optim_config ...
    );
tapas_ehgf_binary_plotTraj(est)





