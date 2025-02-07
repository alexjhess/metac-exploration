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

%% dock all figures
set(0,'DefaultFigureWindowStyle','docked')


%% load data (exploration set)

basedir = fullfile('R:METAC', 'behavior', 'raw');
fname = fullfile('behavior', 'task', '*.csv');

id = readtable(fullfile('data', 'metac_ppids_exploration_set.txt'));

u_bin = NaN(80,size(id,1));
y_pred = NaN(80,size(id,1));
u_ymc = NaN(80,size(id,1));

for n = 1:size(id,1)
    ppid = ['TNU_METAC_', num2str(id{n,1})];
    disp(ppid)
    d = dir(fullfile(basedir, ppid, fname));
    for i = 1:size(d,1)
        if contains(d(i).name, 'experiment')
            % specs stored in rows 1-2
            tmp = readtable(fullfile(d(i).folder, d(i).name), 'NumHeaderLines', 2);
        end
    end
    u_bin(:,n) = tmp.jSuccess(1:80);
    y_pred(:,n) = tmp.prediction(1:80);
    try
        y_mc(:,n) = tmp.control(1:80);
    catch
        y_mc(:,n) = tmp.expl_control(1:80);
    end
end

u_pe = u_bin - y_pred;

%% create data struct
dat.u_bin = u_bin;
dat.y_pred = y_pred;
dat.y_mc = y_mc;
dat.u_pe = u_pe;

%% save struct
save('data\pilots_tmp.mat', 'dat', '-mat');

%% load 
load('data\pilots_tmp.mat')


%% plot success, avg success rate and avg prediction

figure;
for n = 1:size(dat.u_bin,2)
    subplot(ceil(size(dat.u_bin,2)/2), 2, n)
    plot(dat.u_bin(:,n), '.')
    title(['sub ', num2str(n)])
    hold on;
    plot(mean(dat.y_pred,2))
    plot(mean(dat.u_bin,2))
    hold off;
end
legend('u_bin', 'avg y_{int}', 'avg u_{bin}')


%% HGF MAB example (GitHub Issues)
% https://github.com/translationalneuromodeling/tapas/issues/243

% example data
load('data\hgf_binary_mab_example\hgf_binary_mab.mat')

% example fit
fit = tapas_fitModel(y,...
u,...
'tapas_hgf_binary_mab_config',...
'tapas_bayes_optimal_binary_config',...
'tapas_quasinewton_optim_config');
tapas_hgf_binary_mab_plotTraj(fit)











% %% fit interoceptive HGF_mab (Bayes Optimal)
% 
% hgf_mab_config = tapas_hgf_binary_mab_config;
% hgf_mab_config.n_bandits = 2;
% 
% for n = 1%:12
%     bo_est = tapas_fitModel(u_mab,...
%         dat.u_bin(:,n),...
%         hgf_mab_config,...
%         tapas_bayes_optimal_binary_config,...
%         tapas_quasinewton_optim_config ...
%         );
%     tapas_hgf_binary_mab_plotTraj(bo_est)
% end
% 
% 
% %% _______________
% % sim interoceptive HGF_mab
% 
% for n = 1%:12
%     sim = tapas_simModel([u_bin(:,n) u_mab],...
%         'tapas_hgf_binary_mab',...
%         [NaN 0 1 NaN 0.1000 1 NaN 0 0 1 1 NaN -8.1733 -6.0410],... default [NaN 0 1 NaN 0.1 1 NaN 0 0 1 1 NaN -2 -6],...
%         'tapas_beta_obs',...
%         log(128),...
%         12345 ...
%         );
%     tapas_hgf_binary_mab_plotTraj(sim)
%     hold on;
%     plot(sim.y)
%     hold off;
% end



%% CREATE SEPARATE INPUTS FOR EACH CONDITION

% pre-allocation
u_sep = NaN(size(exp_seq,1), 2);

n=1;
% different input seq for wind / no wind
for k = 1:size(exp_seq,1)
    if exp_seq.wind(k)==0
        u_sep(k,1) = u_bin(k,n);
    elseif exp_seq.wind(k)==1 || exp_seq.wind(k)==-1
        u_sep(k,2) = u_bin(k,n);
    end
end

figure
plot(u_sep, '.')
ylim([-0.1 1.1])
legend('no wind', 'wind')

%% fit individual models for each condition

n=1;

% no wind
est1 = tapas_fitModel(y_pred(:,n), ...
    u_sep(:,1), ...
    tapas_rw_binary_config, ...
    tapas_beta_obs_config, ...
    tapas_quasinewton_optim_config ...
    );

tapas_rw_binary_plotTraj(est1)

% wind
est2 = tapas_fitModel(y_pred(:,n), ...
    u_sep(:,2), ...
    tapas_rw_binary_config, ...
    tapas_beta_obs_config, ...
    tapas_quasinewton_optim_config ...
    );

tapas_rw_binary_plotTraj(est2)
hold on;
plot(est1.y, '.')
plot(est1.traj.v)

%%
n=1;
figure
plot(est1.traj.da+est2.traj.da)
hold on
plot(u_pe(:,n))
legend('RW inferred PE', 'raw PE')

%%
figure
bar([est1.optim.negLl est2.optim.negLl est1.optim.negLl+est2.optim.negLl])

%%
figure
bar(est1.optim.trialLogLlsplit)
hold on;
bar(est2.optim.trialLogLlsplit)
hold off;








