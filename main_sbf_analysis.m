%% [METAC] exploration data set - SBF Analysis


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

%% load task data (EXPLORATION)
dataset = 1; % 1=discovery set
% dat = load_task_data(dataset);
load(fullfile('data', 'discovery_set_tmp.mat'));
ds.dat = dat;

%% preprocess task data
ds.pdat = metac_preproc(ds.dat);
    
%% create model space
[ds.ip.mod, bo] = metac_create_model_space(1); % 1=logit space

%% grid search on param values
n_pe = 1; % PE traj from pilot sub n_pe
metac_init_priors_grid_search(ds.dat, ds.pdat, ds.ip.mod, n_pe);

%% fit behavioural models & extract other task features
priors = 1; % initial priors
[ds.gest, ds.mod] = metac_fit_behaviour(ds.dat, ds.pdat, ds.ip.mod, priors);
save(fullfile('data', 'discovery_set_fits_tmp.mat'), 'ds', '-mat');

% load(fullfile('data', 'discovery_set_fits_tmp.mat'));

%% avg overall control, tolerance, ...
ds.dat.task.avg_c = mean(ds.dat.y_c,'omitnan')';
ds.dat.task.avg_tol = mean(ds.dat.y_tol,'omitnan')';
ds.dat.task.avg_av = mean(ds.dat.y_av,'omitnan')';

%% plot ASE-inspired readouts
metac_plot_ase_readouts_raw(ds.dat);

%% load quest data (discovery set)
dataset = 1; % 1=discovery set
% quest = load_quest_data(dataset);
load(fullfile('data', 'discovery_set_quest_tmp.mat'));
ds.quest = quest;

%% create table for winning mod 1 (DISCOVERY SET)
m = ds.gest.ffx.idx;
ds_pars = ds.gest.param_mat';
if m == 4 && size(ds_pars,2) == 5
    ds.tab = table(ds.quest.y_fas, ds.quest.age, ds.quest.gender, ...
        ds_pars(:,2), ds_pars(:,3), ds_pars(:,4), ds_pars(:,5), ...
        ds.dat.task.avg_c, ds.dat.task.avg_tol, ds.dat.task.avg_av, ...
        ds.quest.y_mfis);
    ds.tab.Properties.VariableNames = {'FAS', 'age', 'gender',...
        'gamma', 'shift', 'scale', 'w',...
        'contr', 'tol', 'av',...
        'MFIS'};
end
% normalize values
ds.tab_norm = normalize(ds.tab);

%% write csv file with data for Bayesian ANCOVA
save_path = fullfile('data', ['tmp_ds_data.csv']);
writetable(ds.tab, save_path);

save_path2 = fullfile('data', ['tmp_norm_ds_data.csv']);
writetable(ds.tab_norm, save_path2);

%% Correlations (DISCOVERY SET)
% plot
figure
% plotmatrix(table2array(tab))
gplotmatrix(table2array(ds.tab),[],[],[],[],[],[],[],ds.tab.Properties.VariableNames)
% calc
[coef, pval] = corr(table2array(ds.tab));
pval<0.05
figdir = fullfile('figures', ['scatterplotmat_ds']);
print(figdir, '-dpng');
close;


%% estimate empirical priors
[vs.mod, ep.mod] = metac_est_ep(ds.dat, ds.mod);

%% create synthetic data
n_pe = 1; % PE traj from pilot sub n_pe
n_sim = 100;
[ep.sim] = metac_sim_ep(ds.dat, ds.pdat, vs.mod, n_pe, n_sim);

%% fit simulated data

% loop over synthetic subjects & model space
for n = 1:n_sim
    for m = 1:size(ds.mod,2) % sim model
        fprintf('current iteration: n=%1.0f, m=%1.0f \n', n,m);
        for i = 1:size(ds.mod,2) %est model
            rec.est(m,n,i).data = tapas_fitModel(ep.sim.sub(n,m).data.y,...
                [ds.dat.u_bin(:,n_pe) ds.pdat.u_pe(:,n_pe)],...
                vs.mod(m).prc_config,...
                vs.mod(m).obs_config,...
                tapas_quasinewton_optim_config ...
                );
            % store LME in matrix
            rec.model(m).LME(n,i) = rec.est(m,n,i).data.optim.LME;
        end
        % store parameter values (sim & est)
        % rec.param.prc(m).sim(n,:) = ep.sim.sub(n,m).input.prc.transInp(vs.mod(m).prc_idx);
        rec.param.obs(m).sim(n,:) = ep.sim.sub(n,m).input.obs.transInp(vs.mod(m).obs_idx);
        % rec.param.prc(m).est(n,:) = rec.est(m,n,m).data.p_prc.ptrans(vs.mod(m).prc_idx);
        rec.param.obs(m).est(n,:) = rec.est(m,n,m).data.p_obs.ptrans(vs.mod(m).obs_idx);
    end

end

%% recovery analysis
[ep.rec] = metac_rec_analysis(vs, rec, n_sim);
save(fullfile('data', 'empirical_priors_rec.mat'), 'ep', '-mat');

%% load validation set task
dataset = 2; % 2=validation set
% dat = load_task_data(dataset);
load(fullfile('data', 'validation_set_tmp.mat'));
vs.dat = dat;

%% preprocess task data
vs.pdat = metac_preproc(vs.dat);

%% fit behavioural models using empirical priors
priors = 2; % empirical priors
[vs.gest, vs.mod] = metac_fit_behaviour(vs.dat, vs.pdat, vs.mod, priors);
% dat = metac_fit_behaviour2(dat);
save(fullfile('data', 'validation_set_fits_tmp.mat'), 'vs', '-mat');

%% avg overall control, tolerance, ...
vs.dat.task.avg_c = mean(vs.dat.y_c,'omitnan')';
vs.dat.task.avg_tol = mean(vs.dat.y_tol,'omitnan')';
vs.dat.task.avg_av = mean(vs.dat.y_av,'omitnan')';

%% plot ASE-inspired readouts
metac_plot_ase_readouts_raw(vs.dat);

%% load quest data
dataset = 2; % 2=validation set
% quest = load_quest_data(dataset);
load(fullfile('data', 'validation_set_quest_tmp.mat'));
vs.quest = quest;


%% create table for winning mod 1 (VALIDATION SET)
m = vs.gest.ffx.idx;
vs_pars = vs.gest.param_mat';
if m == 4 && size(vs_pars,2) == 5
    vs.tab = table(vs.quest.y_fas, vs.quest.age, vs.quest.gender, ...
        vs_pars(:,2), vs_pars(:,3), vs_pars(:,4), vs_pars(:,5), ...
        vs.dat.task.avg_c, vs.dat.task.avg_tol, vs.dat.task.avg_av, ...
        vs.quest.y_mfis);
    vs.tab.Properties.VariableNames = {'FAS', 'age', 'gender',...
        'gamma', 'shift', 'scale', 'w',...
        'contr', 'tol', 'av',...
        'MFIS'};
end
% normalize values
vs.tab_norm = normalize(vs.tab);


%% write csv file with data for Bayesian ANCOVA
save_path = fullfile('data', ['tmp_vs_data.csv']);
writetable(vs.tab, save_path);

save_path2 = fullfile('data', ['tmp_norm_vs_data.csv']);
writetable(vs.tab_norm, save_path2);



%% Correlations (VALIDATION SET)

% plot
figure
% plotmatrix(table2array(tab))
gplotmatrix(table2array(vs.tab),[],[],[],[],[],[],[],vs.tab.Properties.VariableNames)
% calc
[coef, pval] = corr(table2array(vs.tab));
pval<0.05
% save
figdir = fullfile('figures', ['scatterplotmat_vs']);
print(figdir, '-dpng');
close;



% %% tol vs av
% figure
% scatter(ds.tab.tol, ds.tab.av)
% xlabel('tolerance')
% ylabel('aversiveness')
% 
% %% FAS vs MFIS
% figure
% plot(ds.tab.MFIS, ds.tab.FAS, '.')
% xlabel('MFIS')
% ylabel('FAS')