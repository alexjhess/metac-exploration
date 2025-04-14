function [] = meatc_job_3_ep(EULER)
% [] = meatc_job_3_ep(EULER)
%
% Estimates sufficient statistics of empirical prior distributions based on MAP
% estimates obtained from model inversion on held-out (discovery) data set.
%
% INPUT
%   EULER        binary           Binary indicator variable
%
%   OPTIONAL:
%
% OUTPUT    
%   argout       type
%
% _________________________________________________________________________
% Author: Alex Hess
%
% Copyright (C) 2025 Translational Neuromodeling Unit
%                    Institute for Biomedical Engineering
%                    University of Zurich & ETH Zurich
%
% This file is released under the terms of the GNU General Public Licence
% (GPL), version 3. You can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% This file is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <https://www.gnu.org/licenses/>.
% _________________________________________________________________________

%% setup path
sdir = metac_paths(EULER);

%% load analysis specifications
opts = load(fullfile(sdir, 'results', ['options']));

%% load data & init model space
ds = load(fullfile(sdir, 'results', 'ds', ['init_mod_space']));

%% load results from model inversion
gest.F = NaN(size(ds.ip.mod,2),size(ds.dat.u_bin,2));
gest.Ll = NaN(size(ds.ip.mod,2),size(ds.dat.u_bin,2));
gest.comp = NaN(size(ds.ip.mod,2),size(ds.dat.u_bin,2));

for m = 1:size(ds.ip.mod,2)
    % init obs param_mat (size n_obs_pars x N)
    ds.ip.mod(m).param_mat = NaN(size(ds.ip.mod(m).obs_idx,2),...
            size(ds.dat.u_pe,2));
    for n = 1:opts.ds.nS
        fprintf('load fits (iteration: n=%1.0f, m=%1.0f) \n', n,m);
        ds.ip.mod(m).sub(n).est = load(fullfile(sdir, 'results',...
            'ds', ['sub', num2str(n)], ['est_mod', num2str(m)]));

        gest.F(m,n) = ds.ip.mod(m).sub(n).est.optim.LME;
        gest.Ll(m,n) = ds.ip.mod(m).sub(n).est.optim.accu;
        gest.comp(m,n) = ds.ip.mod(m).sub(n).est.optim.comp;
        ds.ip.mod(m).param_mat(:,n) = ds.ip.mod(m).sub(n).est.p_obs.p(ds.ip.mod(m).obs_idx);
             
    end
end

%% DS: model comp

% FFX BMS
[gest.ffx.sumLME, gest.ffx.pp, gest.ffx.GBF, gest.ffx.ABF] = metac_FFXBMS(gest.F');
[gest.ffx.val, gest.ffx.idx] = max(gest.ffx.sumLME);
disp('FFX BMS results: ')
fprintf('winning model %i \n', gest.ffx.idx)
if ~isempty(gest.ffx.GBF)
    fprintf('GBF %i \n', gest.ffx.GBF)
end
if ~isempty(gest.ffx.ABF)
    fprintf('ABF %i \n', gest.ffx.ABF)
end

% RFX BMS
[gest.rfx.posterior, gest.rfx.out] = VBA_groupBMC(gest.F);
[gest.rfx.val, gest.rfx.idx] = max(gest.rfx.out.pxp);
disp('RFX BMS results: ')
fprintf('winning model %i \n', gest.rfx.idx)
fprintf('PXP %.2f \n', gest.rfx.out.pxp(gest.rfx.idx))
fprintf('Ef %.2f \n', gest.rfx.out.Ef(gest.rfx.idx))
% save fig
figdir = fullfile('figures', 'ds', 'ip',...
    ['rfx_bms']);
print(figdir, '-dpng');
close;

% extract est params of winning model (npars x nsub)
gest.param_mat = ds.ip.mod(gest.ffx.idx).param_mat;

%% store vals
ds.gest = gest;
ds.mod = ds.ip.mod;

%% avg + var overall control, tolerance, ...
ds.dat.task.avg_c = mean(ds.dat.y_c,'omitnan')';
ds.dat.task.avg_tol = mean(ds.dat.y_tol,'omitnan')';
ds.dat.task.avg_av = mean(ds.dat.y_av,'omitnan')';
ds.dat.task.var_c = var(ds.dat.y_c,'omitnan')';
ds.dat.task.var_tol = var(ds.dat.y_tol,'omitnan')';

%% load quest data (discovery set)
% dataset = 1; % 1=discovery set
% quest = load_quest_data(dataset);
load(fullfile('data', 'discovery_set_quest_tmp.mat'));
ds.quest = quest;

%% create table for winning mod (DISCOVERY SET)
m = ds.gest.ffx.idx;
ds_pars = ds.gest.param_mat';
if m == 4 && size(ds_pars,2) == 5
    ds.tab = table(ds.quest.y_fas, ds.quest.y_mfis, ds.quest.age, ds.quest.gender, ...
        ds_pars(:,2), ds_pars(:,3), ds_pars(:,4), ds_pars(:,5), ...
        ds.dat.task.avg_c, ds.dat.task.avg_tol, ds.dat.task.avg_av, ...
        ds.dat.task.var_c, ds.dat.task.var_tol, ...
        ds.quest.maia3, ds.quest.maia8, ds.quest.maia38, ds.quest.psqi);
    ds.tab.Properties.VariableNames = {'FAS', 'MFIS', 'age', 'gender',...
        'gamma', 'shift', 'scale', 'w',...
        'contr', 'tol', 'av',...
        'var_c', 'var_tol',...
        'MAIA3', 'MAIA8', 'MAIA38', 'PSQI'};
end

%% normalize values (only features)
ds.tab_norm = ds.tab;
ds.tab_norm(:,3:end) = normalize(ds.tab(:,3:end));

%% write csv file with data for Bayesian ANCOVA
save_path = fullfile('results', 'sq2', ['tmp_ds_data.csv']);
writetable(ds.tab, save_path);

save_path2 = fullfile('results', 'sq2', ['tmp_norm_ds_data.csv']);
writetable(ds.tab_norm, save_path2);

%% Correlations (DISCOVERY SET)
[ds.corr.coef, ds.corr.pval] = corr(table2array(ds.tab));

%% estimate empirical priors
[vs.mod, ep.mod] = metac_est_ep(ds.dat, ds.ip.mod);


%% save results
save(fullfile('results', 'ds', 'discovery_set_fits_tmp.mat'), 'ds', '-mat');
save(fullfile('results', 'vs', 'mod_space.mat'), 'vs', '-mat');

disp('empricial priors successfully estimated.')
disp('ready for simulation analyses and fitting of validation set.')

end