function [] = metac_job_2_ds_fit(EULER, n, m)
% [] = metac_job_2_ds_fit(EULER, n, m)
%
% Inverts model m on data from participant n of the discovery set.
%
% INPUT
%   EULER        binary           Binary indicator variable
%   n            integer          Integer indicating participant index
%   m            integer          Integer indicating model index
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

%% load discovery set & model space
ds = load(fullfile(sdir, 'results', 'ds', ['init_mod_space']));

%% model inversion
% seed for multistart optim
opts.opt_config.seedRandInit = opts.rng.settings.State(opts.rng.idx, 1);
est = tapas_fitModel(ds.pdat.y_mc(:,n),... % mc responses
            [ds.dat.u_bin(:,n) ds.pdat.u_pe(:,n)],... % inputs (succ/fail, PE)
            ds.ip.mod(m).prc_config,... % prc model
            ds.ip.mod(m).obs_config,... % obs model
            opts.opt_config ... % opt algo
            );

%% save model fit as struct
save_path = fullfile(sdir, 'results', 'ds', ['sub', num2str(n)],...
    ['est_mod', num2str(m)]);
save(save_path, '-struct', 'est');

fprintf('save fits (iteration: n=%1.0f, m=%1.0f) \n', n,m);

end