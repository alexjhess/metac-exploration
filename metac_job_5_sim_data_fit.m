function [] = metac_job_5_sim_data_fit(EULER, n, m, i)
% [] = job_runner_5_sim_data_modinv(EULER, n, m, i)
%
% Inverts model i on simulated subject n created with model m under the
% empirical priors.
%
% INPUT
%   EULER        binary           Binary indicator variable
%   n            integer          Integer indicating synthetic subject idx
%   m            integer          Integer for model used to create data
%   i            integer          Integer for model used to fit data
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

%% load ds + vs
tmp1 = load(fullfile('results', 'ds', 'discovery_set_fits_tmp.mat'));
ds = tmp1.ds;
tmp2 = load(fullfile('results', 'vs', 'mod_space.mat'));
vs = tmp2.vs;
ep = load(fullfile(sdir, 'results', 'sim', ['ep_sim_data_set']));

%% model inversion
% seed for multistart optim
opts.opt_config.seedRandInit = opts.rng.settings.State(opts.rng.idx, 1);
% fit model
est = tapas_fitModel(ep.sim.sub(n,m).data.y,... % simulated mc responses
                [ds.dat.u_bin(:,opts.sim.n_pe) ds.pdat.u_pe(:,opts.sim.n_pe)],... % inputs (R/nR + PE)
                vs.mod(m).prc_config,... % prc model
                vs.mod(m).obs_config,... % obs model
                opts.opt_config ... % opt algo
                );

%% save model fit as struct
save_path = fullfile(sdir, 'results', 'sim', ['sub', num2str(n)],...
    ['sim_mod', num2str(m), '_est_mod', num2str(i)]);
save(save_path, '-struct', 'est');

end