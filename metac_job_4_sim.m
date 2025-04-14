function [] = metac_job_4_sim(EULER)
% [] = metac_job_4_sim(EULER)
%
% Create synthetic data for the METAC Task under the empirical priors.
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

%% load ds + vs
tmp1 = load(fullfile('results', 'ds', 'discovery_set_fits_tmp.mat'));
ds = tmp1.ds;
tmp2 = load(fullfile('results', 'vs', 'mod_space.mat'));
vs = tmp2.vs;

%% create synthetic data
disp(ds)
[ep.sim] = metac_sim_ep(ds.dat, ds.pdat, vs.mod, opts.sim.n_pe, opts.sim.n_sim, opts);

%% save simulated data
save_path = fullfile(sdir, 'results', 'sim', ['ep_sim_data_set']);
save(save_path, '-struct', 'ep');

disp('simulated data set successfully created.')

end