function [opts] = metac_opts()
% [opts] = metac_opts()
%
% This function saves all the settings for running the analysis pipeline to
% a struct.
%
% INPUT
%   argin        type           
%
%   OPTIONAL:
%
% OUTPUT    
%   opts         struct       Struct with specification of the analysis
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

% pilot data set
opts.ds.nS = 20;
opts.ds.fname = ['discovery_set_tmp.mat'];

% simulation analysis
opts.sim.n_sim = 100;
opts.sim.n_pe = 1; % PE traj from discovery set sub n_pe

% main data set
opts.vs.nS = 30;
opts.vs.filename = ['metac_vs_',...
    num2str(opts.vs.nS), 'sub_preprocessed_IncludedInAnalysis'];

% posterior predictive checking: number of draws
opts.vs.ppc.nDraws_per_sub = 100;

% optimization algorithm
opts.opt_config = eval('tapas_quasinewton_optim_config');
opts.opt_config.nRandInit = 0%399;

% seed for rng
rng(123, 'twister')
opts.rng.settings = rng;
opts.rng.idx = 1; % Set counter for random number states

% define colors for plotting
opts.col.wh = [1 1 1];
opts.col.gry = [0.5 0.5 0.5];
opts.col.tnub = [0 110 182]/255; 
opts.col.tnuy = [255 166 22]/255;
opts.col.grn = [0 0.6 0];

end