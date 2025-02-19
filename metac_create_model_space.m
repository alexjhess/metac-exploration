function [mod, bo] = metac_create_model_space()
% [ModelSpace] = metac_create_model_space()
%
% Creates and outputs a struct of the model space containing config files 
% for perceptual and response models.
%
% INPUT
%   argin            type           
%
%   OPTIONAL:
%
% OUTPUT    
%   ModelSpace       struct      Struct containing prc & obs model configs
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

%% init struct
mod = struct();

%% Bayes Optimal (dummy obs model)
bo.prc_config = tapas_hgf_binary_mab_2mu0_metac_config;
bo.prc_config.n_bandits = 2;
bo.obs_config = tapas_bayes_optimal_binary_config;

%% Perceptual Model
% 3-level binary HGF mab 2mu0
for m = 1:4
    mod(m).prc_config = tapas_hgf_binary_mab_2mu0_metac_config;
    mod(m).prc_config.n_bandits = 2;
end

%% Response Models

%--------------------------------------------------------------------------
% combining y_int and y_mc

% M1: Null model
mod(1).obs_name = 'null';
mod(1).obs_config = mc_null_inobs_mab_config;

% M2: AR(k) model (Res)
mod(2).obs_name = 'AR(k) Res';
mod(2).obs_config = mc_autoreg_res_inobs_config;

% M3: AR(k) model (PE)
mod(3).obs_name = 'AR(k) PE';
mod(3).obs_config = mc_autoreg_pe_inobs_mab_config;

% M4: AR(k) model (PE + RES)
mod(4).obs_name = 'AR(k) PE + Res';
mod(4).obs_config = mc_autoreg_pe_res_inobs_mab_config;


%% find free obs parameters

% prc
prc_idx = bo.prc_config.priorsas;
prc_idx(isnan(prc_idx)) = 0;
bo.prc_idx = find(prc_idx);

% obs
for m = 1:size(mod,2)
    % prc
    prc_idx = mod(m).prc_config.priorsas;
    prc_idx(isnan(prc_idx)) = 0;
    mod(m).prc_idx = find(prc_idx);
    % obs
    obs_idx = mod(m).obs_config.priorsas;
    obs_idx(isnan(obs_idx)) = 0;
    mod(m).obs_idx = find(obs_idx);
end


end