function [sdir] = metac_paths(EULER)
% [sdir] = metac_paths(EULER)
%
% This function adds all the necessary subfolders to the path and
% outputs the directory where the results & figures folders are stored.
%
% INPUT
%   EULER        binary       Binary indicator variable  
%
%   OPTIONAL:
%
% OUTPUT    
%   sdir         string       Path where res & figs folders are stored
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

%% directory where results & figures folders are stored
sdir = pwd;

%% setup path
if EULER == 1
    % add subfolders
    addpath('comb_obs_models');
    addpath('mc_obs');
    
    % add submodules (toolboxes)
    addpath(genpath('hgf-toolbox'));
    cd('VBA-toolbox');
    VBA_setup();
    cd ..

    sdir= '/cluster/work/tnu/alhess/spirl_slurm_normal';
end

end
