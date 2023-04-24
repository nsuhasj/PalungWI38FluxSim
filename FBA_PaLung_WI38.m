%{ 
Copyright (C) 2023  N. Suhas Jagannathan
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
%}

% This script is used to run flux sampling simulations for the bat PaLung
% nad human WI38 mitochondrial metabolic models using CobraT toolbox v3.

%% Initialize models and simulation parameters 

options.nStepsPerPoint = 200;
options.nPointsReturned = 5000;
options.optPercentage = 0;

b_model = readCbModel('PaLungMetModel.xlsx');
h_model = readCbModel('WI38MetModel.xlsx');
b_model = changeObjective(b_model,'CV_MitoCore');
h_model = changeObjective(h_model,'CV_MitoCore');

%% Control flux sampling simulations (no constraints)

[~, sample_bat_ctrl] = sampleCbModel(b_model,'sample_bat_ctrl.mat','CHRR',options);
[~, sample_human_ctrl] = sampleCbModel(h_model,'sample_human_ctrl.mat','CHRR',options);

%% Set constraints on Complex I and Oxygen intake into mitochondria in both models.

[b_minC1, b_maxC1] = fluxVariability(b_model,0,'max',{'CI_MitoCore'});
b_C1_70prct = b_minC1 + (b_maxC1-b_minC1)*0.7;
b_model = changeRxnBounds(b_model,'CI_MitoCore',b_C1_70prct,'l');

[h_minC1, h_maxC1] = fluxVariability(h_model,0,'max',{'CI_MitoCore'});
h_C1_30prct = h_minC1 + (h_maxC1-h_minC1)*0.3; 
h_model = changeRxnBounds(h_model,'CI_MitoCore',h_C1_30prct,'u'); 

[b_minO2, b_maxO2] = fluxVariability(b_model,0,'max',{'O2tm'});
b_O2_30prct = b_minO2 + (b_maxO2-b_minO2)*0.3; 
b_model = changeRxnBounds(b_model,'O2tm',b_O2_30prct,'u'); 

[h_minO2, h_maxO2] = fluxVariability(h_model,0,'max',{'O2tm'});
h_O2_70prct = h_minO2 + (h_maxO2-h_minO2)*0.7;  
h_O2_70prct = max(h_O2_70prct,b_O2_30prct);
h_model = changeRxnBounds(h_model,'O2tm',h_O2_70prct,'l'); 

%% Flux sampling simulations with constraints

[~, sample_bat_constr] = sampleCbModel(b_model,'sample_bat_constr.mat','CHRR',options);
[~, sample_human_constr] = sampleCbModel(h_model,'sample_human_constr.mat','CHRR',options);
