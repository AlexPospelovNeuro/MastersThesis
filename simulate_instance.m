% This is a "shell" function that does the simulation of the circuit
% activity under the given circumstances.

% INPUTS:

% - An instance of the neural circuit.
% - A maximal duration of the simulation
% - (optional) a rng seed for the simulation
% - An input organiser function, with its own arguments possibly.
% - (optional) An environmental feedback function, with its own arguments possibly.

% OUTPUTS:

% - Raster of activity.
% - The markdown of actually presented inputs
% - (optional) additional outputs and flags such as an early termination flag 


% 02.02.2025 update: Now Raster can be provided as a second input variable
% of the function. In NaN, the function creates a default blanc raster. Can
% be used, for example, to silence the neuron subclasses


function [Raster, input_markdown, extras]=simulate_instance(Inst,Raster,Runtime,tech,maxdur,seed,input,feedback)

input_markdown=[];
extras=[];
% Preparations: creating runtimes, determining initial delay etc


% Preparation: construction of simulation-related variables.

if ~isnan(seed)
    rng(seed, 'twister');
end
tech.totaldur=maxdur+tech.Longest_shape;
if isnan(Raster)
    Raster(1:Inst.N_cells_total,1:tech.totaldur)=NaN;
end
Input_pointer=1;
for I=1:size(input,2)
    Raster(Input_pointer:Input_pointer+Inst.Source.Cells.Input(1,I)-1,tech.Longest_shape+1:tech.totaldur)=input{I};
    Input_pointer=Input_pointer+Inst.Source.Cells.Input(1,I);
end
Raster(:,1:tech.Longest_shape)=0;

Input_vector(1:Inst.N_cells_nonmod)=0; % rewritable placeholder for evaluation of the inputs of the cell
Plasticity_Coeff=1; % A rewritable plasticity coefficient initialization
Mod_Coeff=1; % A rewritable modulation coefficient initialization 

Runtime.DynamicThrsh(1:Inst.N_cells_total,1:tech.totaldur)=NaN;
Runtime.DynamicThrsh(1:Inst.N_cells_total,tech.Longest_shape+1:end)=repmat(Runtime.Thrsh',1,maxdur);
for neu=1:Inst.N_cells_total
    Runtime.DynamicThrsh(neu,tech.Longest_shape+1:end)=Runtime.DynamicThrsh(neu,tech.Longest_shape+1:end)+randn(1,maxdur)*Runtime.ThrshNoise(neu);
end

for t=tech.Longest_shape+1:tech.totaldur
    %[Raster,Runtime]=step_update(Raster,Runtime,Inst,t);
    
    % 22.07.2025 A version for less data trafic between the functions
    Slice=Raster(:,t-tech.Longest_shape:t);
    [Upd,Runtime]=step_update_v2(Slice,Runtime,Inst,t);
    Raster(:,t)=Upd;
    %[Raster]=env_feedback(Raster,func);
end



end