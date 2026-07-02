


% 15.7.2025 An error of modulatory neurons addressation fixed.
% 22.7.2025 A version of the function that takes in only a relevant part f
% the Slice (t-tech.longest:t-1) and returns only 0 or 1. 

function [Out_spike, Runtime]=single_cell_update_v3_1(Slice,Runtime,Inst,cell)
t=size(Slice,2);
Runtime.AbsRefCounter(cell)=max(0,Runtime.AbsRefCounter(cell)-1); % Reset of the Absolute refracterity

cell_cl=Inst.Cellclass(cell); % The subclass of the current (postsynaptic) neuron
Input_vector(1:Inst.N_cells_nonmod)=0; % The placeholder for the input vector of the neuron 
        
% The tracker shifts and output-independent plasticity triggering do not depend on the absolute refracterity, they are doen no matter what.

RR_Coeff=1; % A coefficient for the relative refracterity
for Plst_pat=1:size(Inst.RelRefract{cell_cl},2) % And every plasticity pattern of the type
    Runtime.RefTrack{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)}(:)=[Runtime.RefTrack{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)}(2:end) 0]; % Shifting the tracker of relative refracterity by 1
    RR_Coeff=RR_Coeff*(max(Runtime.RefTrack{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)}(:,1),-1)+1); % The effect of multuple patterns/types is multiplicative

end
for inp=1:Inst.N_cells_nonmod % Sequentually evaluate the possible input from each other neuron (if no connection, it will be zero)
    Plasticity_Coeff=1; % Reset the plasticity coeff to 1: there may be no plasticity at all (by default)
    Mod_Coeff=1; % Reset        
    if  Inst.Connection_matrix(cell,inp)==1
       post_cl=Inst.Synapses.IDs(Inst.ConID(cell,inp),2); % subclass of the postsynaptic cell (should be identical to cell_cl, is it?) 
       pre_cl=Inst.Synapses.IDs(Inst.ConID(cell,inp),3); % subclass of the presynaptic cell
       syn=Inst.Synapses.IDs(Inst.ConID(cell,inp),6); % ID of the synapse within these subclasses' ordered pair
       % Shift of the plasticity trackers, 
       for Plst_type=1:size(Inst.Synapses.Plast{post_cl,pre_cl},2) %For every type of the plasticity in the synapse
           for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type},2) % And every plasticity pattern of the type
               Runtime.PlastTrack{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}=[Runtime.PlastTrack{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}(2:end) 0]; % Make a shift by one in the timescale
               Plasticity_Coeff=Plasticity_Coeff*(max(Runtime.PlastTrack{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}(:,1),-1)+1); % The the effect of multuple patterns/types is multiplicative

           end
       end
       % Modulation 
       for mod_cl=1:Inst.N_mod_subclasses % for each class of the modulation 
           if ~isempty(Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl})
                mod_inputs=Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}((Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(:,1)==Inst.Cell_inclass(cell))&(Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(:,2)==Inst.Cell_inclass(inp)),3)';
                mod_index=Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}((Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(:,1)==Inst.Cell_inclass(cell))&(Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(:,2)==Inst.Cell_inclass(inp)),4)';
                for m_i=1:size(mod_inputs,2)
                    Runtime.ModTrack{post_cl,pre_cl,mod_cl}{mod_index(m_i)}=[Runtime.ModTrack{post_cl,pre_cl,mod_cl}{mod_index(m_i)}(2:end) 0];
                    Mod_Coeff=Mod_Coeff*(max(Runtime.ModTrack{post_cl,pre_cl,mod_cl}{mod_index(m_i)}(:,1),-1)+1);

                end
            end
       end 
       
       Input_vector(inp)=RR_Coeff*Mod_Coeff*Plasticity_Coeff*Inst.Synapses.Powers{post_cl,pre_cl}(syn)*Inst.Connection_matrix(cell,inp)*sum(Slice(inp,t-Runtime.PSPlength{post_cl,pre_cl}{syn}-Inst.Synapses.Delays{post_cl,pre_cl}(syn):t-Inst.Synapses.Delays{post_cl,pre_cl}(syn)-1).*Runtime.PSP{post_cl,pre_cl}{syn});

       % THE OUTPUT-INDEPENDENT PLASTICITY, TRIGGERS NO MATTER IF THE NEURON SPIKES
       for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{1},2) % For each plasticity pattern of the type 1 in the present pair of subclasses
 % >>>>>>>>>>>>>>>>>>>> DELAY                         % For the plasticity of the type 1, I must check the observation frame for spike in the corresponding input neuron, taking the delay into account
           if (Slice(inp,t-Inst.Synapses.Delays{post_cl,pre_cl}(syn))==1)
               Runtime.PlastTrack{post_cl,pre_cl}{1}{Plst_pat}{syn}=tracker_update(Runtime.PlastTrack{post_cl,pre_cl}{1}{Plst_pat}{syn}, Runtime.Plast{post_cl,pre_cl}{1}{Plst_pat}{syn},Runtime.PlastN{post_cl,pre_cl}{1}{Plst_pat}{syn}); % Shift all the available pattens by one, so the "the oldest" expire (if there is an "overflow")
           end
       end
    end
end        

% The output-independent update of trackers is done. Now to decision about spiking
if Runtime.AbsRefCounter(cell)==0 % If the neuron is not refractory 
    if isnan(Slice(cell,t))&&(sum(Input_vector)>Runtime.DynamicThrsh(cell,t)) % The spike, if the triggering condition is fulfilled
        Out_spike=1; % The spike
        Runtime.AbsRefCounter(cell)=Runtime.AbsRefrRef(cell)+1; % Reset absolute refracterity
        % The activation of plasticity patterns if applicable (except type 1). 
        
        % Type 1000, relative refracterity
        for Plst_pat=1:size(Inst.RelRefract{cell_cl},2) % For each plasticity pattern of the type 1000
            Runtime.RefTrack{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)}=tracker_update(Runtime.RefTrack{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)},Runtime.RelRef{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)},Runtime.RefN{cell_cl}{Plst_pat}{Inst.Cell_inclass(cell)}); % Shift all the available pattens by one, so the "the oldest" expire (if there is an "overflow")
        end
        for inp=1:Inst.N_cells_nonmod % Sequentually evaluate the possible input from each other neuron (if no connection, it will be zero)
            if  Inst.Connection_matrix(cell,inp)==1
                post_cl=Inst.Synapses.IDs(Inst.ConID(cell,inp),2);
                pre_cl=Inst.Synapses.IDs(Inst.ConID(cell,inp),3);
                syn=Inst.Synapses.IDs(Inst.ConID(cell,inp),6);
                
                % Type 2, spike-timing plasticity
                for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{2},2) % For each plasticity pattern of the type 2 in the present pair of subclasses
                    % For the plasticity of the type 2, I must check the observation frame for the previous spike of the same neuron
                    if sum(Slice(cell,t-Inst.Synapses.Plast{post_cl,pre_cl}{2}{Plst_pat}{1}(2,syn):t-Inst.Synapses.Plast{post_cl,pre_cl}{2}{Plst_pat}{1}(1,syn)))>0
                        Runtime.PlastTrack{post_cl,pre_cl}{2}{Plst_pat}{syn}=tracker_update(Runtime.PlastTrack{post_cl,pre_cl}{2}{Plst_pat}{syn},Runtime.Plast{post_cl,pre_cl}{2}{Plst_pat}{syn},Runtime.PlastN{post_cl,pre_cl}{2}{Plst_pat}{syn});
                    end
                end
                % Type 3, hebbian plasticity
                for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{3},2) % For each plasticity pattern of the type 3 in the present pair of subclasses
 % >>>>>>>>>>>>>>>>>>>> DELAY                           % For the plasticity of the type 3, I must check the observation frame for spike in the corresponding input neuron, taking the delay into account
                    if sum(Slice(inp,t-Inst.Synapses.Plast{post_cl,pre_cl}{3}{Plst_pat}{1}(2,syn)-Inst.Synapses.Delays{post_cl,pre_cl}(syn):t-Inst.Synapses.Delays{post_cl,pre_cl}(syn)-Inst.Synapses.Plast{post_cl,pre_cl}{3}{Plst_pat}{1}(1,syn)))>0
                        Runtime.PlastTrack{post_cl,pre_cl}{3}{Plst_pat}{syn}=tracker_update(Runtime.PlastTrack{post_cl,pre_cl}{3}{Plst_pat}{syn},Runtime.Plast{post_cl,pre_cl}{3}{Plst_pat}{syn},Runtime.PlastN{post_cl,pre_cl}{3}{Plst_pat}{syn});
                    end
                end
            end
        end
        % The modulation, if applicable
        if cell>Inst.N_cells_vector_incremental(Inst.N_nonmod_subclasses) % The observed neuron is modulatory 
            if ~isempty(Inst.Synapses.Mod_targets{Inst.Cellclass(cell)-Inst.N_nonmod_subclasses}{Inst.Cell_inclass(cell)}) % If the current modulatory neuron has any targets at all
                for mod_target=1:size(Inst.Synapses.Mod_targets{Inst.Cellclass(cell)-Inst.N_nonmod_subclasses}{Inst.Cell_inclass(cell)}(:,:),1)
                    post_cl=Inst.Synapses.Mod_targets{Inst.Cellclass(cell)-Inst.N_nonmod_subclasses}{Inst.Cell_inclass(cell)}(mod_target,1);
                    pre_cl=Inst.Synapses.Mod_targets{Inst.Cellclass(cell)-Inst.N_nonmod_subclasses}{Inst.Cell_inclass(cell)}(mod_target,3);
                    mod_cl=Inst.Cellclass(cell)-Inst.N_nonmod_subclasses;
                    [targ_loc, ~]=ismember([Inst.Synapses.Mod_targets{mod_cl}{Inst.Cell_inclass(cell)}(mod_target,2)  Inst.Synapses.Mod_targets{mod_cl}{Inst.Cell_inclass(cell)}(mod_target,4) Inst.Cell_inclass(cell)], Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(:,1:3), 'rows');
                    Runtime.ModTrack{post_cl,pre_cl,mod_cl}{Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(targ_loc,4)}=tracker_update(Runtime.ModTrack{post_cl,pre_cl,mod_cl}{Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(targ_loc,4)},Runtime.Mod{post_cl,pre_cl,mod_cl}{Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(targ_loc,4)},Runtime.ModN{post_cl,pre_cl,mod_cl}{Inst.Synapses.Mod_source{post_cl,pre_cl,mod_cl}(targ_loc,4)});
                end
            end
        end
    else
        Out_spike=0;
    end 
else
    %Slice(cell,t)=0;
    Out_spike=0;
end    
end


function tracker=tracker_update(tracker,update,N)

    if sum(update)<0 % Negaticve effect on power
        tracker=max(tracker+update,min(update)*N); 
    else % Positive effect on power
        tracker=min(tracker+update,max(update)*N); 
    end

    
end