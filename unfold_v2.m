
% A new iteration of the unfold function.
% The Runtime trackers for modulations and plasticites are significantly
% changed compared to the previous unfold. Runtime generation is transfered
% to a separate function so it can be done independently from unfolding

% 15.7.2025 An error of modulatory neurons addressation fixed.

function [Inst,Runtime,tech]=unfold_v2(Net2,PSP_db,Plast_db,seed)


%%%%%%%%%%%%%%%%% Generation of an instance %%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% The unfolding routine %%%%%%%%%%%%%%%%%
% The variable Inst will contain all the relevant information about the
% instance. 
% The unfolding will happen in stages, to ensure the correct transferring of
% the genetic parameters into the features of the instance. 

% Tentatively, for the routine when one genotype leads to one instance, so
% the distribution parameters for the cell numbers are irrelevant, I simply
% round the numbers and bring below 1 valueas to 1.

if ~isnan(seed)
    rng(seed, 'twister');
end

Net2.Cells.Input(1,:)=max(round(Net2.Cells.Input(1,:)),1);
Net2.Cells.Output(1,:)=max(round(Net2.Cells.Output(1,:)),1);
Net2.Cells.Ion(1,:)=max(round(Net2.Cells.Ion(1,:)),1);
Net2.Cells.Mod(1,:)=max(round(Net2.Cells.Mod(1,:)),1);

% numbers of neurons per subclass. In this test, the one-instance case is
% used, therefore the sigma of cell classes number does not matter (the mean is used always)
Inst.Source=Net2; % The genetic code is saved inside the instance
Inst.N_input_cells=sum(Net2.Cells.Input(1,:));
Inst.N_subclasses=size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2)+size(Net2.Cells.Ion,2)+size(Net2.Cells.Mod,2);
Inst.N_nonmod_subclasses=size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2)+size(Net2.Cells.Ion,2);
Inst.N_mod_subclasses=size(Net2.Cells.Mod,2);
Inst.Mod_pointer=Inst.N_nonmod_subclasses+1:Inst.N_nonmod_subclasses+Inst.N_mod_subclasses;
Inst.N_cells_total=sum(Net2.Cells.Input(1,:))+sum(Net2.Cells.Output(1,:))+sum(Net2.Cells.Ion(1,:))+sum(Net2.Cells.Mod(1,:));
Inst.N_cells_nonmod=sum(Net2.Cells.Input(1,:))+sum(Net2.Cells.Output(1,:))+sum(Net2.Cells.Ion(1,:));
Inst.N_cells_vector=[Net2.Cells.Input(1,:) Net2.Cells.Output(1,:) Net2.Cells.Ion(1,:) Net2.Cells.Mod(1,:)];
Inst.N_cells_vector_incremental=cumsum(Inst.N_cells_vector);
Inst.N_cells_vector_incremental_pointer=[1 Inst.N_cells_vector_incremental+1];
Inst.Cellclass=[];

for cl=1:Inst.N_subclasses
    Inst.Cellclass=[Inst.Cellclass ones(1,Inst.N_cells_vector(cl))*cl];
end
Inst.Connection_matrix=zeros(Inst.N_cells_total,Inst.N_cells_nonmod);

Inst.Cell_inclass=[];
for cl=1:size(Inst.N_cells_vector,2)
    Inst.Cell_inclass=[Inst.Cell_inclass 1:Inst.N_cells_vector(cl)];
end
% Making the templates for the cell parameters for the instance, allocation
% of the distribution-based values

% The standard fro the synaptic parameters:
% For synapses between the ordered pair [Post Pre], the properties (Power,
% delay, plasticity) may depend on both pre and postsynaptic neuron. 
% For Power and delay, the final values for each synapse are a sum of
% values from postsynaptic and presynaptic component. For plasticity,
% presynaptic and postsynaptic types of plasticity form a common list of
% plasticity patterns.
% For the Instance generation, presynaptic and postsynaptic components of
% the Power and delay will be recorded into a cell array where {Post,Pre} 
% coordinates refer to the synapse between the ordered pair, 
% {1/2} refere ro the postsynaptic/presynaptic input, 
% and the numerical vectors within the cell array represent the impact of each neuron. 
% After the connections are established, the properties of each individual
% synapse will be calculated.

% The plasticity simply go into the common lists. 
% To do it more comprehesive and understandable, we need a template with
% the placeholders: not only the numerica featuires, but the nuances of the
% structure of the instance depend on the genotype, so proper features
% allocation may seem very complicated. It is safer to do it more
% explicitly even if it means worse performance. Code optimisation can be
% done in the future versions of the routine. 



% The features of the neurons (not synapses). At this point I use solely
% the continuous values, the adjustment to integers will happen after the
% expression patterns allocation (if applicable)


%Inst.Maxpsp=0; % The maximal length of the psp in the entire circuit
for SubCl_post=1:Inst.N_subclasses
    
    for SubCl_pre=1:Inst.N_nonmod_subclasses
        %Inst.BasicPowers{SubCl_pre,SubCl_post}(1:Inst.N_cells_vector(SubCl_post))=trunc_cont_norm_dist(Inst.N_cells_vector(SubCl_post),Net2.BasicPowersPost(SubCl_post,SubCl_pre,1),Net2.BasicPowersPost(SubCl_post,SubCl_pre,2),NaN);
        %Inst.Delays{SubCl_pre,SubCl_post}(1:Inst.N_cells_vector(SubCl_post))=trunc_cont_norm_dist(Inst.N_cells_vector(SubCl_post),Net2.DelaysPost(SubCl_post,SubCl_pre,1),Net2.DelaysPost(SubCl_post,SubCl_pre,2),0);
        
        Inst.Connection_patternsPost{SubCl_post,SubCl_pre}(1:Inst.N_cells_vector(SubCl_post))=0; % These are not gaussian parameters, these are placeholders for the expression patterns realization
        Inst.Connection_patternsPre{SubCl_pre,SubCl_post}(1:Inst.N_cells_vector(SubCl_pre))=0;
        Inst.Affinity_patternsPost{SubCl_post,SubCl_pre}(1:Inst.N_cells_vector(SubCl_post))=0;
        Inst.Affinity_patternsPre{SubCl_pre,SubCl_post}(1:Inst.N_cells_vector(SubCl_pre))=0;
        

    end
    Inst.Thresholds{SubCl_post}(1:Inst.N_cells_vector(SubCl_post))=trunc_cont_norm_dist(Inst.N_cells_vector(SubCl_post),Net2.Thresholds(1,SubCl_post),Net2.Thresholds(2,SubCl_post),NaN);
    Inst.ThreshNoise{SubCl_post}(1:Inst.N_cells_vector(SubCl_post))=trunc_cont_norm_dist(Inst.N_cells_vector(SubCl_post),Net2.ThreshNoise(1,SubCl_post),Net2.ThreshNoise(2,SubCl_post),NaN);
    Inst.AbsRefract{SubCl_post}(1:Inst.N_cells_vector(SubCl_post))=trunc_cont_norm_dist(Inst.N_cells_vector(SubCl_post),Net2.AbsRefract(1,SubCl_post),Net2.AbsRefract(2,SubCl_post),0);
    Inst.RecurCon{SubCl_post}(1:Inst.N_cells_vector(SubCl_post))=trunc_cont_norm_dist(Inst.N_cells_vector(SubCl_post),Net2.RecurCon(1,SubCl_post),Net2.RecurCon(2,SubCl_post),0);
    Inst.RecurCon{SubCl_post}(Inst.RecurCon{SubCl_post}<-1)=-1;
    Inst.RecurCon{SubCl_post}(Inst.RecurCon{SubCl_post}>1)=1; % rapid brushing, required here for proper connections building
    Inst.RelRefract{SubCl_post}=[];
    % Place for the threshold variation 
    Inst.RelRefract_PlastID{SubCl_post}=[]; % placeholder for the relative refracterity patterns list
    %
    
end

% The hidden parameters drafting and allocation. After this stage, it
% becomes immutable. Also, at this stage, the expression patterns
% allocation to a neurons' features is done, but not to synapse features,
% so the routine will be repeated lated without the hidden values drafting.
if isfield(Net2,'Expression_patterns')
for Pat=1:size(Net2.Expression_patterns,2)
    % Net2.Expression_patterns{1,Pat}{1,1} - the ID of the current subclass
    if Inst.N_cells_vector(Net2.Expression_patterns{1,Pat}{1,1})==1
        Inst.Exp{Pat}=0.5;
    else
        Inst.Exp{Pat}=0:1/(Inst.N_cells_vector(Net2.Expression_patterns{1,Pat}{1,1})-1):1; % The hidden parameter draft for the pattern
    end
    
    order=randperm(size(Inst.Exp{Pat},2));
    Inst.Exp{Pat}=Inst.Exp{Pat}(order); % The randomization of the hidden parameter allocation
    Savedpat_in=1; % Saving of the expressions for the input and output neurons for the case if selective allocation will be necessary. 
    if ismember(Net2.Expression_patterns{Pat}{1},1:size(Net2.Cells.Input,2))
        Inst.In_order{Savedpat_in}{1}=Pat;
        Inst.In_order{Savedpat_in}{2}=order;
        Savedpat_in=Savedpat_in+1;
    end
    Savedpat_out=1;
    if ismember(Net2.Expression_patterns{Pat}{1},size(Net2.Cells.Input,2)+1:size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2))
        Inst.In_order{Savedpat_out}{1}=Pat;
        Inst.In_order{Savedpat_out}{2}=order;
        Savedpat_out=Savedpat_out+1;
    end
    
    
    
    for effect=2:size(Net2.Expression_patterns{1, Pat},2)
        %Pattern_values=Net2.Expression_patterns{1,Pat}{1,effect}{1,2}(1)+Net2.Expression_patterns{1,Pat}{1,effect}{1,2}(2).*Inst.Exp{Pat}+Net2.Expression_patterns{1,Pat}{1,effect}{1,2}(3).*Inst.Exp{Pat}.^2;
        Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
        if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Post connection pattern
            %Inst.Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+cont2int(Pattern_values);
            Inst.Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Pre expression pattern
            Inst.Connection_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Connection_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3) % Post affinity pattern  
            Inst.Affinity_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Affinity_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(1:Inst.N_cells_vector(Net2.Expression_patterns{1,Pat}{1,1}))+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4) % Pre affinity pattern  
            Inst.Affinity_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Affinity_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+Pattern_values;
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Basic Power, postsynaptic (into the current class)
%             Inst.BasicPowers{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}{1}=Inst.BasicPowers{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}{1}+Pattern_values;
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Basic Power, presynaptic (from the current class)
%             Inst.BasicPowers{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}{2}=Inst.BasicPowers{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}{2}+Pattern_values;
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Delay, postsynaptic (into the current class)
%             Inst.Delays{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}{1}=Inst.Delays{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}{1}+cont2int(Pattern_values);
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Delay, presynaptic (from the current class)
%             Inst.Delays{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}{2}=Inst.Delays{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}{2}+cont2int(Pattern_values);
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==4) % Threshold of activation 
            Inst.Thresholds{Net2.Expression_patterns{1,Pat}{1,1}}=Inst.Thresholds{Net2.Expression_patterns{1,Pat}{1,1}}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==5) % Absolute refracterity
            Inst.AbsRefract{Net2.Expression_patterns{1,Pat}{1,1}}=Inst.AbsRefract{Net2.Expression_patterns{1,Pat}{1,1}}+Pattern_values;
        else
            
            %warning(['Unknown adressation of the expression pattern, code ' num2str(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)) ' ' num2str(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2))])
        end
    end
end
end


% Filling in the connection matrix

for post=1:Inst.N_subclasses
    for pre=1:Inst.N_nonmod_subclasses
        
        Inst.Connection_matrix(Inst.N_cells_vector_incremental_pointer(post):Inst.N_cells_vector_incremental_pointer(post+1)-1,Inst.N_cells_vector_incremental_pointer(pre):Inst.N_cells_vector_incremental_pointer(pre+1)-1)=make_connections_v2([Inst.N_cells_vector(post) Net2.Connections(post,pre,1) Net2.Connections(post,pre,2) Inst.N_cells_vector(pre) Net2.Connections(post,pre,3)],Inst.Connection_patternsPost(post,pre),Inst.Connection_patternsPre(pre,post),{},{});
        % Here is the place for the individual recurrent connection exception to be introduced
        
        if pre==post % If the connection is recurrent on the subclass level
            for neu=1:Inst.N_cells_vector(post) % for each neuron of the recurrently connected subclass
                if Inst.RecurCon{post}(neu)>=0 
                    if Inst.Connection_matrix(Inst.N_cells_vector_incremental_pointer(post)-1+neu,Inst.N_cells_vector_incremental_pointer(post)-1+neu)==0 % if the neuron itself is not recurrently connected 
                        if rand<=Inst.RecurCon{post}(neu) % and if it has to be
                            Inst.Connection_matrix(Inst.N_cells_vector_incremental_pointer(post)-1+neu,Inst.N_cells_vector_incremental_pointer(post)-1+neu)=1; % then connect
                        end
                    end
                    
                elseif Inst.RecurCon{post}(neu)<0
                    if Inst.Connection_matrix(Inst.N_cells_vector_incremental_pointer(post)-1+neu,Inst.N_cells_vector_incremental_pointer(post)-1+neu)==1 % if the neuron itself is recurrently connected 
                        if rand<=-Inst.RecurCon{post}(neu) % And it should not be
                            Inst.Connection_matrix(Inst.N_cells_vector_incremental_pointer(post)-1+neu,Inst.N_cells_vector_incremental_pointer(post)-1+neu)=0; % disconnect
                        end
                    end
                    
                end
                
            end
        end
        
      
    end
end




% Making the list of synapses with the individual synapse parameters evaluated
synapse_counter=0;
Inst.actual_synapses_number(1:Inst.N_subclasses,1:Inst.N_subclasses)=0;
Inst.ConID=Inst.Connection_matrix;
for post_class=1:Inst.N_subclasses
    for pre_class=1:Inst.N_nonmod_subclasses
        Inst.Synapses.Sources{post_class,pre_class}=[]; % Initialisation
        Cl_synapse_counter=0;
        for post=1:Inst.N_cells_vector(post_class)
            for pre=1:Inst.N_cells_vector(pre_class)
                if Inst.Connection_matrix(Inst.N_cells_vector_incremental_pointer(post_class)+post-1,Inst.N_cells_vector_incremental_pointer(pre_class)+pre-1)==1
                    synapse_counter=synapse_counter+1;
                    Inst.ConID(Inst.N_cells_vector_incremental_pointer(post_class)+post-1,Inst.N_cells_vector_incremental_pointer(pre_class)+pre-1)=synapse_counter;
                    Cl_synapse_counter=Cl_synapse_counter+1;
                    Inst.Synapses.IDs(synapse_counter,1)=synapse_counter; % synapse ID
                    Inst.Synapses.IDs(synapse_counter,2:3)=[post_class pre_class]; % subclasses IDs
                    Inst.Synapses.IDs(synapse_counter,4:5)=[Inst.N_cells_vector_incremental_pointer(post_class)+post-1 Inst.N_cells_vector_incremental_pointer(pre_class)+pre-1]; % neurons IDs
                    Inst.Synapses.IDs(synapse_counter,6)=Cl_synapse_counter; % ID of the synapse in its class
                    Inst.Synapses.Sources{post_class,pre_class}(Cl_synapse_counter,1)=post; % Postsynaptic class neurons IDs (for pattern indexation)
                    Inst.Synapses.Sources{post_class,pre_class}(Cl_synapse_counter,2)=pre; % Presynaptic class neurons IDs (for pattern indexation)
%                     Inst.Synapses.Powers(synapse_counter,1)=Inst.BasicPowers{post_class,pre_class}{1}(post)*abs(Inst.BasicPowers{post_class,pre_class}{2}(pre));
%                     Inst.Synapses.Delays(synapse_counter,1)=Inst.Delays{post_class,pre_class}{1}(post)+Inst.Delays{post_class,pre_class}{2}(pre);
% 
%                     Inst.Synapses.Shapes{synapse_counter}=Inst.PSPshapes{post_class,pre_class}{post};
                    Inst.actual_synapses_number(post_class,pre_class)=Inst.actual_synapses_number(post_class,pre_class)+1;
                    
                end
            end
        end
        Inst.Synapses_PlastID{post_class,pre_class}{3}=[]; % placeholder for the plasticity patterns reference table
    end
end



        

% Now, in the branch Inst.actual_synapses_number we have the actual number
% of synapses for each ordered pair of subclasses (presynaptic are in
% columns). The normal distributions for the synapse features are to be
% evaluated for these numbers, but they potentially can be affected by the
% expression patterns from both presynaptic and postsynaptic sides. 

% At this point, we can allocate the valies for the synaptic features:
% delay, power, plasticity patterns and psp shapes

for SubCl_post=1:Inst.N_subclasses
    for SubCl_pre=1:Inst.N_nonmod_subclasses
        Inst.Synapses.Delays{SubCl_post,SubCl_pre}=trunc_cont_norm_dist(Inst.actual_synapses_number(SubCl_post,SubCl_pre),Net2.Delays(SubCl_post,SubCl_pre,1),Net2.Delays(SubCl_post,SubCl_pre,2),0); % N-based delays for each synapse
        Inst.Synapses.Powers{SubCl_post,SubCl_pre}=trunc_cont_norm_dist(Inst.actual_synapses_number(SubCl_post,SubCl_pre),Net2.BasicPowers(SubCl_post,SubCl_pre,1),Net2.BasicPowers(SubCl_post,SubCl_pre,2),NaN); % N-based powers for each synapese
        Inst.Synapses.PSPparam{SubCl_post,SubCl_pre}=[trunc_cont_norm_dist(Inst.actual_synapses_number(SubCl_post,SubCl_pre),Net2.PSPshape{SubCl_post,SubCl_pre}(1,1),Net2.PSPshape{SubCl_post,SubCl_pre}(2,1),1)'...
            trunc_cont_norm_dist(Inst.actual_synapses_number(SubCl_post,SubCl_pre),Net2.PSPshape{SubCl_post,SubCl_pre}(1,2),Net2.PSPshape{SubCl_post,SubCl_pre}(2,2),1)' trunc_cont_norm_dist(Inst.actual_synapses_number(SubCl_post,SubCl_pre),Net2.PSPshape{SubCl_post,SubCl_pre}(1,3),Net2.PSPshape{SubCl_post,SubCl_pre}(2,3),1)']; % N-based PSP parameters for each synapse
       
        Inst.Synapses.Plast{SubCl_post,SubCl_pre}{1}=[]; % Placeholders for four types of the plasticity we have in the model. May be more. 
        Inst.Synapses.Plast{SubCl_post,SubCl_pre}{2}=[];
        Inst.Synapses.Plast{SubCl_post,SubCl_pre}{3}=[];
        
        for SubCl_mod=1:Inst.N_mod_subclasses
            %if Inst.actual_synapses_number(SubCl_post,SubCl_pre)>0
                Inst.Mod_Connection_patternsPost{SubCl_post,SubCl_pre,SubCl_mod}(1:Inst.actual_synapses_number(SubCl_post,SubCl_pre))=0; % These are not gaussian parameters, these are placeholders for the expression patterns realization
                Inst.Mod_Affinity_patternsPost{SubCl_post,SubCl_pre,SubCl_mod}(1:Inst.actual_synapses_number(SubCl_post,SubCl_pre))=0;
            %end
            Inst.Mod_Connection_patternsPre{SubCl_post,SubCl_pre,SubCl_mod}(1:Inst.N_cells_vector(SubCl_mod+Inst.N_nonmod_subclasses))=0;
            Inst.Mod_Affinity_patternsPre{SubCl_post,SubCl_pre,SubCl_mod}(1:Inst.N_cells_vector(SubCl_mod+Inst.N_nonmod_subclasses))=0;
        
        end
    end
end





% The modulation connection matrices for the synapses.
% For each ordered pair of subclasses, where presynapse is ionotropic and postsynapse is anything, we have a number of synapses (defined previously
% by the make_connections_v2 routine). Now for each of these synapses and each neuron in the given modulatory class we define presence of the modulation
% connection (one synapse can receive multiple inputs from the same modulatory subclass but not from the same modulatory neuron). The process
% of the connection making works the same way as the connections making for the regular connecton matrix.

% First, application of the expression patterns for the modulatory
% connections. 
if isfield(Net2,'Expression_patterns')
for Pat=1:size(Inst.Exp,2)
    % Net2.Expression_patterns{1,Pat}{1,1} - the ID of the current subclass
    for effect=2:size(Net2.Expression_patterns{1, Pat},2)
        N_syn=NaN;
        % Check if there are any connections at all
        if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==8)&&(ismember(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2),[1 3])) % Modulatory , postsynaptic effect 
            N_syn=size( Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)} );
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==8)&&(ismember(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2),[2 4]))   % Modulatory , presynaptic effect 
            N_syn=size(Inst.Exp{Pat}); 
        end
        if N_syn>0
            
        if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==8)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Modulatory synapse connections number, postsynaptic effect 
            % In this case the connection pattern comes from one of the subclasses that are being modulated. The specific arrangement of the effect depends on the subclass: the expression of a neuron affects all its synapses, so each neuron of the expressing subclass "contributes" as many times as the number of synapses it makes
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the postsynaptic neuron in the ionotropic pair 
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            elseif Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the presynaptic neuron in the ionotropic pair 
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            else error('Wrong pattern allocation link for the modulatory connections number: connections, postsynaptic')
            end
            Inst.Mod_Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}=Inst.Mod_Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==8)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Modulatory synapse connections number, presynaptic effect 
            % In this case the pattern comes from the modulatory input, and every modulatory neuron contributes to every connection it makes
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            else error('Wrong pattern allocation link for the modulatory connections number: connections, presynaptic')    
            end
            Inst.Mod_Connection_patternsPre{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}=Inst.Mod_Connection_patternsPre{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==8)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3) % Modulatory synapse affinity, postsynaptic effect 
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the postsynaptic neuron in the ionotropic pair 
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            elseif Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the presynaptic neuron in the ionotropic pair 
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            else error('Wrong pattern allocation link for the modulatory connections number: affinity, postsynaptic')
            end
            Inst.Mod_Affinity_patternsPost{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}=Inst.Mod_Affinity_patternsPost{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==8)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4) % Modulatory synapse affinity, presynaptic effect     
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            else error('Wrong pattern allocation link for the modulatory connections number: affinity, presynaptic')    
            end
            Inst.Mod_Affinity_patternsPre{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}=Inst.Mod_Affinity_patternsPre{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}+Pattern_values;
        end
        end
    end
end
end

% Allocation of the gaussian features to the modulatory connections. 
%1. Pre-creation of empty tables for the lists of targets

for SubCl_mod=1:Inst.N_mod_subclasses % For each modulatory subclass
    for mod_cell=1:Inst.N_cells_vector(Inst.Mod_pointer(SubCl_mod))
        Inst.Synapses.Mod_targets{SubCl_mod}{mod_cell}=[];
    end
end
    
    

for SubCl_post=1:Inst.N_subclasses
    for SubCl_pre=1:Inst.N_nonmod_subclasses
        for SubCl_mod=1:Inst.N_mod_subclasses
            Inst.Mod_Connection_matrix{SubCl_post,SubCl_pre,SubCl_mod}=make_connections_v2([size(Inst.Synapses.Delays{SubCl_post,SubCl_pre},2) Net2.Connections_Mod(SubCl_post,SubCl_pre,SubCl_mod,1) Net2.Connections_Mod(SubCl_post,SubCl_pre,SubCl_mod,2) Inst.N_cells_vector(SubCl_mod+Inst.N_nonmod_subclasses) Net2.Connections_Mod(SubCl_post,SubCl_pre,SubCl_mod,3)],Inst.Mod_Connection_patternsPost(SubCl_post,SubCl_pre,SubCl_mod),Inst.Mod_Connection_patternsPre(SubCl_post,SubCl_pre,SubCl_mod),Inst.Mod_Affinity_patternsPost(SubCl_post,SubCl_pre,SubCl_mod),Inst.Mod_Affinity_patternsPre(SubCl_post,SubCl_pre,SubCl_mod));
            % The modulation connections counter and assigning
            Mod_inputs_total_number=sum(sum(Inst.Mod_Connection_matrix{SubCl_post,SubCl_pre,SubCl_mod}));
            Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(1,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(1,1),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(2,1),NaN);    % A single trigger effect
            Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(1,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(1,1),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(2,1),NaN);    % A maximal number of triggers
            Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,1),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,1),NaN);    % Delay
            Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,2),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,2),NaN);    % In
            Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(3,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,3),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,3),NaN);    % Plateau
            Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(4,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,4),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,4),NaN);    % Out
            Mod_input_counter_global=0;
            mods_counter=zeros(1,size(Inst.Mod_Connection_matrix{SubCl_post,SubCl_pre,SubCl_mod},2));
            target_ID=1;
            Inst.Synapses.Mod_source{SubCl_post,SubCl_pre,SubCl_mod}=[];
            %Inst.Synapses.Mod_targets{SubCl_mod}={};
            for syn=1:size(Inst.Mod_Connection_matrix{SubCl_post,SubCl_pre,SubCl_mod},1) % For each post<-pre synapse possibly receiving a modulatory input from the subclass SubCl_mod
                %Mod_input_counter_local=0;
                
                for mod_cell=1:size(Inst.Mod_Connection_matrix{SubCl_post,SubCl_pre,SubCl_mod},2) % And for each possible modulatory input source of the subclass SubCl_mod
%                     if ~isfield(Inst.Synapses,'Mod_targets')
%                         Inst.Synapses.Mod_targets{SubCl_mod}{mod_cell}=[];
%                     end
                    if Inst.Mod_Connection_matrix{SubCl_post,SubCl_pre,SubCl_mod}(syn,mod_cell)==1 % if there is a modulatory input
                        Mod_input_counter_global=Mod_input_counter_global+1;
                        mods_counter(mod_cell)=mods_counter(mod_cell)+1;
                        %Mod_input_counter_local=Mod_input_counter_local+1;
                        % CHECK WHY Inst.Modsyn.Features and Inst.Synapses.Mod ARE IDENTICAL. IS THERE SOME IDEA HERE???
%                         Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(1,Mod_input_counter_global)=Inst.Modsyn.Features{SubCl_post,SubCl_pre,SubCl_mod}{1}(1,Mod_input_counter_global); % A single trigger effect
%                         Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(1,Mod_input_counter_global)=Inst.Modsyn.Features{SubCl_post,SubCl_pre,SubCl_mod}{2}(1,Mod_input_counter_global); % A maximal number of triggers
%                         Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,Mod_input_counter_global)=Inst.Modsyn.Features{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,Mod_input_counter_global); % Delay
%                         Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,Mod_input_counter_global)=Inst.Modsyn.Features{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,Mod_input_counter_global); % In
%                         Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(3,Mod_input_counter_global)=Inst.Modsyn.Features{SubCl_post,SubCl_pre,SubCl_mod}{3}(3,Mod_input_counter_global); % Plateau
%                         Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(4,Mod_input_counter_global)=Inst.Modsyn.Features{SubCl_post,SubCl_pre,SubCl_mod}{3}(4,Mod_input_counter_global); % Out
                        if (Mod_input_counter_global>1)&&(  ~isequal([Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,1) Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,2)],[Inst.Synapses.Mod_source{SubCl_post,SubCl_pre,SubCl_mod}(Mod_input_counter_global-1,1:2)]))
                            target_ID=target_ID+1;
                        end

                        Inst.Synapses.Mod_source{SubCl_post,SubCl_pre,SubCl_mod}(Mod_input_counter_global,:)=[Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,1) Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,2) mod_cell Mod_input_counter_global target_ID]; % A triplet of neurons (intra-subclass indexation) forming the modulated synapse
                        % Format: [Intraclass_ID_of_post Intraclass_ID_of_pre Intraclass_ID_of_mod Global_ID_of_mod_contact Local_ID of ratget synapse]
                        %mods_counter
                        %Inst.Synapses.Mod_targets{SubCl_mod}{mod_cell}(mods_counter(mod_cell),:)=[SubCl_post Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,1) SubCl_pre Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,2)];
                        Inst.Synapses.Mod_targets{SubCl_mod}{mod_cell}(size(Inst.Synapses.Mod_targets{SubCl_mod}{mod_cell},1)+1,:)=[SubCl_post Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,1) SubCl_pre Inst.Synapses.Sources{SubCl_post,SubCl_pre}(syn,2)];
%                     else
%                         
%                         Inst.Synapses.Mod_targets{SubCl_mod}{mod_cell}=[];
                    end
                end
            end
        end
    end 
end





% The plasticity parameters for the synapses

% 
if isfield(Net2,'Plasticity')
for plast=1:size(Net2.Plasticity,2)
    if Net2.Plasticity{plast}{1}==1000  % relative refracterity
        current_plN=size(Inst.RelRefract{Net2.Plasticity{plast}{3}},1); % The current number of relative refracterities in the current subclass of neurons
        Inst.RelRefract_PlastID{Net2.Plasticity{plast}{3}(1)}=[Inst.RelRefract_PlastID{Net2.Plasticity{plast}{3}(1)} plast];
%         Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{1}(1,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{2}(1,1),Net2.Plasticity{plast}{2}(2,1),0); % The beginning of the trigger observation interval
%         Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{1}(2,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{2}(1,2),Net2.Plasticity{plast}{2}(2,2),0); % The end of the trigger observation interval
        Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{2}(1,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{4}{1}(1),Net2.Plasticity{plast}{4}{1}(2),NaN); % The single trigger effect
        Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{3}(1,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{4}{2}(1),Net2.Plasticity{plast}{4}{2}(2),NaN); % Maximal number of simultaneous triggers
        Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{4}(1,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{4}{3}(1,1),Net2.Plasticity{plast}{4}{3}(2,1),NaN); % Trigger delay
        Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{4}(2,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{4}{3}(1,2),Net2.Plasticity{plast}{4}{3}(2,2),NaN); % In
        Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{4}(3,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{4}{3}(1,3),Net2.Plasticity{plast}{4}{3}(2,3),NaN); % Plato
        Inst.RelRefract{Net2.Plasticity{plast}{3}(1)}{current_plN+1}{4}(4,:)=trunc_cont_norm_dist(Inst.N_cells_vector(Net2.Plasticity{plast}{3}),Net2.Plasticity{plast}{4}{3}(1,4),Net2.Plasticity{plast}{4}{3}(2,4),NaN); % Out
        
    elseif Net2.Plasticity{plast}{1}==1 % pair-pulse facilitation/inhibition
        current_plN=size(Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1},1);
        Inst.Synapses_PlastID{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}=[Inst.Synapses_PlastID{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1} plast];
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}{current_plN+1}{2}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{1}(1),Net2.Plasticity{plast}{4}{1}(2),NaN); % The single trigger effect
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}{current_plN+1}{3}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{2}(1),Net2.Plasticity{plast}{4}{2}(2),NaN); % Maximal number of simultaneous triggers
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}{current_plN+1}{4}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,1),Net2.Plasticity{plast}{4}{3}(2,1),NaN); % Trigger delay
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}{current_plN+1}{4}(2,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,2),Net2.Plasticity{plast}{4}{3}(2,2),NaN); % In
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}{current_plN+1}{4}(3,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,3),Net2.Plasticity{plast}{4}{3}(2,3),NaN); % Plato
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{1}{current_plN+1}{4}(4,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,4),Net2.Plasticity{plast}{4}{3}(2,4),NaN); % Out

    elseif Net2.Plasticity{plast}{1}==2 % spike-timing plasticity
        current_plN=size(Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2},1);
        Inst.Synapses_PlastID{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}=[Inst.Synapses_PlastID{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2} plast];
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{1}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{2}(1,1),Net2.Plasticity{plast}{2}(2,1),0); % The beginning of the trigger observation interval
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{1}(2,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{2}(1,2),Net2.Plasticity{plast}{2}(2,2),0); % The end of the trigger observation interval
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{2}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{1}(1),Net2.Plasticity{plast}{4}{1}(2),NaN); % The single trigger effect
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{3}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{2}(1),Net2.Plasticity{plast}{4}{2}(2),NaN); % Maximal number of simultaneous triggers
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{4}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,1),Net2.Plasticity{plast}{4}{3}(2,1),NaN); % Trigger delay
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{4}(2,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,2),Net2.Plasticity{plast}{4}{3}(2,2),NaN); % In
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{4}(3,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,3),Net2.Plasticity{plast}{4}{3}(2,3),NaN); % Plato
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{2}{current_plN+1}{4}(4,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,4),Net2.Plasticity{plast}{4}{3}(2,4),NaN); % Out
        
    elseif Net2.Plasticity{plast}{1}==3 % hebbian plasticity
        current_plN=size(Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3},1); % The current number of the plasticity patterns of the given type in this synapse 
        Inst.Synapses_PlastID{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}=[Inst.Synapses_PlastID{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3} plast];
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{1}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{2}(1,1),Net2.Plasticity{plast}{2}(2,1),0); % The beginning of the trigger observation interval
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{1}(2,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{2}(1,2),Net2.Plasticity{plast}{2}(2,2),0); % The end of the trigger observation interval
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{2}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{1}(1),Net2.Plasticity{plast}{4}{1}(2),NaN); % The single trigger effect
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{3}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{2}(1),Net2.Plasticity{plast}{4}{2}(2),NaN); % Maximal number of simultaneous triggers
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{4}(1,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,1),Net2.Plasticity{plast}{4}{3}(2,1),NaN); % Trigger delay
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{4}(2,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,2),Net2.Plasticity{plast}{4}{3}(2,2),NaN); % In
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{4}(3,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,3),Net2.Plasticity{plast}{4}{3}(2,3),NaN); % Plato
        Inst.Synapses.Plast{Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)}{3}{current_plN+1}{4}(4,:)=trunc_cont_norm_dist(Inst.actual_synapses_number(Net2.Plasticity{plast}{3}(1),Net2.Plasticity{plast}{3}(2)),Net2.Plasticity{plast}{4}{3}(1,4),Net2.Plasticity{plast}{4}{3}(2,4),NaN); % Out

    end

end
end











% The expression pattern allocation for the synaptic features
if isfield(Net2,'Expression_patterns')
for Pat=1:size(Inst.Exp,2)
    % Net2.Expression_patterns{1,Pat}{1,1} - the ID of the current subclass
    for effect=2:size(Net2.Expression_patterns{1, Pat},2)
        %Pattern_values=Net2.Expression_patterns{1,Pat}{1,effect}{1,2}(1)+Net2.Expression_patterns{1,Pat}{1,effect}{1,2}(2).*Inst.Exp{Pat}+Net2.Expression_patterns{1,Pat}{1,effect}{1,2}(3).*Inst.Exp{Pat}.^2;
        %Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
%         if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Post connection pattern
%             Inst.Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Connection_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+cont2int(Pattern_values);
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Pre expression pattern
%             Inst.Connection_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Connection_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+cont2int(Pattern_values);
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3) % Post affinity pattern  
%             Inst.Affinity_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Affinity_patternsPost{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(1:Inst.N_cells_vector(Net2.Expression_patterns{1,Pat}{1,1}))+Pattern_values;
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4) % Pre affinity pattern  
%             Inst.Affinity_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Affinity_patternsPre{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+Pattern_values;
        N_syn=NaN;

        if (ismember(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1),[2 3 6 7]))
            %(~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}))
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Basic Power, postsynaptic (into the current class)
                if (Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)},1);
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Basic Power, presynaptic (from the current class)  
                if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}},1);
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Delay, postsynaptic (into the current class)
                if (Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)},1);
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Delay, presynaptic (from the current class)  
                if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}},1);
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1)&&(ismember(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1},[1 2 3])) % plasticity, postsynaptic 
                if (Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)},1);
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2)&&(ismember(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1},[1 2 3])) % plasticity, presynaptic  
                if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}},1);
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==7)
                if (Net2.Expression_patterns{1,Pat}{1,1}<=size(Inst.Synapses.Sources,1))&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)<=size(Inst.Synapses.Sources,2))
                    N_syn=size(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)},1);
                end
            end
            %Тут нужна проверка наличия синапсов с учётом пре и постсинаптической адресации. Всё вываливается из-за пустого массива. Проверка выше годится только для постсинаптических свойств
                
            
        if N_syn>0
            
        if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Basic Power, postsynaptic (into the current class)
            %Post_exp_values=Inst.Exp{Pat}([Inst.Synapses.Sources{post_class,pre_class}(:,1)])
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});

%             [Net2.Expression_patterns{1,Pat}{1,1} Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)]
%             Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}
%             Inst.Exp{Pat}
%             Pattern_values
            Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+Pattern_values;
%             Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}

        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Basic Power, presynaptic (from the current class)
            %Pre_exp_values=Inst.Exp{Pat}([Inst.Synapses.Sources{post_class,pre_class}(:,2)])
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});

%             [Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3) Net2.Expression_patterns{1,Pat}{1,1}]
%             Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}
%             Inst.Exp{Pat}
%             Pattern_values
            Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}=Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}+Pattern_values;
%             Inst.Synapses.Powers{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}

        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % Delay, postsynaptic (into the current class)
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            Inst.Synapses.Delays{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}=Inst.Synapses.Delays{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}+Pattern_values;
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % Delay, presynaptic (from the current class)
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            Inst.Synapses.Delays{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}=Inst.Synapses.Delays{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,1}}+Pattern_values;
            
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==4) % Threshold of activation 
%             Inst.Thresholds{Net2.Expression_patterns{1,Pat}{1,1}}=Inst.Thresholds{Net2.Expression_patterns{1,Pat}{1,1}}+Pattern_values;
%         elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==5) % Absolute refracterity
%             Inst.AbsRefract{Net2.Expression_patterns{1,Pat}{1,1}}=Inst.AbsRefract{Net2.Expression_patterns{1,Pat}{1,1}}+cont2int(Pattern_values);


        %%%% The expression patterns applied to the plasticity
        %%%%%%%%%%% SINGLE TRIGGER EFFECT %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, single trigger effect, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{2}(1,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{2}(1,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {2}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{2}(1,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{2}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, single trigger effect, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {2}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{2}(1,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{2}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
            
            
        %%%%%%%%%%% MAXIMAL NUMBER OF EFFECTS %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{3}(1,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{3}(1,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {3}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{3}(1,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{3}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {3}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{3}(1,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{3}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        
        
        %%%%%%%%%%% THE BEGINNING OF OBSERVATION INTERVAL %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==7)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{1}(1,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{1}(1,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)
                    error('The observation interval for the plasticity of type 1 is not defined')
                end
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {1}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(1,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==7)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)
                    error('The observation interval for the plasticity of type 1 is not defined')
                end
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {1}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(1,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
            
            
        %%%%%%%%%%% THE END OF OBSERVATION INTERVAL %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==8)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{1}(2,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{1}(2,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)
                    error('The observation interval for the plasticity of type 1 is not defined')
                end
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {1}(2,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(2,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(2,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==8)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)
                    error('The observation interval for the plasticity of type 1 is not defined')
                end
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {1}(2,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(2,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{1}(2,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
            
            
        %%%%%%%%%%% THE INITIAL DELAY OF PLASTICITY %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(1,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(1,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(1,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(1,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(1,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(1,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        
        
        %%%%%%%%%%% THE ONSET OF PLASTICITY %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(2,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(2,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(2,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(2,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(2,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(2,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(2,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(2,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
            
            
        %%%%%%%%%%% THE PLATEAU OF PLASTICITY %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==5)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(3,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(3,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(3,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(3,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(3,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==5)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(3,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(3,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(3,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
            
            
        %%%%%%%%%%% THE OFFSET OF PLASTICITY %%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==1) % plasticity, maximal number of effects, postsynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                Pattern_values=pattern_allocation(Inst.Exp{Pat},Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(4,:)=Inst.RelRefract{Net2.Plasticity{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}}{3}(1)}{size(Pattern_values,2)+1}{4}(4,:)+Pattern_values;
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Post_cl,Pre_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Post_cl,Pre_cl}  =  {Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(4,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Expression_patterns{1,Pat}{1}, Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(4,:) = Inst.Synapses.Plast{Net2.Expression_patterns{1,Pat}{1},  Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(2)}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(4,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==2) % plasticity, maximal number of effects, presynaptic
            if (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1000) % relative refracterity
                error('Effect of the expression on the relative refracterity can be only from the postsynaptic neuron')
            elseif (Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==1)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==2)||(Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}{1}==3) % other type of plasticity
                Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1),  Net2.Expression_patterns{1,Pat}{1}}(:,2)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                %Inst.Synapses.Plast{Pre_cl,Post_cl}{Plast_type}{Specific plast pattern} [The parameter]
                % {Pre_cl,Post_cl}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}    
                % {Plast_type}  =  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}
                % {Specific plast pattern}  =  find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                % [The parameter]  =  {4}(4,:)
                Plast_reference=find(Inst.Synapses_PlastID{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1), Net2.Expression_patterns{1,Pat}{1} }  {Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}  ==  Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)  );
                Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(4,:) = Inst.Synapses.Plast{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{3}(1)  , Net2.Expression_patterns{1,Pat}{1}}{Net2.Plasticity{Net2.Expression_patterns{1,Pat}{1,effect}{1}(4)}{1}}{Plast_reference}{4}(4,:)+Pattern_values;
            else
                error('pattern allocation to an unknown type of plasticity, single trigger effect')
            end
            
       %%%%%%%%%%%%%%%%%%%%%%%%% THE PSP PARAMETERS %%%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==7)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1) % PSP parameters, in
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            Inst.Synapses.PSPparam{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)=Inst.Synapses.PSPparam{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)+Pattern_values';
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==7)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2) % PSP parameters, plateau
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            Inst.Synapses.PSPparam{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,2)=Inst.Synapses.PSPparam{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,2)+Pattern_values';
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==7)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3) % PSP parameters, out
            Pattern_values=pattern_allocation(Inst.Exp{Pat}([Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,1)]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
            Inst.Synapses.PSPparam{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,3)=Inst.Synapses.PSPparam{Net2.Expression_patterns{1,Pat}{1,1},Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)}(:,3)+Pattern_values';
            
        else
            
            %warning(['Unknown adressation of the expression pattern, code ' num2str(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)) ' ' num2str(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2))])
        end
        end
        end
    end
end
end



% The expression pattern allocation for the modulation. 


if isfield(Net2,'Expression_patterns')
for Pat=1:size(Inst.Exp,2)
    % Net2.Expression_patterns{1,Pat}{1,1} - the ID of the current subclass
    for effect=2:size(Net2.Expression_patterns{1, Pat},2)
            
        % [9 1 x y z 1/2] - A singular modulaton effect. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod). All for synapse x<-y
        % [9 2 x y z 1/2] - Max number of triggers. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
        % [9 3 x y z 1/2] - Delay. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
        % [9 4 x y z 1/2] - Onset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
        % [9 5 x y z 1/2] - Plateau duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
        % [9 6 x y z 1/2] - Offset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
        %Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(1,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(1,1),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{1}(2,1),NaN);    % A single trigger effect
        %Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(1,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(1,1),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{2}(2,1),NaN);    % A maximal number of triggers
        %Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,1),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,1),NaN);    % Delay
        %Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,2),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,2),NaN);    % In
        %Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(3,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,3),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,3),NaN);    % Plateau
        %Inst.Synapses.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(4,:)=trunc_cont_norm_dist(Mod_inputs_total_number,Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(1,4),Net2.Mod{SubCl_post,SubCl_pre,SubCl_mod}{3}(2,4),NaN);    % Out
        
        if (ismember(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1),[9]))&&(~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})) % If the link is synaptic there are any synapses at all within the ordered pair and
      
        
        %%%%%%%%%%%%%%%%%%%%%%%%% THE SINGLE EFFECT %%%%%%%%%%%%%%%%%%%%%
        if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==1) % modulation, a single effect, postsynaptic
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1, Pat}{1}) % postsynaptic effect, expression in ionotropic postsynapse
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1, Pat}{1})      % postsynaptic effect, expression in ionotropic presynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            else error('Wrong pattern allocation link for the modulatory feature: postsynaptic single effect')            
            end
            if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}) 
                Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{1}(1,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{1}(1,:)+Pattern_values([Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,5)']);
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==1)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==2) % modulation, a single effect, presynaptic
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([ Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,3)']),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                    Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{1}(1,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{1}(1,:)+Pattern_values;
                end
            else error('Wrong pattern allocation link for the modulatory feature: presynaptic single effect')    
            end
                    
        %%%%%%%%%%%%%%%%%%%%%%%%% MAXIMAL NUMBER OF TRIGGERS %%%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==1) % modulation, a maximal number of triggers, postsynaptic
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1, Pat}{1}) % postsynaptic effect, expression in ionotropic postsynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1, Pat}{1})      % postsynaptic effect, expression in ionotropic presynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            else error('Wrong pattern allocation link for the modulatory feature: postsynaptic single effect')            
            end
            if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}) 
                Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{2}(1,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{2}(1,:)+Pattern_values([Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,5)']);
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==2)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==2) % modulation, a maximal number of triggers, presynaptic
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([ Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,3)']),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                    Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{2}(1,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{2}(1,:)+Pattern_values;
                end
            else error('Wrong pattern allocation link for the modulatory feature: presynaptic single effect')    
            end
                    
        %%%%%%%%%%%%%%%%%%%%%%%%% DELAY OF MODULATION EFFECT %%%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==1) % modulation, delay, postsynaptic
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1, Pat}{1}) % postsynaptic effect, expression in ionotropic postsynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1, Pat}{1})      % postsynaptic effect, expression in ionotropic presynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            else error('Wrong pattern allocation link for the modulatory feature: postsynaptic single effect')            
            end
            if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}) 
                Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(1,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(1,:)+Pattern_values([Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,5)']);
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==3)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==2) % modulation, delay, presynaptic
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([ Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,3)']),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                    Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(1,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(1,:)+Pattern_values;
                end
            else error('Wrong pattern allocation link for the modulatory feature: presynaptic single effect')    
            end
            
            
                    
        %%%%%%%%%%%%%%%%%%%%%%%%% ONSET OF MODULATION EFFECT %%%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==1) % modulation, onset, postsynaptic
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1, Pat}{1}) % postsynaptic effect, expression in ionotropic postsynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1, Pat}{1})      % postsynaptic effect, expression in ionotropic presynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            else error('Wrong pattern allocation link for the modulatory feature: postsynaptic single effect')            
            end
            if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}) 
                Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(2,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(2,:)+Pattern_values([Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,5)']);
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==4)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==2) % modulation, onset, presynaptic
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([ Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,3)']),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                    Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(2,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(2,:)+Pattern_values;
                end
            else error('Wrong pattern allocation link for the modulatory feature: presynaptic single effect')    
            end
            %Pattern_values
            %Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(2,:)
                   
        %%%%%%%%%%%%%%%%%%%%%%%%% PLATEAU OF MODULATION EFFECT %%%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==5)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==1) % modulation, onset, postsynaptic
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1, Pat}{1}) % postsynaptic effect, expression in ionotropic postsynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1, Pat}{1})      % postsynaptic effect, expression in ionotropic presynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            else error('Wrong pattern allocation link for the modulatory feature: postsynaptic single effect')            
            end
            if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}) 
                Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(3,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(3,:)+Pattern_values([Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,5)']);
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==5)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==2) % modulation, onset, presynaptic
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([ Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,3)']),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                    Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(3,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(3,:)+Pattern_values;
                end
            else error('Wrong pattern allocation link for the modulatory feature: presynaptic single effect')    
            end
            
        %%%%%%%%%%%%%%%%%%%%%%%%% OFFSET OF MODULATION EFFECT %%%%%%%%%%%%%%%%%%%%%
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==1) % modulation, onset, postsynaptic
            if (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)==Net2.Expression_patterns{1, Pat}{1}) % postsynaptic effect, expression in ionotropic postsynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,1)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)==Net2.Expression_patterns{1, Pat}{1})      % postsynaptic effect, expression in ionotropic presynapse 
                if ~isempty(Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([  Inst.Synapses.Sources{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3)    ,    Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4)}(:,2)    ]),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                end
            else error('Wrong pattern allocation link for the modulatory feature: postsynaptic single effect')            
            end
            if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}) 
                Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(4,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(4,:)+Pattern_values([Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,5)']);
            end
        elseif (Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(1)==9)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(2)==6)&&(Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(6)==2) % modulation, onset, presynaptic
            if Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)==Net2.Expression_patterns{1,Pat}{1} % If the link refers to the modulatory neuron in the triplet
                if ~isempty(Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses})
                    Pattern_values=pattern_allocation(Inst.Exp{Pat}([ Inst.Synapses.Mod_source{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}(:,3)']),Net2.Expression_patterns{1,Pat}{1,effect}{1,2});
                    Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(4,:)=Inst.Synapses.Mod{Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(3),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(4),Net2.Expression_patterns{1,Pat}{1,effect}{1,1}(5)-Inst.N_nonmod_subclasses}{3}(4,:)+Pattern_values;

                end
            else error('Wrong pattern allocation link for the modulatory feature: presynaptic single effect')    
            end
                    

            
        end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Brushing: after the continuous values for the features are obtained from the normal distributions and the expression patterns, some (many) of them need to be turned into integers and truncated again if necessary. 

%reserve=Inst.AbsRefract{4};


% Neuron parameters
% Threshold: no brushing required
for sb_cl=1:size(Inst.AbsRefract,2) % Absolute refracterity
    Inst.AbsRefract{sb_cl}=cont2int(Inst.AbsRefract{sb_cl});
end
for sb_cl=1:size(Inst.RelRefract,2) % Relative refracterity
    for plast_pat=1:size(Inst.RelRefract{sb_cl},2)
%         Inst.RelRefract{sb_cl}{plast_pat}{1}(1,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{1}(1,:)); % Beginning of observation interval (CHECK AND REMOVE, RELATIVE REFRACTERITY DOES NOT HAVE OBSERVATION PERIOD)
%         Inst.RelRefract{sb_cl}{plast_pat}{1}(1,Inst.RelRefract{sb_cl}{plast_pat}{1}(1,:)<1)=1; % truncation
%         Inst.RelRefract{sb_cl}{plast_pat}{1}(2,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{1}(2,:)); % End of observation interval
%         Inst.RelRefract{sb_cl}{plast_pat}{1}(2,Inst.RelRefract{sb_cl}{plast_pat}{1}(2,:)<1)=1; % truncation
%         Inst.RelRefract{sb_cl}{plast_pat}{1}(2,Inst.RelRefract{sb_cl}{plast_pat}{1}(2,:)<Inst.RelRefract{sb_cl}{plast_pat}{1}(1,:))=Inst.RelRefract{sb_cl}{plast_pat}{1}(1,Inst.RelRefract{sb_cl}{plast_pat}{1}(2,:)<Inst.RelRefract{sb_cl}{plast_pat}{1}(1,:));  % CHECK!
        % Single trigger effect: no brushing required
        Inst.RelRefract{sb_cl}{plast_pat}{2}(1,Inst.RelRefract{sb_cl}{plast_pat}{2}(1,:)<-1)=-1; % Truncation: the single trigger effect cannot be less than -1
        Inst.RelRefract{sb_cl}{plast_pat}{3}(1,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{3}(1,:)); % maximal simultaneous number of triggers
        Inst.RelRefract{sb_cl}{plast_pat}{3}(1,Inst.RelRefract{sb_cl}{plast_pat}{3}(1,:)<1)=1; % truncation. Cannot be less than 1, otherwise no plasticity at all (Think, maybe 0 is a viable option: some symapes of a certain type have plasticity, some do not)
        Inst.RelRefract{sb_cl}{plast_pat}{4}(1,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{4}(1,:)); % delay
        Inst.RelRefract{sb_cl}{plast_pat}{4}(1,Inst.RelRefract{sb_cl}{plast_pat}{4}(1,:)<1)=1; % The "time" fearures of the plasticity are indices for the shapes list, they have to be positive integers. 
        Inst.RelRefract{sb_cl}{plast_pat}{4}(2,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{4}(2,:)); % onset
        Inst.RelRefract{sb_cl}{plast_pat}{4}(2,Inst.RelRefract{sb_cl}{plast_pat}{4}(2,:)<1)=1;
        Inst.RelRefract{sb_cl}{plast_pat}{4}(3,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{4}(3,:)); % plateau
        Inst.RelRefract{sb_cl}{plast_pat}{4}(3,Inst.RelRefract{sb_cl}{plast_pat}{4}(3,:)<1)=1;
        Inst.RelRefract{sb_cl}{plast_pat}{4}(4,:)=cont2int(Inst.RelRefract{sb_cl}{plast_pat}{4}(4,:)); % offset
        Inst.RelRefract{sb_cl}{plast_pat}{4}(4,Inst.RelRefract{sb_cl}{plast_pat}{4}(4,:)<1)=1;
    end
end

% Synapse parameters

for post_cl=1:size(Inst.Synapses.Delays,1) 
    for pre_cl=1:size(Inst.Synapses.Delays,2) 
        Inst.Synapses.Delays{post_cl,pre_cl}=cont2int(Inst.Synapses.Delays{post_cl,pre_cl}); % Delays
        Inst.Synapses.Delays{post_cl,pre_cl}(Inst.Synapses.Delays{post_cl,pre_cl}(:)<0)=0; % Truncation, delay cannot be less than zero
        % Powers are real-valued and do not have require truncation
        Inst.Synapses.PSPparam{post_cl,pre_cl}(:,1)=cont2int(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,1));
        Inst.Synapses.PSPparam{post_cl,pre_cl}(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,1)<1,1)=1;
        Inst.Synapses.PSPparam{post_cl,pre_cl}(:,2)=cont2int(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,2));
        Inst.Synapses.PSPparam{post_cl,pre_cl}(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,2)<1,2)=1;
        Inst.Synapses.PSPparam{post_cl,pre_cl}(:,3)=cont2int(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,3));
        Inst.Synapses.PSPparam{post_cl,pre_cl}(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,3)<1,3)=1;

        for Plst_type=1:size(Inst.Synapses.Plast{post_cl,pre_cl},2)
            for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type},2)
                if Plst_type~=1
                    Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,:)); % Beginning of the observation interval
                    Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,:)<1)=1; 
                    Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,:)); % End of the observation interval
                    Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,:)<1)=1; 
                    Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,:)<Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,:))=Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(2,:)<Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{1}(1,:));
                end
                % Single trigger effect: no brushing required
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{2}(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{2}(1,:)<-1)=-1; % Truncation: the single trigger effect cannot be less than -1
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{3}(1,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{3}(1,:)); % maximal simultaneous number of triggers
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{3}(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{3}(1,:)<1)=1; 

                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(1,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(1,:)); % delay
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(1,:)<1)=1; 
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,:)); % onset
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,:)<1)=1; 
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,:)); % plateau
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,:)<1)=1; 
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,:)=cont2int(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,:)); % offset
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,:)<1)=1; 


            end
        end
        % modulation
        for mod_cl=1:Inst.N_mod_subclasses % for each modulatory class
            for syn=1:size(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl},2)
                for mod_inp=1:size(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{syn},2)
                    %Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1}(1,:)=cont2int(Inst.Modsyn.Features{post_cl,pre_cl,mod_cl}{1}(1,:)); % A single trigger effect
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1}(1,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1}(1,:)<-1)=-1;
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{2}(1,:)=cont2int(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{2}(1,:)); % A maximal number of triggers
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{2}(1,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{2}(1,:)<1)=1;
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(1,:)=cont2int(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(1,:)); % Delay
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(1,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(1,:)<1)=1;
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,:)=cont2int(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,:)); % In
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,:)<1)=1;
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,:)=cont2int(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,:)); % Plateau
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,:)<1)=1;
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,:)=cont2int(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,:)); % Out
                    Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,:)<1)=1;

                    
                end
            end
        end
    end
end


% The PSP shapes
tech.Longest_shape=0;
% Runtime.Thrsh=[]; %TRANSFERED INTO A NEW FUNCTION, REMOVE AFTER TESTS
% Runtime.ThrshNoise=[];
% Runtime.AbsRefrRef=[];
for post_cl=1:size(Inst.Synapses.PSPparam,1) 
    for pre_cl=1:size(Inst.Synapses.PSPparam,2) 
        %Since the set of the possible PSP shapes is finite, the actual shapes parameters must be truncated from the top. It was not done during the instance construction because potentially the change of the file with PSP shapes is possible, so the upper limits will change
        Inst.Synapses.PSPparam{post_cl,pre_cl}(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,1)>size(PSP_db.PSP.in,2),1)=size(PSP_db.PSP.in,2); % Onset
        Inst.Synapses.PSPparam{post_cl,pre_cl}(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,2)>size(PSP_db.PSP.mid,2),2)=size(PSP_db.PSP.mid,2); % Plateau
        Inst.Synapses.PSPparam{post_cl,pre_cl}(Inst.Synapses.PSPparam{post_cl,pre_cl}(:,3)>size(PSP_db.PSP.out,2),3)=size(PSP_db.PSP.out,2); % Offset
        
        for syn=1:size(Inst.Synapses.PSPparam{post_cl,pre_cl},1)
            PSPlength=size(fliplr([PSP_db.PSP.in{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,1)} PSP_db.PSP.mid{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,2)} PSP_db.PSP.out{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,3)}]),2);
%             Runtime.PSP{post_cl,pre_cl}{syn}=fliplr([PSP_db.PSP.in{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,1)} PSP_db.PSP.mid{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,2)} PSP_db.PSP.out{Inst.Synapses.PSPparam{post_cl,pre_cl}(syn,3)}]);
%             Runtime.PSPlength{post_cl,pre_cl}{syn}=size(Runtime.PSP{post_cl,pre_cl}{syn},2);%TRANSFERED INTO A NEW FUNCTION, REMOVE AFTER TEST          
  
            if tech.Longest_shape<PSPlength+Inst.Synapses.Delays{post_cl,pre_cl}(syn)
                tech.Longest_shape=PSPlength+Inst.Synapses.Delays{post_cl,pre_cl}(syn);
            end
            for plast_type=1:3 % For every type of synaptic plasticity
                if ~isempty(Inst.Synapses.Plast{post_cl,pre_cl}{plast_type}) % If the contact between subclasses has this type of plasticity
                    for pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{plast_type},1) 
                        if ~isempty(Inst.Synapses.Plast{post_cl,pre_cl}{plast_type}{pat}) 
                            if ~isempty(Inst.Synapses.Plast{post_cl,pre_cl}{plast_type}{pat}{1}) % If it has the observation period
                            if tech.Longest_shape<=max(Inst.Synapses.Plast{post_cl,pre_cl}{plast_type}{pat}{1}(2,:)+Inst.Synapses.Delays{post_cl,pre_cl}) % If the observation period is longer than the current "blanc" raster
                                tech.Longest_shape=max(Inst.Synapses.Plast{post_cl,pre_cl}{plast_type}{pat}{1}(2,:)+Inst.Synapses.Delays{post_cl,pre_cl}); % Increase the "blanc" raster by it. 
                            end
                            end
                        end
                    end
                end
            end
        end
    end 
%     Runtime.Thrsh=[Runtime.Thrsh Inst.Thresholds{post_cl}]; % The thresholds in a long list   %TRANSFERED INTO A NEW FUNCTION, REMOVE AFTER TESTS
%     Runtime.ThrshNoise=[Runtime.ThrshNoise Inst.ThreshNoise{post_cl}]; % The threshold noises in a long list
%     Runtime.AbsRefrRef=[Runtime.AbsRefrRef Inst.AbsRefract{post_cl}]; % The absolute refracterity reference
end


for post_cl=1:size(Inst.Synapses.Plast,1) 
    for pre_cl=1:size(Inst.Synapses.Plast,2) 
        % The plasticity shapes and the track tables (Except Relative refracterity, type 1000)
        for Plst_type=1:size(Inst.Synapses.Plast{post_cl,pre_cl},2)
            for Plst_pat=1:size(Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type},2)
                % timescale truncation from the top
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,:)>size(Plast_db.plast.in,2))=size(Plast_db.plast.in,2); 
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,:)>size(Plast_db.plast.mid,2))=size(Plast_db.plast.mid,2); 
                Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,:)>size(Plast_db.plast.out,2))=size(Plast_db.plast.out,2); 
%                 for syn=1:size(Inst.Synapses.Powers{post_cl,pre_cl},2) % For each synapse between between an ordered pair of subclasses  %TRANSFERED INTO A NEW FUNCTION, REMOVE AFTER TESTS
%                     Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}=[zeros(1,Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(1,syn)) Plast_db.plast.in{Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(2,syn)} ...
%                         Plast_db.plast.mid{Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(3,syn)}  Plast_db.plast.out{Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{4}(4,syn)} ]; % A shape template (plateau=1)
%                     Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}=Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}*Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{2}(1,syn); % A final shape (amplitude considered)
%                     Runtime.Track{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn}(1:Inst.Synapses.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{3}(1,syn),1:size(Runtime.Plast{post_cl,pre_cl}{Plst_type}{Plst_pat}{syn},2))=0;
%                    
%                 end
            end 
        end
        % The modulation shapes and track tables
        for mod_cl=1:size(Inst.Synapses.Mod,3) 
            if ~isempty(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl})
                Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,:)>size(Plast_db.plast.in,2))=size(Plast_db.plast.in,2);
                Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,:)>size(Plast_db.plast.mid,2))=size(Plast_db.plast.mid,2);
                Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,:)>size(Plast_db.plast.out,2))=size(Plast_db.plast.out,2);
        
%                 for syn=1:size(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1},2) %TRANSFERED INTO A NEW FUNCTION, REMOVE AFTER TESTS
%                     Runtime.Mod{post_cl,pre_cl,mod_cl}{syn}=[zeros(1,Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(1,syn)) Plast_db.plast.in{Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(2,syn)}...
%                         Plast_db.plast.mid{Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(3,syn)} Plast_db.plast.out{Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{3}(4,syn)}];
%                     Runtime.Mod{post_cl,pre_cl,mod_cl}{syn}=Runtime.Mod{post_cl,pre_cl,mod_cl}{syn}*Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{1}(1,syn);
%                     Runtime.ModTrack{post_cl,pre_cl,mod_cl}{syn}(Inst.Synapses.Mod{post_cl,pre_cl,mod_cl}{2}(1,syn), 1:size(Runtime.Mod{post_cl,pre_cl,mod_cl}{syn},2))=0;
%                 end
            end
        end

    end 
    for Plst_pat=1:size(Inst.RelRefract{post_cl},2)
        Inst.RelRefract{post_cl}{Plst_pat}{4}(2,Inst.RelRefract{post_cl}{Plst_pat}{4}(2,:)>size(Plast_db.plast.in,2))=size(Plast_db.plast.in,2); 
        Inst.RelRefract{post_cl}{Plst_pat}{4}(3,Inst.RelRefract{post_cl}{Plst_pat}{4}(3,:)>size(Plast_db.plast.mid,2))=size(Plast_db.plast.mid,2); 
        Inst.RelRefract{post_cl}{Plst_pat}{4}(4,Inst.RelRefract{post_cl}{Plst_pat}{4}(4,:)>size(Plast_db.plast.out,2))=size(Plast_db.plast.out,2); 
%         for cell=1:Inst.N_cells_vector(post_cl) % For each cell of the postsynaptic subclass
%             Runtime.RelRef{post_cl}{Plst_pat}{cell}=[zeros(1,Inst.RelRefract{post_cl}{Plst_pat}{4}(1,cell)) Plast_db.plast.in{Inst.RelRefract{post_cl}{Plst_pat}{4}(2,cell)} ...
%                 Plast_db.plast.mid{Inst.RelRefract{post_cl}{Plst_pat}{4}(3,cell)}  Plast_db.plast.out{Inst.RelRefract{post_cl}{Plst_pat}{4}(4,cell)} ]; % A shape template (plateau=1)
%             Runtime.RelRef{post_cl}{Plst_pat}{cell}=Runtime.RelRef{post_cl}{Plst_pat}{cell}*Inst.RelRefract{post_cl}{Plst_pat}{2}(1,cell); % A final shape (amplitude considered)
%             Runtime.RefTrack{post_cl}{Plst_pat}{cell}(1:Inst.RelRefract{post_cl}{Plst_pat}{3}(1,cell),1:size(Runtime.RelRef{post_cl}{Plst_pat}{cell},2))=0;
%                    
%         end
    end
end

%Runtime.AbsRefCounter=zeros(1,Inst.N_cells_total);

Runtime=make_runtime(Inst,PSP_db,Plast_db);




end