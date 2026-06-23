

% The Testing routine functions script turned into a function

function out=full_evo_sim(run,Point_Mut)


%%%%%%%%%%%%%%%%% genetic code generation %%%%%%%%%%%%%%%%%%

% All the genetic parameters will be stored in the Net2 variable, which is
% the main product of the routine





% The composition of the circuit
Net2.Cells.Input=[[1 1 1];[0 0 0]]; % SD has to be zero for the reserved classes
Net2.Cells.Output=[[1];[0]]; % SD has to be zero for the reserved classes
Net2.Cells.Ion=[[15 10];[2 1]]; % The upper row - means, the lower row - sigmas 
Net2.Cells.Mod=[[3 3];[2 1]];

n_subclasses=size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2)+size(Net2.Cells.Ion,2)+size(Net2.Cells.Mod,2);
n_nonmod=size(Net2.Cells.Input,2)+size(Net2.Cells.Output,2)+size(Net2.Cells.Ion,2);
% The number of connections between the classes, postsynapses in rows, pretsynapses in columns, order: [input output ion mod]

Net2.Connections(:,:,1)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[1 1 1 1 2 3];[1 1 1 2 3 5];[1 1 1 1 2 3];[1 1 1 1 2 1.5];[1 1 1 1 2 1]]; % Mean
Net2.Connections(:,:,2)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 2 1];[0 0 1 2 1 1];[0 0 0 0 5 1];[0 0 0 0 0 1];[0 0 0 0 0 0]]; % sigma of the input
Net2.Connections(:,:,3)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 1 0];[0 0 0 1 1 1];[0 0 0 0 1 2];[0 0 0 0 0 2];[0 0 0 0 0 1]]; % sigma of the output
% Net2.Connections(:,:,2)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0]]; % sigma of the input
% Net2.Connections(:,:,3)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0]]; % sigma of the output

Net2.Connections_Mod(1:8,1:6,1:2,1)=3; % Mean  (post of target, pre of target, modulation source, mean number)
Net2.Connections_Mod(1:8,1:6,1:2,2)=1; % sigma of the input
Net2.Connections_Mod(1:8,1:6,1:2,3)=1; % sigma of the output


% The Basic Powers values for the modulatory effect must have a special format similar to the plasticity single trigger effect 
Net2.BasicPowers(:,:,1)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[1 1 1 1 1 -1];[3 3 3 1 1 -1];[1 1 1 1 2 -1];[1 1 1 1 5 -1];[1 1 1 1 1 -1]]; % Mean
Net2.BasicPowers(:,:,2)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0.1 1];[0 0 0 0 0.2 0];[0 0 0 0 0 0]]; % sigma

% The Delays for the modulatory effect potentially could be stored in the same format at the Connections, but it is easier to store them as the delay of the effect itself (similar to the plasticity)
Net2.Delays(:,:,1)=[[2 2 2 2 2 2];[2 2 2 2 2 2];[2 2 2 2 2 2];[2 2 2 2 2 2];[2 2 2 2 2 2];[2 2 2 2 2 2];[2 2 2 2 2 2];[2 2 2 2 2 2]]; % Mean
Net2.Delays(:,:,2)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 1 1];[0 0 0 0 0.5 0];[0 0 0 0 0.1 0];[0 0 0 0 0 0]]; % sigma
%Net2.Delays(:,:,2)=[[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0]]; % sigma

for post=1:n_subclasses
    for pre=1:n_nonmod % The PSP is nonapplicable for the modulatore "presynapse"
        Net2.PSPshape{post,pre}=[[2 2 4];[0 0 0]]; % A template for the PSP shape, same for all synapses at this point. The second row is sigma
    end
end
Net2.PSPshape{3,4}=[[1 1 3];[0 0.1 1]]; % Couple of deviating PSP shapes to control the flow. 
Net2.PSPshape{5,1}=[[3 3 5];[1 1 1]];

% The parameters of the neurons themselves are applicable for the
% modulation
Net2.Thresholds=[[0 0 0 0 1 1 1 0];[0 0 0 0 0 0 0 0]]; % second row is sigma
Net2.AbsRefract=[[0 0 0 0 5 1 1 0];[0 0 0 0 1 0.2 0 0]]; % second row is sigma
Net2.RecurCon=[[0 0 0 0.1 -0.5 0 0 0];[0 0 0 0 0 0 0 0]]; % An individual recurrent connection affinity. Value [-1 1]. [0 1] = p(switch 0 individual recurrent connection to 1 if 0). [-1 0] = -p(switch 1 individual recurrent connection to 0). Inapplicable for modulatory neurons 
Net2.ThreshNoise=[[0 0 0 0.1 1 0.1 0.1 1];[0 0 0 0.1 0.1 0.1 0 0]]; % A temporal sigma of the gaussian noise of the activation threshold (the mean of the noise is always zero, the steady alternation of threshold is in Net2.Thresholds)

modpattern1={[0.4;0.1],[10;1],[[5 2 5 15];[1 0 1 1]]}; % A template just for tests. First - effect of a single trigger, second - max number of triggers, third - [delay, in, plateau, out] all with sigmas
modpattern2={[-0.2;0.1],[10;1],[[5 2 5 15];[1 0 1 1]]};

for post=1:n_subclasses
    for pre=1:n_nonmod 
        Net2.Mod{post,pre,1}=modpattern1;
        Net2.Mod{post,pre,2}=modpattern2;
    end
end


% The patterns and the coordinate system
% Code: each expression pattern refers to a subclass of neurons (a scalar
% in the first cell). Enumeration of the subclasses: [In Out Ion Mod]
% Other elements in the cells array (consisting of two vectors) contain the
% coordinate of the parameter (first vector) and the coefficients (second
% parameter).
% Coordinates:
% Class-specific targets
% [0 x] - The number of neurons of the class x (mostly applicable to the
% multiple instances case); TO BE ADDED LATER

% Neuron-specific targets
% [1 1 y] - Post_connection_pattern from the class y
% [1 2 y] - Pre_connection_pattern to the class y
% [1 3 y] - Post_affinity_pattern  from the class y
% [1 4 y] - Pre_affinity_pattern to the class y
% [4] - Thresholds pattern 
% [5] - AbsRefract pattern 
% [10] - RecurCon pattern
% [11] - ThreshNoise pattern

% Synapse-specific targets
% [2 1 y] - Basic_Powers pattern from the class y
% [2 2 y] - Basic_Powers pattern to the class y
% [3 1 y] - Delay pattern from the class y
% [3 2 y] - Delay pattern to the class y 
% [6 1 1/2 t] - For the plasticity pattern with ID t, a single act maximal effect. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 2 1/2 t] - For the plasticity pattern with ID t, Maximal simultaneous actions. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 3 1/2 t] - For the plasticity pattern with ID t, a delay. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 4 1/2 t] - For the plasticity pattern with ID t, an onset duration. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 5 1/2 t] - For the plasticity pattern with ID t, a plateau duration. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 6 1/2 t] - For the plasticity pattern with ID t, an offset duration. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 7 1/2 t] - For the plasticity pattern with ID t, a triggering condition frame beginning. Third place: 1 for postsynaptic, 2 for presynaptic
% [6 8 1/2 t] - For the plasticity pattern with ID t, a triggering condition frame end. Third place: 1 for postsynaptic, 2 for presynaptic
% [7 1 y] - For the PSP from the class y, an onset duration. 
% [7 2 y] - For the PSP from the class y, a plateau duration. 
% [7 3 y] - For the PSP from the class y, an offset duration. 

% Modulation targets 
% Modulation can be affected by the genes of all three participants: post, pre_ion and pre_mod
% [8 1 x y z] - Post_connection_pattern for modulation of the synapse x<-y and the modulatory input z 
% [8 2 x y z] - Pre_connection_pattern for modulation of the synapse x<-y and the modulatory input z 
% [8 3 x y z] - Post_affinity_pattern for modulation of the synapse x<-y and the modulatory input z 
% [8 4 x y z] - Pre_affinity_pattern for modulation of the synapse x<-y and the modulatory input z 
% [9 1 x y 1/2] - A singular modulaton effect. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod). All for synapse x<-y
% [9 2 x y 1/2] - Max number of triggers. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 3 x y 1/2] - Delay. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 4 x y 1/2] - Onset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 5 x y 1/2] - Plateau duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% [9 6 x y 1/2] - Offset duration. The fifth place: 1 - postsynapse (any subclass, either x or y, modulatory y is forbidden), 2 presynapse (mod)
% 
% Net2.Expression_patterns{1}={[5],{[1 1 5],[0 1 3]},{[1 2 6],[1 1 0]},{[2 2 6],[1 1 1]},{[2 1 6],[1 2 3]}};
% Net2.Expression_patterns{2}={[5],{[1 1 5],[0 1 2]},{[3 1 5],[3 -1 -2]},{[3 2 5],[3 -1 -2]},{[2 2 5],[3 -1 -2]}};
% Net2.Expression_patterns{3}={[6],{[2 2 6],[0 1 -3]},{[9 1 6 6 7 1],[1 5 5]},{[8 1 6 5 8],[3 1 1]}};
% Net2.Expression_patterns{4}={[4],{[3 1 3],[0 0 -10]},{[2 1 3],[0 5 0]},{[6 2 1 4],[1 2 3]}};
% Net2.Expression_patterns{5}={[4],{[6 1 1 4],[0 2 3]},{[6 3 1 4],[1 4 8]}}; % The pattern to test the plasticity routine (postsynaptic)
% Net2.Expression_patterns{6}={[4],{[6 1 2 6],[1 2 0]},{[2 2 4],[0 0 -10]},{[6 4 1 4],[1 4 8]}}; % The pattern to test the plasticity routine (presynaptic)
% Net2.Expression_patterns{7}={[7],{[9 1 6 6 7 2],[1 10 10]},{[9 2 8 6 7 2],[1 10 10]},{[8 2 8 6 7],[0 0 -10]},{[9 1 7 5 7 1],[3 3 3]}}; % To test the modulation
% 
% 
% %The plasticity pattern for the synapse is legal only if the presynaptic neuron is non-modulatory
% Net2.Plasticity{1}={1,[NaN],[5 6],{[-0.3;0.1],[10;1],[[0 1 10 5];[0 0 2 1]]}}; % The coordinate in the outer cell array is the pattern's ID (immutable, reassignment is possible if the pattern in eliminated by mutation)
% Net2.Plasticity{2}={1000,[NaN],[6],{[0.4;0.1],[10;1],[[2 3 20 10];[0 0 2 1]]}};
% Net2.Plasticity{3}={2,[3 30;1 1],[6 5],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{4}={3,[3 6;1 1],[6 6],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{5}={3,[3 6;1 1],[8 5],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{6}={2,[3 6;1 1],[6 6],{[0.4;0.1],[10;1],[[2 2 2 15];[1 0 2 1]]}};
% Net2.Plasticity{7}={1000,[NaN],[6],{[0.4;0.1],[10;1],[[2 3 20 10];[0 0 2 1]]}};
clear Net2
%load('MyNet2_I_1_10_O_1_NotEvolved.mat')
%load('MyNet_Evolved_I_1_1_1_O_1_probabilistic.mat');
%load('MyNet2_I_111_O_1111_NotEvolved.mat')


%for outeriter=4:4 % The envelope function to run the routine multiple times
%clear


Save_mode=2;
% Save_mode = 1 - save all the genotypes, phenotypes and rasters, leads to a very big export variables, but direct interaction with the results later on results
% Save_mode = 2 - save only the genotypes. Everything else can be resimulated using the seeds
% Save_mode = 3 - don't save even the genotypes, only the initial settings and the seeds. Smallest files, but require full resimulation of the evolution

load('MyNet2_I_111_O_1_NotEvolved.mat');
%Net2=sleep_circuit_build;
%load('MyNet2_I_111111_O_111_NotEvolved.mat');
%load('D:\Iteration4\New_model\MyNet2_I_111_O_11_NotEvolved.mat');
%load('C:\Users\Zavrd\Documents\MyNet2_I_111_O_1_diffOR_Evolved_1.mat');
%load('MyNet2_I_111_O_1111_NotEvolved.mat');
%load('MyNet2_I_11_O_11_NotEvolved.mat');
%load('MyNet2_I_1_10_O_1_NotEvolved.mat');
%load('MyNet2_I_1_20_O_1_Quant_Evolved_modified.mat')
Zoo.Settings.Core_Genotype=Net2;
Zoo.Settings.PreMutationsN=50;
PSP_db=load('PSP_shapes1.mat');
Plast_db=load('plast_shapes1.mat');
%load('point_mutation_template2.mat');
struct_mut=[0.2 0.025 0.2 0.05 0.2 0.05 0.2 0.05]; %[p_cl_div p_cl_rem p_plast_add p_plast_remove p_exp_add p_exp_remove p_link_add p_link_remove]


Gen_N=300; % Number of generations
Inst_N=20; % Number of instances per generation
for init=1:Inst_N
    for q=1:Zoo.Settings.PreMutationsN
        [Net2, ~]=mutation_routine(Net2,Point_Mut,struct_mut,NaN);
    end
    legality(Net2)
    Primogen(init)={Net2};
end
%Zoo.Ancestor=Net2;

% The input

% Input=rand(3,500);
% Threshold=0.8;
% Input(Input>=Threshold)=1;
% Input(Input<Threshold)=0;



Net2=rmfield(Net2,'Info'); % remove the Info field to not save the copy of it thousand times. It is applicable only for the initial genotype anyway. 

maxfit(1:Gen_N)=NaN;
Zoo.Settings.Save_mode=Save_mode; % The saving mode used to create an export file
Zoo.Settings.Ancestor=Net2;
Zoo.Settings.PSP=PSP_db;
Zoo.Settings.Plast=Plast_db;
Zoo.Settings.Point_mutation_template=Point_Mut;
Zoo.Settings.Struct_mutation_probs=struct_mut;
Zoo.Parents(1:Gen_N,1:Inst_N)=NaN;
Zoo.Mutation_seeds=randi(1000000,Gen_N,Inst_N); % Seeds for the RNG for the mutation
Zoo.Unfolding_seeds=randi(1000000,Gen_N,Inst_N); % Seeds for the RNG for the instance generation
Zoo.Simulation_seeds=randi(1000000,Gen_N,Inst_N); % Seeds for the RNG for the instance simulation
Protocol_seeds=randi(1000000,Gen_N,Inst_N);

%triplet(1:Gen_N,1:Inst_N,1:3)=NaN; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for G=1:Gen_N
    Mut_seed=Zoo.Mutation_seeds(G,:); % The workers for the parallel processes work in misterious ways, this slicing should speed up the procces a lot
    Unfold_seed=Zoo.Unfolding_seeds(G,:);
    Sim_seed=Zoo.Simulation_seeds(G,:);
    Prot_seed=Protocol_seeds(G,:);
    parfor C=1:Inst_N
    %for C=1:Inst_N
        %Mut_seed(C)
        [NewNet2{C}, Transition_tracker{C}]=mutation_routine(Primogen{C},Point_Mut,struct_mut,Mut_seed(C));
        
        %Protocol{C}=Protocol_maker1(Prot_seed(C));
        %Protocol{C}=Protocol_maker_differential(Prot_seed(C));
        %Protocol{C}=Protocol_maker_states(Prot_seed(C));
        %Protocol{C}=Protocol_maker_prob(Prot_seed(C));
        %Protocol{C}=Protocol_maker_NumericalInputAssesment(Prot_seed(C));
        %Protocol{C}=Protocol_maker_Independence(Prot_seed(C));
        %Protocol{C}=Protocol_maker_modstates(Prot_seed(C));
        %Protocol{C}=Protocol_maker_4subclassCPG(Prot_seed(C));
        %Protocol{C}=Protocol_maker_duallogic(Prot_seed(C))
        %Protocol{C}=Protocol_maker_prob(Prot_seed(C),0.2)
        %Protocol{C}=Protocol_maker_Sleep(Prot_seed(C))
        %Protocol{C}=Protocol_maker_order(Prot_seed(C));
        %Protocol{C}=Protocol_maker_long(Prot_seed(C));
        Protocol{C}=Protocol_maker_simplepattern(Prot_seed(C));
        %Protocol{C}=Protocol_maker_simplepattern_off(Prot_seed(C));
        %Protocol{C}=Protocol_maker_NumericalInputAssesment_20inp(Prot_seed(C))
%         Input=zeros(3,size(Protocol.Input{1},2));
%         for Q=1:3
%             Input(Q,:)=Protocol.Input{Q};
%         end
        Input{C}=Protocol{C}.Input;
        dur(C)=size(Input{C}{1},2);
    
        %triplet1(:,C)=[Mut_seed(C) Unfold_seed(C) Sim_seed(C)];
        [MyInst{C},Runtime{C},tech{C}]=unfold_v2(NewNet2{C},PSP_db,Plast_db,Unfold_seed(C)); % Unfolding genotype into an instance

        [Raster{C}, input_markdown, extras]=simulate_instance(MyInst{C},NaN,Runtime{C},tech{C},dur(C),Sim_seed(C),Input{C},[]);

        [Desired_output{C}, weights{C}]=make_desired_output(Protocol{C},tech{C});
%         for out=1:size(Protocol{C}.Output,2)
%             Desired_output{C}{out}=[zeros(size(Protocol{C}.Output{out},1),tech{C}.Longest_shape) Protocol{C}.Output{out}]; 
%         end
        
        % the extraction of the real output
        Real_output{C}=extract_real_output(Raster{C},MyInst{C});
%         class_pointer(C)=1;
%         for out_cl=1:size(MyInst{C}.Source.Cells.Output,2) % for each output subclass
%             Real_output{C}{out_cl}=Raster{C}(sum(MyInst{C}.Source.Cells.Input(1,:))+class_pointer(C):sum(MyInst{C}.Source.Cells.Input(1,:))+class_pointer(C)-1+MyInst{C}.Source.Cells.Output(1,out_cl),:);
%             class_pointer(C)=class_pointer(C)+MyInst{C}.Source.Cells.Output(1,out_cl);
%         end
        
        
        %Real_output=Raster{C}(12,:);
        
        fitness(C)=fit_func1(Desired_output{C},Real_output{C},weights{C},MyInst{C},Raster{C});
        
        
        
    end
    
    %triplet(G,:,:)=triplet1';
    % The 
    Zoo.Protocols(G,:)=Protocol;
    Zoo.Genotypes(G,:)=NewNet2;
    Zoo.Inst(G,:)=MyInst;
    Zoo.Rasters(G,:)=Raster;
    Zoo.Inputs(G,:)=Input; 
    Zoo.Outputs(G,:)=Desired_output;
    Zoo.Fitnesses(G,:)=fitness;
    Zoo.Weights(G,:)=weights;
    Zoo.Mutations(G,:)=Transition_tracker;
    
    % Fitness scramble 
    %fitness=fitness(randperm(Inst_N)); % Uncomment to randomly reassign fitness function values within the generation and thus scramble the selection process
    
    winner_order=sortrows([fitness' [1:Inst_N]'],'descend');
    maxfit(G)=winner_order(1,1);
    winner(G)=winner_order(1,2);
   
    
    %winner(G)=find(fitness==maxfit(G),1);
%     for P=1:Winners_N % Equal offsprings from top N performers
%         Primogen((P-1)*N_offspring+1:P*N_offspring)=Zoo.Genotypes(G,winner_order(P,2));
%         Zoo.Parents(G,(P-1)*N_offspring+1:P*N_offspring)=winner_order(P,2);
%     end
    %offsprings=frac_offspring(winner_order); % returns a vector with IDs of the members of the current generation that will leave the offspring that proceeds into the next generation 
    %offsprings=offspring_allocation(winner_order,"winner",1);
    offsprings=offspring_allocation(winner_order,"proportional");

    Zoo.Parents(G,:)=offsprings; % The indices for the offsprings 
    for C=1:Inst_N
        Primogen(C)=Zoo.Genotypes(G,offsprings(C));
    end

    title=['Generation ' num2str(G) ', fitness ' num2str(fitness(winner(G)))];
    if or(mod(G,100)==0,G==1)
        figure
        ax1=subplot(2,2,1:2);
        Raster_plot1(ax1,Zoo.Inst{G,winner(G)},Zoo.Rasters{G,winner(G)},Zoo.Outputs(G,winner(G)),Zoo.Weights(G,winner(G)),title)
        ax2=subplot(2,2,3);
        instance_vis(ax2,Zoo.Inst{G,winner(G)})
        disp(['Generation ' num2str(G)])
%         toc 
%         tic
    end
    
    
    Time(G)=toc;
    if mod(G,10)==0
        disp(['Elapsed time, 10 generations: ' num2str(sum(Time(G-9:G)))])
    end
    tic
end
figure
subplot(2,1,1)
plot(maxfit,'color','r','linewidth',3);
hold on
for G=1:Gen_N
   plot(G+randn(1,Inst_N-1)/5,Zoo.Fitnesses(G,[1:winner(G)-1   winner(G)+1:end]),'.k') 
end
axis([0 G+1 -1.1 1.1])
subplot(2,1,2)
plot(Time,'color','k','linewidth',3);
hold on
axis([0 G+1 0 max(Time)*1.1])
Zoo.Time=Time;
%triplet(end,:,:)

%find(Zoo.Parents(19,:)==winner(20),1)

% STANDALONE CALLS FOR PARTS OF THE ROUTINE: TO MAKE SEPARATE
% FUNCTIONS/SCRIPTS
% [TestNet2_1, Test_Transition_tracker_1]=mutation_routine(Zoo.Genotypes{Gen_N-1,Zoo.Parents(Gen_N-1,winner(Gen_N))},Point_Mut,struct_mut,Zoo.Mutation_seeds(Gen_N,winner(Gen_N)));
% [TestMyInst_1,Test_Runtime_1,Test_tech_1]=unfold(TestNet2_1,PSP_db,Plast_db,Zoo.Unfolding_seeds(Gen_N,winner(Gen_N)));
% Test_Input_1=Zoo.Inputs(Gen_N,winner(Gen_N));
% Test_dur_1=size(Test_Input_1{1}{1},2);
% [Test_Raster_1, Test_input_markdown_1, Test_extras_1]=simulate_instance(TestMyInst_1,Test_Runtime_1,Test_tech_1,Test_dur_1,Zoo.Simulation_seeds(Gen_N,winner(Gen_N)),Test_Input_1{1},[]);
% [Test_Desired_output_1, Test_weights_1]=make_desired_output(Zoo.Protocols{Gen_N,winner(Gen_N)},Test_tech_1);
% Test_Real_output_1=extract_real_output(Test_Raster_1,TestMyInst_1);
% title=['Reproducibility test 1'];
% figure
%         ax1=subplot(2,2,1:2);
%         Raster_plot1(ax1,TestMyInst_1,Test_Raster_1,{Test_Desired_output_1},{Test_weights_1},title)
%         ax2=subplot(2,2,3);
%         instance_vis(ax2,TestMyInst_1)
%         
% [TestNet2_2, Test_Transition_tracker_2]=mutation_routine(Zoo.Genotypes{Gen_N-1,Zoo.Parents(Gen_N-1,winner(Gen_N))},Point_Mut,struct_mut,Zoo.Mutation_seeds(Gen_N,winner(Gen_N)));
% [TestMyInst_2,Test_Runtime_2,Test_tech_2]=unfold(TestNet2_2,PSP_db,Plast_db,Zoo.Unfolding_seeds(Gen_N,winner(Gen_N)));
% Test_Input_2=Zoo.Inputs(Gen_N,winner(Gen_N));
% Test_dur_2=size(Test_Input_2{1}{1},2);
% [Test_Raster_2, Test_input_markdown_2, Test_extras_2]=simulate_instance(TestMyInst_2,Test_Runtime_2,Test_tech_2,Test_dur_2,Zoo.Simulation_seeds(Gen_N,winner(Gen_N)),Test_Input_2{1},[]);
% [Test_Desired_output_2, Test_weights_2]=make_desired_output(Zoo.Protocols{Gen_N,winner(Gen_N)},Test_tech_2);
% Test_Real_output_2=extract_real_output(Test_Raster_2,TestMyInst_2);
% title=['Reproducibility test 2'];
% figure
%         ax1=subplot(2,2,1:2);
%         Raster_plot1(ax1,TestMyInst_2,Test_Raster_2,{Test_Desired_output_2},{Test_weights_2},title)
%         ax2=subplot(2,2,3);
%         instance_vis(ax2,TestMyInst_2)
%         
%         
%       
%         
%         
%         
% 
% [TestMyInst_3,Test_Runtime_3,Test_tech_3]=unfold(Zoo.Genotypes{Gen_N,winner(Gen_N)},PSP_db,Plast_db,Zoo.Unfolding_seeds(Gen_N,winner(Gen_N)));
% Test_Input_3=Zoo.Inputs(Gen_N,winner(Gen_N));
% Test_dur_3=size(Test_Input_3{1}{1},2);
% [Test_Raster_3, Test_input_markdown_3, Test_extras_3]=simulate_instance(TestMyInst_3,Test_Runtime_3,Test_tech_3,Test_dur_3,Zoo.Simulation_seeds(Gen_N,winner(Gen_N)),Test_Input_3{1},[]);
% [Test_Desired_output_3, Test_weights_3]=make_desired_output(Zoo.Protocols{Gen_N,winner(Gen_N)},Test_tech_3);
% Test_Real_output_3=extract_real_output(Test_Raster_3,TestMyInst_3);
% title=['Reproducibility test 3'];
% figure
%         ax1=subplot(2,2,1:2);
%         Raster_plot1(ax1,TestMyInst_3,Test_Raster_3,{Test_Desired_output_3},{Test_weights_3},title)
%         ax2=subplot(2,2,3);
%         instance_vis(ax2,TestMyInst_3)
%         
        
        
        
% Saving
% Since the whole protocol structures are saved in the Zoo, separate inputs and putputs are not required
Zoo=rmfield(Zoo,'Inputs');
Zoo=rmfield(Zoo,'Outputs');
if Save_mode==2
    Zoo=rmfield(Zoo,'Inst');
    Zoo=rmfield(Zoo,'Rasters');
elseif Save_mode==3
    Zoo=rmfield(Zoo,'Inst');
    Zoo=rmfield(Zoo,'Rasters');
    Zoo=rmfield(Zoo,'Genotypes');
end
out=Zoo;
%savename=['C:\Users\Zavrd\Documents\Zoo_NXOR_mode2_500gen_20Inst_' num2str(outeriter) '.mat'];
%savename=['C:\Users\Zavrd\Documents\Zoo_Separ_mode2_300gen_10Inst_' num2str(run) '.mat'];
%savename=['C:\Users\Zavrd\Documents\Zoo_DiffOR_mode2_1000gen_10Inst_1.mat'];
%save(savename,'Zoo','-v7.3')        
% %end       




end