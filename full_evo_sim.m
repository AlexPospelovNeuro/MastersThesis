
% The Testing routine functions script turned into a function

function out=full_evo_sim(run,Point_Mut)

Save_mode=2;
% Save_mode = 1 - save all the genotypes, phenotypes and rasters, leads to a very big export variables, but direct interaction with the results later on results
% Save_mode = 2 - save only the genotypes. Everything else can be resimulated using the seeds
% Save_mode = 3 - don't save even the genotypes, only the initial settings and the seeds. Smallest files, but require full resimulation of the evolution

load('MyNet2_I_111_O_1_NotEvolved.mat');
%load('MyNet_Evolved_I_1_1_1_O_1_probabilistic.mat');
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
struct_mut=[0.2 0.025 0.2 0.05 0.2 0.05 0.2 0.05]; %[p_cl_div p_cl_rem p_plast_add p_plast_remove p_exp_add p_exp_remove p_link_add p_link_remove]


Gen_N=100; % Number of generations
Inst_N=10; % Number of instances per generation
for init=1:Inst_N
    for q=1:Zoo.Settings.PreMutationsN
        [Net2, ~]=mutation_routine(Net2,Point_Mut,struct_mut,NaN);
    end
    legality(Net2)
    Primogen(init)={Net2};
end

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

tic
for G=1:Gen_N
    Mut_seed=Zoo.Mutation_seeds(G,:); % The workers for the parallel processes work in misterious ways, this slicing should speed up the procces a lot
    Unfold_seed=Zoo.Unfolding_seeds(G,:);
    Sim_seed=Zoo.Simulation_seeds(G,:);
    Prot_seed=Protocol_seeds(G,:);
    parfor C=1:Inst_N
    %for C=1:Inst_N   % uncomment for the troubleshooting or if parfor is unavailable. Reciprocally, comment the previous line. 

        [NewNet2{C}, Transition_tracker{C}]=mutation_routine(Primogen{C},Point_Mut,struct_mut,Mut_seed(C));
        
        Protocol{C}=Protocol_maker_bin(Prot_seed(C));
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
        %Protocol{C}=Protocol_maker_simplepattern(Prot_seed(C));
        %Protocol{C}=Protocol_maker_simplepattern_off(Prot_seed(C));
        %Protocol{C}=Protocol_maker_NumericalInputAssesment_20inp(Prot_seed(C))

        Input{C}=Protocol{C}.Input;
        dur(C)=size(Input{C}{1},2);
    

        [MyInst{C},Runtime{C},tech{C}]=unfold_v2(NewNet2{C},PSP_db,Plast_db,Unfold_seed(C)); % Unfolding genotype into an instance

        [Raster{C}, input_markdown, extras]=simulate_instance(MyInst{C},NaN,Runtime{C},tech{C},dur(C),Sim_seed(C),Input{C},[]);

        [Desired_output{C}, weights{C}]=make_desired_output(Protocol{C},tech{C});

        
        % the extraction of the real output
        Real_output{C}=extract_real_output(Raster{C},MyInst{C});

        fitness(C)=fit_func1(Desired_output{C},Real_output{C},weights{C},MyInst{C},Raster{C});
        
        
        
    end
    
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
   

    %offsprings=offspring_allocation(winner_order,"winner",1); % Alternative offspring allocation mode, to show the interface
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



figure
plot(maxfit,'color','r','linewidth',3);
hold on
for G=1:Gen_N
   plot(G+randn(1,Inst_N-1)/5,Zoo.Fitnesses(G,[1:winner(G)-1   winner(G)+1:end]),'.k') 
end
axis([0 G+1 -1.1 1.1])


        
% Saving
% Since the whole protocol structures are saved in the Zoo, separate inputs and putputs are not required
% Currenly, the full_evo_sim_envelope does the saving after getting Zoo variable from here. But the code below is an alternative for this approach. 
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

end