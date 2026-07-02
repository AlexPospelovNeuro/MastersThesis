

% The script that extracts and saves a circuit by resimulation for the
% evolutionary Zoo saved
clear 

load('C:\Users\Zavrd\Documents\Zoo_bigquantconcave_mode2_100gen_10Inst_1.mat');

Net2=Zoo.Genotypes{end,Zoo.Fitnesses(end,:)==max(Zoo.Fitnesses(end,:))};





% 
Net2.Info="IO structure: I 1 20 O 1, 20-neuron input size evaluation circuitm concave profile, source tree Zoo_bigquantconcave_mode2_100gen_10Inst_1";
save('MyNet2_I_1_20_O_1_QuantConcave_Evolved.mat','Net2')

PSP_db=load('PSP_shapes1.mat');
Plast_db=load('plast_shapes1.mat');
% Resimulation

Protocol{1}=Zoo.Protocols{end,Zoo.Fitnesses(end,:)==max(Zoo.Fitnesses(end,:))};
Input{1}=Protocol{1}.Input;
dur(1)=size(Input{1}{1},2);
[MyInst{1},Runtime{1},tech{1}]=unfold_v2(Net2,PSP_db,Plast_db,Zoo.Unfolding_seeds(end,Zoo.Fitnesses(end,:)==max(Zoo.Fitnesses(end,:)))); % Unfolding genotype into an instance
[Raster{1}, input_markdown, extras]=simulate_instance(MyInst{1},Runtime{1},tech{1},dur(1),Zoo.Simulation_seeds(end,Zoo.Fitnesses(end,:)==max(Zoo.Fitnesses(end,:))),Input{1},[]);
[Desired_output{1}, weights{1}]=make_desired_output(Protocol{1},tech{1});
Real_output{1}=extract_real_output(Raster{1},MyInst{1});

fitness_resim(1)=fit_func1(Desired_output{1},Real_output{1},weights{1},MyInst{1},Raster{1});

title=['Resimulation, fitness ' num2str(fitness_resim(1))];
figure
ax1=subplot(2,2,1:2);
Raster_plot1(ax1,MyInst{1},Raster{1},Desired_output(1),weights(1),title)
ax2=subplot(2,2,3);
instance_vis(ax2,MyInst{1})


% Phenotypic variability: same genotype, different phenotypes

N_I=10; % Number of instances 
Protocol_Phenvar{1}=Protocol_maker_NumericalInputAssesment(NaN); % Make sure that the protocol is correct and uses the same parameters. Otherwise, make a new protocol form the ones saved in Zoo.
Input_Phenvar{1}=Protocol_Phenvar{1}.Input;
dur(1)=size(Input_Phenvar{1}{1},2);
for I=1:N_I
    
    
    [MyInst_Phenvar{I},Runtime_Phenvar{I},tech_Phenvar{I}]=unfold_v2(Net2,PSP_db,Plast_db,NaN); % Unfolding genotype into an instance with rando rng seeds]
    [Raster_Phenvar{I}, input_markdown_Phenvar, extras_Phenvar]=simulate_instance(MyInst_Phenvar{I},Runtime_Phenvar{I},tech_Phenvar{I},dur(1),NaN,Input_Phenvar{1},[]);
    [Desired_output_Phenvar{I}, weights_Phenvar{I}]=make_desired_output(Protocol_Phenvar{1},tech_Phenvar{1});
    Real_output_Phenvar{I}=extract_real_output(Raster_Phenvar{I},MyInst_Phenvar{I});
    fitness_Phenvar(I)=fit_func1(Desired_output_Phenvar{I},Real_output_Phenvar{I},weights_Phenvar{I},MyInst_Phenvar{I},Raster_Phenvar{I});
    
    
    
    title=['Phenotypic variability, fitness ' num2str(fitness_Phenvar(I))];
    figure
    ax1=subplot(2,2,1:2);
    Raster_plot1(ax1,MyInst_Phenvar{I},Raster_Phenvar{I},Desired_output_Phenvar(I),weights_Phenvar(I),title)
    ax2=subplot(2,2,3);
    instance_vis(ax2,MyInst_Phenvar{I})
    
    
end
clear MyInst_Phenvar Runtime_Phenvar tech_Phenvar Raster_Phenvar input_markdown_Phenvar extras_Phenvar Desired_output_Phenvar weights_Phenvar Real_output_Phenvar


% The "follow curves"
% From each generation, we choose 5 genotypes that occupy [1 6 11 16 20]
% places in the fitness ranking. We resimulate each of them, and also do 10
% phenotypic variability simulations for each. 

ranks=[1 6 11 16 20]; % The ranks of the genotypes to resimulate
N_Var=100; % The number of phenotypic varians for each genotype to run
winner_order(1:size(Zoo.Genotypes,1),1:size(Zoo.Genotypes,2))=NaN;
fitness_resim(1:size(Zoo.Genotypes,1),1:size(ranks,2))=NaN;
fitness_Phenvar_evo(1:size(Zoo.Genotypes,1),1:size(ranks,2),1:N_Var)=NaN;
for G=1:size(Zoo.Genotypes,1) % for each generation
    wo=sortrows([Zoo.Fitnesses(G,:)' [1:size(Zoo.Genotypes,2)]'],'descend');
    winner_order(G,:)=wo(:,2);
    for inst=1:size(ranks,2) % For each one of five selected instances
        Net2=Zoo.Genotypes{G,winner_order(G,ranks(inst))};
        Protocol{1}=Zoo.Protocols{G,winner_order(G,ranks(inst))};
        Input{1}=Protocol{1}.Input;
        dur(1)=size(Input{1}{1},2);
        [MyInst{1},Runtime{1},tech{1}]=unfold_v2(Net2,PSP_db,Plast_db,Zoo.Unfolding_seeds(G,winner_order(G,ranks(inst)))); % Unfolding genotype into an instance
        [Raster{1}, input_markdown, extras]=simulate_instance(MyInst{1},Runtime{1},tech{1},dur(1),Zoo.Simulation_seeds(G,winner_order(G,ranks(inst))),Input{1},[]);
        [Desired_output{1}, weights{1}]=make_desired_output(Protocol{1},tech{1});
        Real_output{1}=extract_real_output(Raster{1},MyInst{1});
        fitness_resim(G,inst)=fit_func1(Desired_output{1},Real_output{1},weights{1},MyInst{1},Raster{1});
        
        
        
        % Phenotypic variability simulation with the same input that was used for the resimulation (and the original simulation)
        
        parfor v=1:N_Var
            
            [MyInst_Phenvar{v},Runtime_Phenvar{v},tech_Phenvar{v}]=unfold_v2(Net2,PSP_db,Plast_db,NaN); % Unfolding genotype into an instance with rando rng seeds]
            [Raster_Phenvar{v}, input_markdown_Phenvar, extras_Phenvar]=simulate_instance(MyInst_Phenvar{v},Runtime_Phenvar{v},tech_Phenvar{v},dur(1),NaN,Input{1},[]);
            [Desired_output_Phenvar{v}, weights_Phenvar{v}]=make_desired_output(Protocol{1},tech_Phenvar{v});
            Real_output_Phenvar{v}=extract_real_output(Raster_Phenvar{v},MyInst_Phenvar{v});
            fitness_Phenvar_evo_slice(v)=fit_func1(Desired_output_Phenvar{v},Real_output_Phenvar{v},weights_Phenvar{v},MyInst_Phenvar{v},Raster_Phenvar{v});
            
            
        end
        fitness_Phenvar_evo(G,inst,:)=fitness_Phenvar_evo_slice;
    end
    disp(['Generation ' num2str(G)])
end


colors={[0 1 0],[0 0.75 0.25],[0 0.5 0.5],[0 0.25 0.75],[0 0 1]};
figure
plot(Zoo.Fitnesses,'.k')
hold on
for inst=1:size(ranks,2)
    plot(fitness_resim(:,inst),'linewidth',1,'color',colors{inst},'linestyle','--')
    
    %fitness_Phenvar_evo_forplot=fitness_Phenvar_evo(:,inst,:);
    plot(squeeze(fitness_Phenvar_evo(:,inst,:)),'linestyle','none','Marker','.','color',colors{inst},'Markersize',10)
end


fitness_Phenvar_quant=quantile(fitness_Phenvar_evo,[0.1 0.9],3);
figure
plot(Zoo.Fitnesses,'.k')
hold on
for inst=1:size(ranks,2)
    plot(fitness_resim(:,inst),'linewidth',2,'color',colors{inst},'linestyle','-')
    
    %fitness_Phenvar_evo_forplot=fitness_Phenvar_evo(:,inst,:);
    %plot(squeeze(fitness_Phenvar_evo(:,inst,:)),'linestyle','none','Marker','.','color',colors{inst},'Markersize',10)
    %fill([1:size(Zoo.Genotypes,1) fliplr(1:size(Zoo.Genotypes,1))]',[fitness_Phenvar_quant(:,inst,1)' fliplr(fitness_Phenvar_quant(:,inst,2))']','FaceColor',colors{inst},'EdgeColor',colors{inst},'FaceAlpha',0.2,'EdgeAlpha',1)
    fill([1:size(Zoo.Genotypes,1) fliplr(1:size(Zoo.Genotypes,1))]',[fitness_Phenvar_quant(:,inst,1)' fliplr(fitness_Phenvar_quant(:,inst,2)')]','r','FaceColor',colors{inst},'EdgeColor',colors{inst},'FaceAlpha',0.2,'EdgeAlpha',1)
end



 Fit_var(1:size(Zoo.Genotypes,1),1:size(ranks,2))=NaN;
for G=1:size(Zoo.Genotypes,1) % for each generation
    for inst=1:size(ranks,2) % For each one of five selected instances
        Fit_var(G,inst)=var(fitness_Phenvar_evo(G,inst,:));
    end
end

figure
plot(Zoo.Fitnesses,'.k')
hold on
for inst=1:size(ranks,2)
    plot(fitness_resim(:,inst),'linewidth',2,'color',colors{inst},'linestyle','-')
    
    %fitness_Phenvar_evo_forplot=fitness_Phenvar_evo(:,inst,:);
    %plot(squeeze(fitness_Phenvar_evo(:,inst,:)),'linestyle','none','Marker','.','color',colors{inst},'Markersize',10)
    %fill([1:size(Zoo.Genotypes,1) fliplr(1:size(Zoo.Genotypes,1))]',[fitness_Phenvar_quant(:,inst,1)' fliplr(fitness_Phenvar_quant(:,inst,2))']','FaceColor',colors{inst},'EdgeColor',colors{inst},'FaceAlpha',0.2,'EdgeAlpha',1)
    fill([1:size(Zoo.Genotypes,1) fliplr(1:size(Zoo.Genotypes,1))]',[fitness_resim(:,inst)'-Fit_var(:,inst)' fliplr(fitness_resim(:,inst)'+Fit_var(:,inst)')]','r','FaceColor',colors{inst},'EdgeColor',colors{inst},'FaceAlpha',0.2,'EdgeAlpha',1)
end



figure
plot(fitness_resim,Fit_var,'.k')
fitness_track.orig=fitness_resim;
fitness_track.variability=fitness_Phenvar_evo;
save('C:\Users\Zavrd\Documents\Fenvar_Quant1_mode2_300gen_20Inst_4.mat','fitness_track','-v7.3')
