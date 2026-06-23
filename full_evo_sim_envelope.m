clear
load('point_mutation_template2.mat');
for run=3:5

    zoo=full_evo_sim(run,Point_Mut);
    %zoo=full_evo_learn_sim(run,Point_Mut);
    %zoo=full_evo_sim_heritable_diversity(run,Point_Mut);
    savename=['C:\Users\apospelo\Documents\Zoo_simplepattern_' num2str(run) '.mat'];
    save(savename,'zoo','-v7.3')

end

