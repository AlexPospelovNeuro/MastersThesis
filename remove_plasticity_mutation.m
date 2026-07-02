% A function that removes plasticity with probability p
% INPUTS:
% - A properly structured Net2 genotype
% - A probability p, applied to each available plasticity
% OUTPUTS:
% - A properly structured Net2 genotype



function [NewNet2,transition_tracker]=remove_plasticity_mutation(Net2,p,transition_tracker)

NewNet2=Net2;  
if isfield(Net2,'Plasticity')
    N_plast=size(NewNet2.Plasticity,2);
    for plast=1:N_plast
        if rand(1)<p
            [NewNet2,transition_tracker]=remove_plasticity_determ(NewNet2,N_plast-plast+1,transition_tracker); % The deletion is done in the reversed order so shifting if the indexing after the deletion is not required
            if isstruct(transition_tracker)
                transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)==N_plast-plast+1)=NaN; % "void" the removed plasticity
                transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)>N_plast-plast+1)=transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)>N_plast-plast+1)-1; % reduce the indices of the plasticities with higher index
            end
        end
    end
else


%legality(NewNet2)
end