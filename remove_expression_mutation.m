% A function that removes expressions with probability p
% INPUTS:
% - A properly structured Net2 genotype
% - A probability p, applied to each available expression
% OUTPUTS:
% - A properly structured Net2 genotype


function [NewNet2,transition_tracker]=remove_expression_mutation(Net2,p,transition_tracker)

NewNet2=Net2;  
if isfield(Net2,'Expression_patterns')
    N_exp=size(NewNet2.Expression_patterns,2);
    for exp=1:N_exp
        if rand(1)<p
            NewNet2=remove_expression_determ(NewNet2,N_exp-exp+1); % The deletion is done in the reversed order so shifting if the indexing after the deletion is not required
            if isstruct(transition_tracker)
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)==N_exp-exp+1)=NaN;
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)=transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)-1;
            end
        end
    end
else
if isempty(NewNet2,'Expression_patterns')
    NewNet2=rmfield(NewNet2,'Expression_patterns');
end
%legality(NewNet2)
end