% A function that removes a selected plasticity pattern and all the
% expression patterns that refer to it (so, it is a complete removal of the plasticity)

% INPUTS:
% - A properly structured Net2 genotype
% - An index of a plasticiyu pattern
% OUTPUTS:
% - A properly structured Net2 genotype




function [NewNet2,transition_tracker]=remove_plasticity_determ(Net2,plast,transition_tracker)

if ~isfield(Net2,'Plasticity')
    error('There are no plasticities in the genotype')
end
if plast>size(Net2.Plasticity,2)
     error('The target plasticity does not exist')
end
NewNet2=Net2;  
NewNet2.Plasticity(plast)=[];

if isfield(NewNet2,'Expression_patterns') % If there are expression patterns, some of them may refer to the deleted plasticity
    N_exp=size(NewNet2.Expression_patterns,2);
    for exp=1:N_exp

        N_link=size(NewNet2.Expression_patterns{N_exp-exp+1},2); % going from the last to the first to not have to adjust indexing
        for pat=2:N_link
            if NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(1)==6 % if the link refers to a plasticity
                %disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                if NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)==plast
                    %disp('expression pattern removal')
                    NewNet2=remove_link_determ(NewNet2,N_exp-exp+1,N_link-pat+1);
                    if isstruct(transition_tracker)
                        transition_exp_index=find(transition_tracker.expressions(2,:)==N_exp-exp+1);
                        transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)==N_link-pat+2)=NaN;
                        transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)=transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)-1;
                    end
                elseif NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)>plast
                    NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)=NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)-1; % adjusting the index of plasticity with respect to a deleded one
                end
            end
        end
        
        if size(NewNet2.Expression_patterns{N_exp-exp+1},2)==1 % If there are no links left
            NewNet2=remove_expression_determ(NewNet2,N_exp-exp+1); % remove the whole expression
            if isstruct(transition_tracker)
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)==N_exp-exp+1)=NaN; 
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)=transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)-1;
            end
            %exps_affected=exps_affected-1;
        end
    end
end
if isempty(NewNet2.Plasticity)
    NewNet2=rmfield(NewNet2,'Plasticity');
end
%legality(NewNet2)
end