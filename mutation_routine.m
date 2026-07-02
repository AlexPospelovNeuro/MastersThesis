% This function applies a full mutation routine to the genotype. The
% probabilities are hardcoded at this point, beware

% The order of mutation applications:
% - Divide class
% - Delete class
% - Duplicate plasticity
% - Delete plasticity
% - Create expression
% - Remove expression
% - Create link
% - Remove link
% - Apply point mutations

% INPUTS:
% - A properly structured Net2 genotype
% OUTPUTS:
% - A properly structured Net2 genotype



function [NewNet2, transition_tracker]=mutation_routine(Net2,Point_Mut,struct_mut,seed)

if ~isnan(seed)
    rng(seed, 'twister');
end
NewNet2=Net2;
% The transition_table variable contains the track of all "structural" mutations that happened during the routine in form of the index changes. For each index that is in the "input" genotype (for subclass, plasticity, expression or link), it will contain a new index of the same entity. 
% If the Index stayed the same, the "new" index will be equal to the old one. If the entity was removed, the new index will be NaN. If the entity appeared de novo, the "old" index will be NaN. If the entity appeared due to a duplication, there will be few (2 or 4) "new" indices for one "old"
% for some of the entities. Following the transition tables by a backtracing will alow to track the same entity from generation to generation even if its index changed during the evolutionary routine. 

transition_tracker=make_transition_tracker(NewNet2);

% Division of one random subclass
p_div=struct_mut(1);
if rand<p_div
    target_class=randi([size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+1 size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]);
    [NewNet2,transition_tracker]=divide_subclass_determ(NewNet2,target_class,transition_tracker);
end


% Deletion of subclasses (applied to each one)
p_del=struct_mut(2); % If p_del is at least as big as p_div/2 (given that smallest number of interneuron subclasses is two, one ionotropic and one modulatory), the overall tendency of random mutation process will be to delete the subclasses (on average). This process will be kept in check only by the negative selection.
N_max=(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2));
N_reserved=(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2));
for target_class=(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+1):N_max
    %target_class
    if rand<p_del
        [NewNet2,transition_tracker]=remove_subclass_determ(NewNet2,N_max-target_class+N_reserved+1,transition_tracker); % Going through the subclasses from the last to the first to reduce the reindexation
    end
end

% Creation of a plasticity (random type, random synapse)
p_newplast=struct_mut(3); % Add one plasticity pattern. Since the plasticity is available for every synapse (so, number of subclasses squared), it makes sense to keep probability of a new plasticity creation slightly higher
[NewNet2,transition_tracker]=add_plasticity(NewNet2,p_newplast,transition_tracker); 


% Removal of a plasticity
p_removeplast=struct_mut(4); % Again, the removal is applied to each plasticity, so the overall tendency is to reduce the number of plasticityes at random
if isfield(NewNet2,'Plasticity')
    [NewNet2,transition_tracker]=remove_plasticity_mutation(NewNet2,p_removeplast,transition_tracker); % The mutation function applies the removal to each available plasticity
end

% Creation of the expression
p_newexpression=struct_mut(5); % Add one expression
[NewNet2,transition_tracker]=add_expression(NewNet2,p_newexpression,transition_tracker);

% Removal of the expression
p_removeexpression=struct_mut(6); % Remove each expression
if isfield(NewNet2,'Expression_patterns')
    [NewNet2,transition_tracker]=remove_expression_mutation(NewNet2,p_removeexpression,transition_tracker);
end

% Link creation. Adds one link with decided probability to each available
% expression pattern
p_newlink=struct_mut(7);
if isfield(NewNet2,'Expression_patterns')
    N_exp=size(NewNet2.Expression_patterns,2);
    for exp=1:N_exp
        [NewNet2,flag]=add_link(NewNet2,exp,p_newlink);
        if flag==1
            if isstruct(transition_tracker)
                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(1,end+1)=NaN;
                transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,end)=max(transition_tracker.links{transition_tracker.expressions(2,:)==exp}(2,:))+1;
            end
        end
    end
end

% Link removal 
p_removelink=struct_mut(8);
if isfield(NewNet2,'Expression_patterns')
    N_exp=size(NewNet2.Expression_patterns,2);
    for exp=1:N_exp
        N_pattern=size(NewNet2.Expression_patterns{N_exp-exp+1},2)-1;
        for pat=1:N_pattern
            if rand<p_removelink
                NewNet2=remove_link_determ(NewNet2,N_exp-exp+1,N_pattern-pat+1); % reversed order for both expression and pattern
                if isstruct(transition_tracker)
                    transition_exp_index=find(transition_tracker.expressions(2,:)==N_exp-exp+1);
                    transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)==N_pattern-pat+2)=NaN;
                    transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_pattern-pat+2)=transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_pattern-pat+2)-1;
                end
            end
        end
        if size(NewNet2.Expression_patterns{N_exp-exp+1},2)==1 % If there are no links left
            NewNet2=remove_expression_determ(NewNet2,N_exp-exp+1); % remove the whole expression
            if isstruct(transition_tracker)
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)==N_exp-exp+1)=NaN;
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)=transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)-1;
            end
        end
    end
    if isempty(NewNet2.Expression_patterns)
        NewNet2=rmfield(NewNet2,'Expression_patterns');
    end
end

% apply point mutations. The probabilities of mutations and their
% magnitudes are in the Point_Mut structure
NewNet2=point_mutation(NewNet2,Point_Mut);
%legality(NewNet2)
end