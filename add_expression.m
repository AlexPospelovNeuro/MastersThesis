% Function adds a new expression to a Net2 structure and then - one link to
% this new expression

% INPUTS:
% - A properly structured Net2 genotype
% - A probability p (each existing pattern will have exactly one new link with probability p)
% OUTPUTS:
% - A new Net2 genotype, structured similarly to the input (same set of the interneuron classes, plasticity patterns and expression patterns)


function [NewNet2,transition_tracker]=add_expression(Net2,p,transition_tracker)
NewNet2=Net2;
if rand(1)<p
    if isfield(NewNet2,'Expression_patterns')
        Index=size(NewNet2.Expression_patterns,2)+1;
    else
        Index=1;
    end
    NewNet2.Expression_patterns{Index}{1}=randi([1 size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2)]); % Expression pattern can appear in any subclass including in and out;
    NewNet2=add_link(NewNet2,Index,1);
    if isstruct(transition_tracker)
        if isempty(transition_tracker.expressions)
            transition_tracker.expressions(1,1)=NaN;
            transition_tracker.expressions(2,1)=1;
            transition_tracker.links{1}=[NaN;2];
        else
            transition_tracker.expressions(1,end+1)=NaN;
            transition_tracker.expressions(2,end)=max(transition_tracker.expressions(2,:))+1;
            %transition_tracker.links{max(transition_tracker.expressions(2,:))}=[NaN;2];
            transition_tracker.links{size(transition_tracker.expressions,2)}=[NaN;2]; % Since the links from the removed expressiosn are kept in the transition table, the new ones have to be appended to the "exhausive" list
        end
    end
end
%legality(NewNet2)    
    
end