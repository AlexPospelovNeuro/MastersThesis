

% A function that removes links from expressions with probability p (deletion routine is applied for every link of every expression)
% INPUTS:
% - A properly structured Net2 genotype
% - A probability p, applied to each available expression
% OUTPUTS:
% - A properly structured Net2 genotype





function [NewNet2, transition_tracker]=remove_link_mutation(Net2,p,transition_tracker)
%report.links={};
%report.exps=[];
NewNet2=Net2;  
if isfield(NewNet2,'Expression_patterns')
    N_exp=size(NewNet2.Expression_patterns,2);
    for exp=1:N_exp
%         local_counter=1;
%         exps_affected=0;
        N_link=size(NewNet2.Expression_patterns{N_exp-exp+1},2); % going from the last to the first to not have to adjust indexing
        for pat=2:N_link
            if rand(1)<p
                NewNet2=remove_link_determ(NewNet2,N_exp-exp+1,N_link-pat+1);
                if isstruct(transition_tracker)
                    transition_exp_index=find(transition_tracker.expressions(2,:)==N_exp-exp+1);
                    transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)==N_link-pat+2)=NaN;
                    transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)=transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)-1;
                end
                
                
%                 if local_counter==1
%                     report{exps_affected+1}{1}=exp;
%                 end
%                 report{exps_affected+1}{2}(local_counter)=pat;
%                 local_counter=local_counter+1;
            end
            
        end
%         if local_counter>1
%             exps_affected=exps_affected+1;
%         end
        if size(NewNet2.Expression_patterns{N_exp-exp+1},2)==1 % If there are no links left
            NewNet2=remove_expression_determ(NewNet2,N_exp-exp+1); % remove the whole expression
            if isstruct(transition_tracker)
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)==N_exp-exp+1)=NaN; 
                transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)=transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)-1;
            end
%             report{exps_affected}={};
%             exps_affected=exps_affected-1;
        end
    end
end
%legality(Net2)
end