% The function that forms the tracker variable that keeps track of the indexadions' change due to the mutations (except the point mutations). Each branch of the tracker variable consists of two rows: the top row
% does not change unless appended with a new item (in this case, the new
% value in the top row is NaN for de novo item generation, or the parent
% index for the division). The second row is being changed according to the
% indexation change. This way, in the end of the mutation routine, the
% second row will contain the new indices of items in the new genotype,
% while the old - the indices of the same items in the old genotype (or
% NaNs, if there were no such items). Sanity check is: the indexation of
% any type of items in the new genotype should be valid: start with 1 and
% and with total number of items of the type, one number per item. Their
% order in theory can be any.



function transition_tracker=make_transition_tracker(Net2)
transition_tracker.subclasses=[1:size(Net2.Connections,1);1:size(Net2.Connections,1)];
if isfield(Net2,'Plasticity')
    transition_tracker.plasticities=[1:size(Net2.Plasticity,2);1:size(Net2.Plasticity,2)];
else
    transition_tracker.plasticities=[[];[]];
end
if isfield(Net2,'Expression_patterns')
    transition_tracker.expressions=[1:size(Net2.Expression_patterns,2);1:size(Net2.Expression_patterns,2)];
    for exp=1:size(Net2.Expression_patterns,2)
        transition_tracker.links{exp}=[2:size(Net2.Expression_patterns{exp},2);2:size(Net2.Expression_patterns{exp},2)];
    end
else
    transition_tracker.expressions=[[];[]];
    transition_tracker.links={};
end



end