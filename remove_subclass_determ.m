% A function that removes a selected interneuron subclass and everything
% that refers to it: the features, the expression patterns and the
% plasticities

% INPUTS:
% - A properly structured Net2 genotype
% - An index of a subclass
% OUTPUTS:
% - A properly structured Net2 genotype


function [NewNet2,transition_tracker]=remove_subclass_determ(Net2,cl,transition_tracker)
NewNet2=Net2; 

% removed.subclasses=cl;
% removed.plast=[];
% removed.expressions=[];
% removed.links={};

if cl<=(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2))
    error('deletion of a reserved input and output neurons is forbidden')
end
if cl>(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2)+size(NewNet2.Cells.Mod,2))
    error('deletion of a non-existing subclass')
end
if cl>(size(NewNet2.Cells.Input,2)+size(NewNet2.Cells.Output,2)+size(NewNet2.Cells.Ion,2))
    type=2; % modulatory class
    Index=cl-size(NewNet2.Cells.Input,2)-size(NewNet2.Cells.Output,2)-size(NewNet2.Cells.Ion,2); % Index inside the class
else
    type=1; % ionotropic class
    Index=cl-size(NewNet2.Cells.Input,2)-size(NewNet2.Cells.Output,2); % Index inside the class
end

if ((type==1)&&(size(NewNet2.Cells.Ion,2)>1)) || ((type==2)&&(size(NewNet2.Cells.Mod,2)>1)) % The mutation cannot remove the last subclass of the class
    transition_tracker.subclasses(2,transition_tracker.subclasses(2,:)==cl)=NaN;
    transition_tracker.subclasses(2,transition_tracker.subclasses(2,:)>cl)=transition_tracker.subclasses(2,transition_tracker.subclasses(2,:)>cl)-1;
    if type==1
        NewNet2.Cells.Ion(:,Index)=[]; % removing the cell numbers
    elseif type==2
        NewNet2.Cells.Mod(:,Index)=[]; % removing the cell numbers
    else
        error('Unknown subclass to remove')
    end
    NewNet2.Connections(cl,:,:)=[]; % removing the ionotropic connection
    if type==1
        NewNet2.Connections(:,cl,:)=[];
    end
    NewNet2.Connections_Mod(cl,:,:,:)=[]; % removing the modulatory inputs
    if type==1
        NewNet2.Connections_Mod(:,cl,:,:)=[];
    elseif type==2
        NewNet2.Connections_Mod(:,:,Index,:)=[];
    end
    NewNet2.BasicPowers(cl,:,:)=[]; % Removing basic powers
    if type==1
        NewNet2.BasicPowers(:,cl,:)=[];
    end
    NewNet2.Delays(cl,:,:)=[]; % Removing delays
    if type==1
        NewNet2.Delays(:,cl,:)=[];
    end
    NewNet2.PSPshape(cl,:)=[]; % Removing PSP shapes
    if type==1
        NewNet2.PSPshape(:,cl)=[];
    end
    NewNet2.Thresholds(:,cl)=[]; % removing thresholds
    NewNet2.AbsRefract(:,cl)=[]; % removing absolute refracterity
    NewNet2.ThreshNoise(:,cl)=[]; % removing threshold noise
    NewNet2.RecurCon(:,cl)=[]; % removing recurrent connections
    
    NewNet2.Mod(cl,:,:)=[]; % Removing modulation parameters
    if type==1
        NewNet2.Mod(:,cl,:)=[];
    else
        NewNet2.Mod(:,:,Index)=[];
    end
    % Removal of the plasticity
    if isfield(NewNet2,'Plasticity') % If there are plasticity patterns in the genotype
        N_plast=size(NewNet2.Plasticity,2);
        for plast=1:N_plast
            affected=NewNet2.Plasticity{N_plast-plast+1}{3};
            if ismember(cl,affected)
                [NewNet2,transition_tracker]=remove_plasticity_determ(NewNet2,N_plast-plast+1,transition_tracker); % Removal of the related expression patterns is embedded into the plasticity removal
                
                if isstruct(transition_tracker)
                    transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)==N_plast-plast+1)=NaN; % "void" the removed plasticity
                    transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)>N_plast-plast+1)=transition_tracker.plasticities(2,transition_tracker.plasticities(2,:)>N_plast-plast+1)-1; % reduce the indices of the plasticities with higher index
                end
%                 removed.plast=[removed.plast N_plast-plast+1];
%                 removed.links={removed.links{:} report.links{:}};
%                 removed.expressions=[removed.expressions report.expressions];
%                 removed
%                 removed.links{:}
            else
                for class=1:size(affected,2)
                    if NewNet2.Plasticity{N_plast-plast+1}{3}(class)>cl % if the plasticity pattern refers to a subclass with bigger index than the removed subclass
                        NewNet2.Plasticity{N_plast-plast+1}{3}(class)=NewNet2.Plasticity{N_plast-plast+1}{3}(class)-1; % Reduce the index of this subclass in the reference
                    end
                end
            end
        end
    end
    % Removal of the expression patterns (unrelated to plasticity)
    if isfield(NewNet2,'Expression_patterns') % If there are any expressions
        N_exp=size(NewNet2.Expression_patterns,2);
        for exp=1:N_exp % for every pattern
            if NewNet2.Expression_patterns{N_exp-exp+1}{1}(1)==cl % If the pattern is expressed by the subclass that is being removed
                NewNet2=remove_expression_determ(NewNet2,N_exp-exp+1); % remove it
                if isstruct(transition_tracker)
                    transition_tracker.expressions(2,transition_tracker.expressions(2,:)==N_exp-exp+1)=NaN;
                    transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)=transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)-1;
                end
                %removed.expressions=[removed.expressions N_exp-exp+1];
            else % If the pattern is not expressed by the removed subclass, it can still refer to it in the patterns of other expressions
                N_link=size(NewNet2.Expression_patterns{N_exp-exp+1},2); % going from the last to the first to not have to adjust indexing
                for pat=2:N_link
                    if ismember(NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(1),[1 2 3 7]) % If the referred subclass of the link is in the third position
                        if NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)==cl % and it is the class that is being removed
                            NewNet2=remove_link_determ(NewNet2,N_exp-exp+1,N_link-pat+1); % remove the link
                            if isstruct(transition_tracker)
                                transition_exp_index=find(transition_tracker.expressions(2,:)==N_exp-exp+1);
                                transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)==N_link-pat+2)=NaN;
                                transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)=transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)-1;
                            end
                        elseif NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)>cl % If the referred class has a bigger index
                            NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)=NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)-1; % reduce the index
                        end
                    elseif ismember(NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(1),[8 9]) % If the referred subclass of the link is in the third or forth position
                        if ismember(cl,NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3:5)) % If one of the referred subclasses in the link matches the one being removed
                            NewNet2=remove_link_determ(NewNet2,N_exp-exp+1,N_link-pat+1); % remove the link
                            if isstruct(transition_tracker)
                                transition_exp_index=find(transition_tracker.expressions(2,:)==N_exp-exp+1);
                                transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)==N_link-pat+2)=NaN;
                                transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)=transition_tracker.links{transition_exp_index}(2,transition_tracker.links{transition_exp_index}(2,:)>N_link-pat+2)-1;
                            end
                        else % If there is no such match
                            if NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)>cl % If the referred class has a bigger index
                                NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)=NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(3)-1; % reduce the index
                            end
                            if NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)>cl % If the referred class has a bigger index
                                NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)=NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(4)-1; % reduce the index
                            end
                            if NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(5)>cl % If the referred class has a bigger index
                                NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(5)=NewNet2.Expression_patterns{N_exp-exp+1}{N_link-pat+2}{1}(5)-1; % reduce the index
                            end
                        end
                        
                    end
                    % The threshold and absolute refracterity do not "expand" beyond the class of expression, the plasticity-referring links were processed already.
                    
                    
                end
                if NewNet2.Expression_patterns{N_exp-exp+1}{1}(1)>cl % If the pattern is expressed by a subclass with a higher index than the removed one
                    NewNet2.Expression_patterns{N_exp-exp+1}{1}(1)=NewNet2.Expression_patterns{N_exp-exp+1}{1}(1)-1; % Reduce that index
                end
                if size(NewNet2.Expression_patterns{N_exp-exp+1},2)==1 % If there are no links left
                    NewNet2=remove_expression_determ(NewNet2,N_exp-exp+1); % remove the whole expression
                    if isstruct(transition_tracker)
                        transition_tracker.expressions(2,transition_tracker.expressions(2,:)==N_exp-exp+1)=NaN; 
                        transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)=transition_tracker.expressions(2,transition_tracker.expressions(2,:)>N_exp-exp+1)-1;
                    end
                    
                    
                    %removed.expressions=[removed.expressions N_exp-exp+1];
                end
            end
            
        end
    end
end

%legality(NewNet2)
end



