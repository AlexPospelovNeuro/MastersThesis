% This is the function that makes the pairing of connections in the ordered
% pair of neuron subclasses. Visualization included. 

% INPUTS:
%  [gauss] - a vector with genetically coded properties of random
% connection, in the format of: [Post_N Post_mu Post_sigma Pre_N Pre_sigma]
%  {Post_connection_pattern} - a 1xk cell array, where k is number of
%  patterns that affect the number of connections from the postsynaptic perspective. Each pattern has the from
%  of {ID,[values]}, where ID is the ID of the pattern and the values are
%  calculated from the [a b c] coefficients (elsewhere)
%  {Pre_connection_pattern} - a 1xk cell array, where k is number of
%  patterns that affect the number of connections from the presynaptic perspective. Each pattern has the from
%  of {ID,[values]}, where ID is the ID of the pattern and the values are
%  calculated from the [a b c] coefficients (elsewhere)
%  {Post_preference_pattern} - a 1xk cell array, where k is number of
%  patterns that affect the connections preference from the postsynaptic perspective. Each pattern has the from
%  of {ID,[values]}, where ID is the ID of the pattern and the values are
%  calculated from the [a b c] coefficients (elsewhere)
%  {Pre_preference_pattern} - a 1xk cell array, where k is number of
%  patterns that affect the connections preference from the presynaptic perspective. Each pattern has the from
%  of {ID,[values]}, where ID is the ID of the pattern and the values are
%  calculated from the [a b c] coefficients (elsewhere)

% The patterns are continuous, not integer

% OUTPUT:
% Post_N x Pre_N connection matrix



function connections=make_connections_v2(gauss,Post_connection_pattern,Pre_connection_pattern,Post_affinity_pattern,Pre_affinity_pattern)





if gauss(1)==0
    connections=[];
else
connections(1:gauss(1),1:gauss(4))=0;
Post_values=trunc_int_norm_dist(gauss(1),gauss(2),gauss(3),0); % Generate the gaussian-derived values
for pat=1:size(Post_connection_pattern,1)
    Post_values=Post_values+Post_connection_pattern{pat}; % Add the connection patterns for postsynapse;
end
PostV_continuous=Post_values; % Non-rounded values are saved to use for the affinity calculation (since affinity is continuous)
Post_values=cont2int(Post_values); % After all the numbers are summed up, apply rounding
Post_values(Post_values<0)=0; % remove negative values
Total_connections=sum(Post_values); % calculate total number of connections


Pre_values(1,1:gauss(4))=0;    % Number of presynaptic conections is zero at this point
for pat=1:size(Pre_connection_pattern,1)
    Pre_values=Pre_values+Pre_connection_pattern{pat}; % Only the pattern-derived connections are considered;
end
Pre_gauss_connection_N=Total_connections-sum(Pre_values); % Calculate the number of gaussian-derived connections by subtracting the pattern-derived from the total number
if Pre_gauss_connection_N>0
    Pre_gauss_values=trunc_int_norm_dist(gauss(4),Pre_gauss_connection_N/gauss(4),gauss(5),0); % calculate the numbers of gaussian-derived presynaptic connections per cell
else
    Pre_gauss_values=zeros([1 gauss(4)]);
end
Pre_values=Pre_gauss_values+Pre_values; % Add up pattern-derived and gaussian-derived connections
PreV_continuous=Pre_values; % Non-rounded values are saved to use for the affinity calculation (since affinity is continuous)
Pre_values=cont2int(Pre_values);
Pre_values(Pre_values<0)=0; % remove negative values

pre_ID=1:gauss(4);
pre_table=[pre_ID' Pre_values' Pre_values' PreV_continuous'];
for pat=1:size(Pre_affinity_pattern,1)
    pre_table(:,4)=pre_table(:,4)+Pre_affinity_pattern{pat}'; % Add the affinity patterns presynapse;
end
pre_table_sorted=sortrows(pre_table,4,"descend");
pre_table_sorted(:,5)=pre_table_sorted(:,4)./pre_table_sorted(:,3);  % impact of each synapse of each neuron to the affinity of this neuron


post_ID=1:gauss(1);

post_table=[post_ID' Post_values' Post_values' PostV_continuous'];  % The second column is for the record. The third is for calculations. The fours is for sorting.
for pat=1:size(Post_affinity_pattern,1)
    post_table(:,4)=post_table(:,4)+Post_affinity_pattern{pat}'; % Add the connection patterns presynapse;
end

post_table_sorted=sortrows(post_table,4,"descend");


for n=1:gauss(1) % for each postsynaptic neuron
    if post_table_sorted(n,2)>0 % If it has any input connections
        sources=find(pre_table_sorted(:,3)>0); % The number of potential source neurons for the given post synaptic neuron 
        synapses=sources(choose_con(post_table_sorted(n,3),pre_table_sorted(sources,4))); % get the list of connections
        connections(post_table_sorted(n,1),pre_table_sorted(synapses,1))=1; % update the connection matrix
        pre_table_sorted(synapses,3)=pre_table_sorted(synapses,3)-1; % Reduce the number of available connections for the next neuron
        pre_table_sorted(synapses,4)=pre_table_sorted(synapses,4)-pre_table_sorted(synapses,5); % Reduce the affinity for the sorting before the next neuron
        clear sources synapses
    end
    pre_table_sorted=sortrows(pre_table_sorted,4,"descend");   % Redo the sorting based on the updated affinity, I don't know why it is necessary but it is
end
    
end

end