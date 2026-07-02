% A small function to create a "desired_output" variable within the parfor
% loop

function [Desired_output, weights]=make_desired_output(Protocol,tech)
    for out=1:size(Protocol.Output,2)
        Desired_output{out}=[zeros(size(Protocol.Output{out},1),tech.Longest_shape) Protocol.Output{out}]; 
        weights{out}=[zeros(size(Protocol.Weights{out},1),tech.Longest_shape) Protocol.Weights{out}];
    end
end


