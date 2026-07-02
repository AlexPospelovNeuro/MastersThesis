
% 6.5.25 probability is added as a parameter instead of being hardcoded inside the function. More handy.  

function Protocol=Protocol_maker_prob(seed,p)

if ~isnan(seed)
    rng(seed, 'twister');
end

Inputs=[1 1 1];
Output=1;  % Input-output composition of the simulated circuit.

%p=0.2; % The probability that the reaction to the stimulus is "correct" (and also that the absence of reaction to the absence of stimulus is incorrect. )



% Correct presence of the response to the input
Protocol.Primitives{1}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{1}{1}{2}=[0 0 0 0 0 1 0 0 0 0];
Protocol.Primitives{1}{1}{3}=[0 0 0 0 0 0 0 0 0 0];

%Protocol.Primitives{1}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % The expected output from the single output channels
Protocol.Primitives{1}{2}{1}=[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0];


% Correct absence of the response to the absence of  input
Protocol.Primitives{2}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{2}{1}{2}=[0 0 0 0 0 0 0 0 0 0];
Protocol.Primitives{2}{1}{3}=[0 0 0 0 0 0 0 0 0 0];

Protocol.Primitives{2}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % The expected output from the single output channels


% Correct absence of the response to the input
Protocol.Primitives{3}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{3}{1}{2}=[0 0 0 0 0 1 0 0 0 0];
Protocol.Primitives{3}{1}{3}=[0 0 0 0 0 0 0 0 0 0];

Protocol.Primitives{3}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % The expected output from the single output channels



% Correct presence of the response to the  absence input
Protocol.Primitives{4}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{4}{1}{2}=[0 0 0 0 0 0 0 0 0 0];
Protocol.Primitives{4}{1}{3}=[0 0 0 0 0 0 0 0 0 0];

%Protocol.Primitives{4}{2}{1}=[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0]; % The expected output from the single output channels
Protocol.Primitives{4}{2}{1}=[0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0];

Weight_template=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

Delay_distrib=[20 40]; % The distance between the inputs as min and max of a uniform distribution of integers.
Tail=30; % The duration of the silent period in the end.
Repetitions=round([20*p 20*p 20*(1-p) 20*(1-p)]); % How many times to apply each primitive. 
Order=[];
for n=1:size(Repetitions,2)
    Order=[Order repmat(n,1,Repetitions(n))];
end
Order=Order(randperm(size(Order,2)));

Delays=Delay_distrib(1)-1+randi(Delay_distrib(2)-Delay_distrib(1),1,sum(Repetitions));

Protocol.Input{1}=[];
Protocol.Input{2}=[];
Protocol.Input{3}=[];
Protocol.Output{1}=[];
Protocol.Weights{1}=[];
for n=1:sum(Repetitions)
    Protocol.Input{1}=[Protocol.Input{1} zeros(1,Delays(n)) Protocol.Primitives{Order(n)}{1}{1}];
    Protocol.Input{2}=[Protocol.Input{2} zeros(1,Delays(n)) Protocol.Primitives{Order(n)}{1}{2}];
    Protocol.Input{3}=[Protocol.Input{3} zeros(1,Delays(n)) Protocol.Primitives{Order(n)}{1}{3}];
    if n==1
        Protocol.Output{1}=[Protocol.Output{1} zeros(1,Delays(n)) Protocol.Primitives{Order(n)}{2}{1}];
        Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Delays(n)) Weight_template];
    else
        Protocol.Output{1}=[Protocol.Output{1} zeros(1,Delays(n)-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2))) Protocol.Primitives{Order(n)}{2}{1}];
        Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Delays(n)-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2))) Weight_template];
    end 
end

Protocol.Input{1}=[Protocol.Input{1} zeros(1,Tail)];
Protocol.Input{2}=[Protocol.Input{2} zeros(1,Tail)];
Protocol.Input{3}=[Protocol.Input{3} zeros(1,Tail)];
Protocol.Output{1}=[Protocol.Output{1} zeros(1,Tail-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2)))];
Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Tail-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2)))];

end