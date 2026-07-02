function Protocol=Protocol_maker_NumericalInputAssesment(seed)

if ~isnan(seed)
    rng(seed, 'twister');
end

% Here the input 1 is the synchro, and the input 2 consists of 10 neurons.
% The task is to distinguish the amount of simultaneously activated inputs.


Inputs=[1 10];
Output=1;  % Input-output composition of the simulated circuit.


Protocol.Primitives{1}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{1}{1}{2}=[[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0]];
%Protocol.Primitives{6}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; % Inverted
Protocol.Primitives{1}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 

Protocol.Primitives{2}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{2}{1}{2}=[[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0]];
%Protocol.Primitives{5}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0]; % Inverted
Protocol.Primitives{2}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0];

Protocol.Primitives{3}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{3}{1}{2}=[[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0]];
%Protocol.Primitives{4}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0]; % Inverted 
Protocol.Primitives{3}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0];

Protocol.Primitives{4}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{4}{1}{2}=[[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0]];
%Protocol.Primitives{3}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0]; % Inverted
Protocol.Primitives{4}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0];

Protocol.Primitives{5}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{5}{1}{2}=[[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 0];[0 0 0 0 0 0 0 0 0 0]];
%Protocol.Primitives{2}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0]; % Inverted
Protocol.Primitives{5}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0 0];

Protocol.Primitives{6}{1}{1}=[0 0 0 0 0 1 0 0 0 0]; % The inputs into three input channels
Protocol.Primitives{6}{1}{2}=[[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1];[0 0 0 0 0 0 0 0 0 1]];
%Protocol.Primitives{1}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1]; % Inverted
Protocol.Primitives{6}{2}{1}=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1];

Weight_template=[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

Delay_distrib=[20 40]; % The distance between the inputs as min and max of a uniform distribution of integers.
Tail=30; % The duration of the silent period in the end.
Repetitions=[1 1 1 1 1 1]; % How many times to apply each primitive. 
Order=[];
for n=1:size(Repetitions,2)
    Order=[Order repmat(n,1,Repetitions(n))];
end
Order=Order(randperm(size(Order,2)));

Delays=Delay_distrib(1)-1+randi(Delay_distrib(2)-Delay_distrib(1),1,sum(Repetitions));


Protocol.Input{1}=[];
Protocol.Input{2}=[];

Protocol.Output{1}=[];
Protocol.Weights{1}=[];
for n=1:sum(Repetitions)
    Protocol.Input{1}=[Protocol.Input{1} zeros(size(Protocol.Primitives{1}{1}{1},1),Delays(n)) Protocol.Primitives{Order(n)}{1}{1}];
    Protocol.Input{2}=[Protocol.Input{2} zeros(size(Protocol.Primitives{1}{1}{2},1),Delays(n)) Protocol.Primitives{Order(n)}{1}{2}];

    if n==1
        Protocol.Output{1}=[Protocol.Output{1} zeros(size(Protocol.Primitives{1}{2}{1},1),Delays(n)) Protocol.Primitives{Order(n)}{2}{1}];
        Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Delays(n)) Weight_template];
    else
        Protocol.Output{1}=[Protocol.Output{1} zeros(size(Protocol.Primitives{1}{2}{1},1),Delays(n)-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2))) Protocol.Primitives{Order(n)}{2}{1}];
        Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Delays(n)-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2))) Weight_template];
    end 
end

Protocol.Input{1}=[Protocol.Input{1} zeros(size(Protocol.Primitives{1}{1}{1},1),Tail)];
Protocol.Input{2}=[Protocol.Input{2} zeros(size(Protocol.Primitives{1}{1}{2},1),Tail)];

Protocol.Output{1}=[Protocol.Output{1} zeros(1,Tail-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2)))];
Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Tail-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2)))];

end