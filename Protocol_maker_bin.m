% A first protocol constructor for the Net2 type of circuits. To be heavily
% expanded to accomodate the feedback protocols. Later. 
% Just uncomment the %Protocol.Primitives{_}{2}{1} lines that labelled with
% the function you want to model, and comment all the rest. 


function Protocol=Protocol_maker_bin(seed)

if ~isnan(seed)
    rng(seed, 'twister');
end
Inputs=[1 1 1];
Output=1;  % Input-output composition of the simulated circuit.


Protocol.Primitives{1}{1}{1}=logical([0 0 0 0 0 1 0 0 0 0]); % The inputs into three input channels
Protocol.Primitives{1}{1}{2}=logical([0 0 0 0 0 0 0 0 0 0]);
Protocol.Primitives{1}{1}{3}=logical([0 0 0 0 0 0 0 0 0 0]);
%Protocol.Primitives{1}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]); % OR AND XOR
Protocol.Primitives{1}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0]); % NOR NAND NXOR



Protocol.Primitives{2}{1}{1}=logical([0 0 0 0 0 1 0 0 0 0]); % The inputs into three input channels
Protocol.Primitives{2}{1}{2}=logical([0 0 0 0 0 1 0 0 0 0]);
Protocol.Primitives{2}{1}{3}=logical([0 0 0 0 0 0 0 0 0 0]); 
%Protocol.Primitives{2}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]); % NOR NXOR AND
Protocol.Primitives{2}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0]); % NAND OR XOR

Protocol.Primitives{3}{1}{1}=logical([0 0 0 0 0 1 0 0 0 0]); % The inputs into three input channels
Protocol.Primitives{3}{1}{2}=logical([0 0 0 0 0 0 0 0 0 0]);
Protocol.Primitives{3}{1}{3}=logical([0 0 0 0 0 1 0 0 0 0]); 
%Protocol.Primitives{3}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]); % NOR NXOR AND
Protocol.Primitives{3}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0]); % NAND OR XOR

Protocol.Primitives{4}{1}{1}=logical([0 0 0 0 0 1 0 0 0 0]); % The inputs into three input channels
Protocol.Primitives{4}{1}{2}=logical([0 0 0 0 0 1 0 0 0 0]);
Protocol.Primitives{4}{1}{3}=logical([0 0 0 0 0 1 0 0 0 0]);
Protocol.Primitives{4}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]); % XOR NAND NOR
%Protocol.Primitives{4}{2}{1}=logical([0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0]); % NXOR AND OR


Weight_template=[0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

Delay_distrib=[20 40]; % The distance between the inputs as min and max of a uniform distribution of integers.
Tail=30; % The duration of the silent period in the end.
Repetitions=[2 2 2 2]; % How many times to apply each primitive. 
Order=[];
for n=1:size(Repetitions,2)
    Order=[Order repmat(n,1,Repetitions(n))];
end
Order=Order(randperm(size(Order,2)));

Delays=Delay_distrib(1)-1+randi(sum(Repetitions),1,Delay_distrib(2)-Delay_distrib(1));

Protocol.Input{1}=[];
Protocol.Input{2}=[];
Protocol.Input{3}=[];
Protocol.Output{1}=[];
Protocol.Weights{1}=[];
for n=1:sum(Repetitions)
    Protocol.Input{1}=[Protocol.Input{1} false(1,Delays(n)) Protocol.Primitives{Order(n)}{1}{1}];
    Protocol.Input{2}=[Protocol.Input{2} false(1,Delays(n)) Protocol.Primitives{Order(n)}{1}{2}];
    Protocol.Input{3}=[Protocol.Input{3} false(1,Delays(n)) Protocol.Primitives{Order(n)}{1}{3}];
    if n==1
        Protocol.Output{1}=[Protocol.Output{1} false(1,Delays(n)) Protocol.Primitives{Order(n)}{2}{1}];
        Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Delays(n)) Weight_template];
    else
        Protocol.Output{1}=[Protocol.Output{1} false(1,Delays(n)-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2))) Protocol.Primitives{Order(n)}{2}{1}];
        Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Delays(n)-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2))) Weight_template];
    end 
end

Protocol.Input{1}=[Protocol.Input{1} false(1,Tail)];
Protocol.Input{2}=[Protocol.Input{2} false(1,Tail)];
Protocol.Input{3}=[Protocol.Input{3} false(1,Tail)];
Protocol.Output{1}=[Protocol.Output{1} false(1,Tail-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2)))];
Protocol.Weights{1}=[Protocol.Weights{1} zeros(1,Tail-(size(Protocol.Primitives{Order(n-1)}{2}{1},2)-size(Protocol.Primitives{Order(n-1)}{1}{1},2)))];

end








