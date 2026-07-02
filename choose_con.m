% This function takes in the vector of affinities (possibly - just numbers
% of connections), assuming that it is sorted in in descending order
% already, and the number of connections to make. If the number of
% connections is bigger or equal to the total number of neurons to connect
% to, the functions returns 1:n sequence ("connect to every available neuron")
% If there are more cells available than connections required, the function
% returns an increasing sequence of N numbers which points to the cells to
% connect to chosen based in ther affinity based on the equal affinity
% distribution algorithm




% N=3;
% affin=randi(15,[30 1]);
% affin=sort(affin,'descend');

function pointers=choose_con(N,affin)
if min(affin)<0
    affin=affin-min(affin);
end
if N>=size(affin,1)
    pointers=1:size(affin,1);
else
    Total=sum(affin);
    Step=Total/N;
    C_affin=cumsum(affin);

    pointers(1:N)=NaN;
    for t=1:N
        pointers(t)=find(C_affin>(t-1)*Step,1);
    end
    for t=2:N
        if pointers(t)<=pointers(t-1)
            pointers(t)=pointers(t-1)+1;
        end
    end


end