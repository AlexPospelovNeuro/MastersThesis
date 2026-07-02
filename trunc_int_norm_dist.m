% The function that produces the truncation of the normal distribution
% and allocates the integer values. 
% INPUTS:
% N - the number of the elements (presumably - cells)
% mu - the mean of the distribution
% sigma - the sigma of the distribution
% cutoff - the threshold value below which the values cannot go (there is no reasons to assume existense of the upper cutoff)
% If cutoff=NaN, do not truncate
% OUTPUTS:
% values - N-by-one vector of integer values distributed according to the
% rectified normal distribution. Order of values is randomized

function values=trunc_int_norm_dist(N,mu,sigma,cutoff)

if sigma==0
    values1(1,1:N)=mu;
    values(1,1:N)=0; % generally speaking, non-integer
    % making it integer
    residual=0;
    for n=1:N
        values(n)=floor(values1(n)+residual);
        residual=values1(n)+residual-values(n);
    end
elseif N==0
    values=[];
else
    
    precision=1000;
    continuous=mu-ceil(log10(N)+2)*sigma:1/precision:mu+ceil(log10(N)+2)*sigma;
%     mu
%     sigma
    y = pdf('Normal',continuous,mu,sigma);
    
    if ~isnan(cutoff)
        y(continuous<cutoff)=0;
    end
    Integers=floor(mu-3*sigma):ceil(mu+3*sigma);
    discrete(1:length(Integers))=NaN;
    for I=1:length(Integers)
        discrete(I)=sum(y((continuous>=Integers(I)-0.5)&(continuous<=Integers(I)+0.5)));
    end
    if ~isnan(cutoff)
        discrete(Integers==cutoff)=discrete(Integers==cutoff)*2;
    end
    if sum(discrete)~=0
        discrete=N*discrete/sum(discrete);
    end
%discrete=discrete*N;

%discrete1(1:length(Integers))=0;
% residual=0;
% count=0;
% if isempty(find(Integers==0, 1))
%     start=1;
% else 
%     start=find(Integers==0);
% end
% for p=start:length(Integers)
%     Int=floor(discrete(p)+residual);
%     residual=discrete(p)+residual-Int;
%     discrete1(p)=Int;
%     count=count+Int;
%     if count>=N
%         discrete1(p)=discrete1(p)-count+N;
%         break
% 
%     end
% 
% end
    discrete1=cont2int(discrete);
% figure
% plot(Integers,discrete,'*')
% hold on
% plot(Integers,discrete1,'*r')

    values(1,1:N)=0;
    pointer=1;
    for n=1:length(discrete1)
        values(pointer:pointer+discrete1(n)-1)=Integers(n);
        pointer=pointer+discrete1(n);
    end
    
end
values=values(randperm(N));
end