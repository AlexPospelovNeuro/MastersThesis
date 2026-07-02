% A function to make a truncated continuous normal distribution and
% allocate the values 
% INPUTS:
% N - the number of the elements (presumably - cells)
% mu - the mean of the distribution
% sigma - the sigma of the distribution
% cutoff - the threshold value below which the values cannot go (there is no reasons to assume existense of the upper cutoff)
% If cutoff=NaN, do not truncate
% OUTPUTS:
% values - N-by-one vector of real values distributed according to the
% truncated normal distribution. Order of values is randomized


function values=trunc_cont_norm_dist(N,mu,sigma,cutoff)

if sigma==0
    if ~isnan(cutoff) && cutoff>mu
        values(1,1:N)=cutoff;
    else
        values(1,1:N)=mu;
    end
else

    precision=1000/sigma;
    
 
    continuous=mu-ceil(log10(N)+2)*sigma:1/precision:mu+ceil(log10(N)+2)*sigma;
    y = cdf('Normal',continuous,mu,sigma);
    
    
    
    %find(continuous>=cutoff,1)
    %y(find(continuous>=cutoff,1))
    if ~isnan(cutoff)
        thresh=find(continuous>=cutoff,1);
        if isempty(thresh) 
            y(:)=0;
        else
            y=(y-y(thresh))/(1-y(thresh));
            y(continuous<cutoff)=0;
        end
    end
    y=N*y;
    pointers=0.5:1:N-0.5;

    values(1:N)=0;
    for C=1:N
        if ~isempty(find(y>pointers(C),1))
            values(C)=continuous(find(y>pointers(C),1));
        end
    end
    values=values(randperm(N));
end

end