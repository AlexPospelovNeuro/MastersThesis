% A general function for the offspring allocation. Replacement for the
% frac_offspring() function. Has its functions and more. 


function offsprings=offspring_allocation(winner_order,mode,varargin)

    switch mode
        case "proportional"
            if ~isempty(varargin)
                error('Proportional offspring allocation does not require extra arguments');
            end
            if var(winner_order(:,1))==0 % If all the circuits are exactly the same
                offsprings=1:size(winner_order,1);
            else
                winner_order(:,1)=winner_order(:,1)-winner_order(end,1);
                winner_order(:,1)=size(winner_order,1)*winner_order(:,1)/sum(winner_order(:,1));
                N_offspring(1:size(winner_order,1))=NaN;
                N_offspring(1)=round(winner_order(1,1));
                residual=winner_order(1,1)-N_offspring(1);
    
                for q=2:size(winner_order,1)
                    N_offspring(q)=round(winner_order(q,1)+residual);
                    residual=winner_order(q,1)-N_offspring(q)+residual;
                end
                offsprings=[];
                for q=1:size(winner_order,1)
                    if N_offspring(q)>0
                        offsprings=[offsprings repmat(winner_order(q,2),1,N_offspring(q))];
                    end
        
                end
            end

        case "winner"
            if numel(varargin) > 1
                error('Winner-favoring offspring allocation requires one or now extra arguments');
            elseif isempty(varargin)
                winners_n={1};
            else
                winners_n=varargin;
            end
            if var(winner_order(:,1))==0 % If all the circuits are exactly the same
                offsprings=1:size(winner_order,1);
            else
                offspring_per_winner=numel(winner_order(:,1))/winners_n{1}; % May be noninteger
                numbers_c=repmat(offspring_per_winner,winners_n{1},1);
                numbers_n=cont2int(numbers_c); % Has to be integer
                numbers_n=numbers_n(randperm(numel(numbers_n))); % To remove unwanted preference amoung the winners
                offsprings=[];
                for q=1:winners_n{1}
                    offsprings=[offsprings repmat(winner_order(q,2),1,numbers_n(q))];

                end
            end

        otherwise
            error('Unknown offspring allocation mode')
    end


end