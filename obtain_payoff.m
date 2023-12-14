function[payoff]=obtain_payoff(npl,u,actions)

payoff=zeros(1,npl);
for k=1:2^npl %Encode actions as binary: each player only has 2!!
    action_string=dec2bin(k-1,npl);
    action=action_string-'0';
    action=action+1;
    action_aux=num2cell(action); %For indexing purposes
    pr=zeros(1,npl);
    for j=1:npl %Update probabilities
        if action(j)==1
            pr(j)=actions(j);
        else
            pr(j)=1-actions(j);
        end
    end
    pr=prod(pr);
    for j=1:npl %Payoff per player
        payoff(j)=payoff(j)+u{j}(action_aux{:})*pr;
    end
end