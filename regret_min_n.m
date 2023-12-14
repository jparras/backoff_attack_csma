function[u,a_out]=regret_min_n(N,u_in,n_iter);

% RM for N players, when there is one of class 1 and N-1 of class 2
% N: total number of players
% u_1: Tensor with the payoffs for player in class 1
% u_2: Tensor with the payoffs for player in class 2
% n_iter: number of iterations of algorithm
% We consider each player has only TWO actions!!
na=2;
W=cell(N,1);
a=zeros(n_iter,N);
u=zeros(n_iter,N);
a_out=zeros(n_iter,na,N);
% Initialize W with the actions of each player set to 0
for i=1:N
    W{i}=zeros(na,1);
end

for t=1:n_iter
    % Obtain actions realizations for the current iteration
    for i=1:N %Obtain actions realizations for each player
        if max(W{i})<=0 %No regret, use random action
            a(t,i)=round(1 + (na-1).*rand(1));
        else
            aux=W{i};
            aux(aux<=0)=0;
            p_aux=aux/sum(aux);
            a_aux=mnrnd(1,p_aux);
            a(t,i)=find(a_aux==1);
        end
        a_out(t,a(t,i),i)=1;
    end
    % Update payoffs
    aux_idx=num2cell(a(t,:));
    for i=1:N
        u(t,i)=u_in{i}(aux_idx{:});
    end
    % Update regrets
    for i=1:N
        for j=1:na
            alt=a(t,:);
            alt(i)=j;
            aux_alt=num2cell(alt);
            W{i}(j)=W{i}(j)+u_in{i}(aux_alt{:})-u_in{i}(aux_idx{:});
        end
    end
end