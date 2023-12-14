function [u]=obtain_u(npl,S_1,S_2,n1,ks,kc,kd);

u=cell(npl,1); % Cell for player's payoffs
% Initialize u
for i=1:npl
    aux=num2cell(2*ones(npl,1));
    u{i}=zeros(aux{:});
end
% Actions code: server -> 1 nd, 2 d; client-> 1 s, 2 ns

% Obtain payoffs

for i=1:2^npl %Encode actions as binary: each player only has 2!!
    action_string=dec2bin(i-1,npl);
    action=action_string-'0';
    action=action+1;
    action_aux=num2cell(action); %For indexing purposes
    % Obtain number of malicious clients
    nmc=sum(action(2:end)==1);
    % Obtain values of throughputs
    Sn_s=S_1(nmc+1); %Throughput in normal station due to malicious stations
    Sc_s=S_2(nmc+1); %Throughput in malicious station
    Sns=S_1(1); % Throughput if no malicious station is present
    % Obtain payoffs for the server
    if action(1)==1 %Server action: nd
        u{1}(action_aux{:})=ks*n1*(Sn_s-Sns);
    else
        u{1}(action_aux{:})=ks*n1*(Sns-Sn_s)-kd;
    end
    for j=2:npl
        if action(j)==2 %This client plays ns
            u{j}(action_aux{:})=0;
        else %This client plays s
            if action(1)==1 %Server does not detect
                u{j}(action_aux{:})=kc*(Sc_s-Sns);
            else %Server detects
                u{j}(action_aux{:})=-kc*Sns;
            end
        end
    end
end