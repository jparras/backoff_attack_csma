%% STATIC GAME LEARNING USING REGRET-MATCHING
% Juan Parras, GAPS-UPM, May 2017
clear all; clc; close all;
%% Initial parameters

n_iter=2e3; %Number of iterations
n_avg=50; %Number of averages per equilibrium

%Save parameter
sa=0;

%Plot parameter
pl=1;

% Load prestored data
load('Data_network_cost_simulations');
N=5; %Number of players
% Cost parameters
ks=1;
kc=1;
kd=0.1;

%% Case 1
display('Case 1')
n2=1;
n1=N-n2; 

npl=1+n2; %Number of players

u=obtain_u(npl,S_1,S_2,n1,ks,kc,kd);

am=-u{1}(1,1);
ac=u{1}(2,1);
af=-u{1}(2,2);
bs=u{2}(1,1);
bc=-u{2}(2,1);

%Theoretical values

yn=bc/(bc+bs);
zn=af/(af+ac+am);
v1n=-af*am/(af+am+ac);
v2n=0;

display('THEORETICAL VALUES');
display(['yn = ' num2str(yn) ' zn = ' num2str(zn)]);
display(['u1 = ' num2str(v1n) ' u2 = ' num2str(v2n)]);

y_avg=zeros(n_avg,1);
z_avg=zeros(n_avg,1);
err_y_avg=zeros(n_avg,1);
err_z_avg=zeros(n_avg,1);
u1=zeros(n_avg,1);
u2=zeros(n_avg,1);
err_u1=zeros(n_avg,1);
err_u2=zeros(n_avg,1);

for i=1:n_avg
    %RM
    [u_out,a]=regret_min_n(npl,u,n_iter);
    % Data save
    y_avg(i)=mean(a(:,1,1));
    z_avg(i)=mean(a(:,1,2));
    err_y_avg(i)=y_avg(i)-yn;
    err_z_avg(i)=z_avg(i)-zn;
    u1(i)=mean(u_out(:,1));
    u2(i)=mean(u_out(:,2));
    err_u1(i)=u1(i)-v1n;
    err_u2(i)=u2(i)-v2n;
end

[mean(err_y_avg) std(err_y_avg) mean(err_z_avg) std(err_z_avg)]
[mean(err_u1) std(err_u1) mean(err_u2) std(err_u2)]

if pl==1
    figure();
    subplot(2,2,1); hist(err_y_avg);title('y_n'); grid on;
    subplot(2,2,2); hist(err_z_avg);title('z_n'); grid on;
    subplot(2,2,3); hist(err_u1);title('u_1');grid on;
    subplot(2,2,4); hist(err_u2);title('u_2');grid on;

    figure();
    plot(1:n_iter,squeeze(cumsum(a(:,1,:)))./repmat((1:n_iter)',1,npl));
    grid on;
end

% [yn mean(y_avg) std(y_avg)]
% [zn mean(z_avg) std(z_avg)]
% [v1n mean(u1) std(u1)]
% [v2n mean(u2) std(u2)]

%% Case 2
actions_avg_cell=cell(3,1);
u_avg_cell=cell(3,1);
idcell=1;
for n2=[2 3 4]
    display(['Case n2 = ' num2str(n2)]);
    n1=N-n2; 

    npl=1+n2; %Number of players
    
    u=obtain_u(npl,S_1,S_2,n1,ks,kc,kd);
    
    actions_avg=zeros(n_avg,npl);
    u_avg=zeros(n_avg,npl);

    for i=1:n_avg
        %RM
        [u_out,a]=regret_min_n(npl,u,n_iter);
        % Data save
        u_avg(i,:)=mean(u_out);
        a_aux=squeeze(mean(a));
        actions_avg(i,:)=a_aux(1,:);
    end
    
    actions_avg_cell{idcell}=actions_avg;
    u_avg_cell{idcell}=u_avg;
    idcell=idcell+1;
    
    if pl==1
        figure();
        plot(1:n_iter,squeeze(cumsum(a(:,1,:)))./repmat((1:n_iter)',1,npl));
        grid on;
        legend('Server','Client 1','Client 2','Client 3','Client 4');
    end
end
%% Values output
prec=4; %Number of precission decimals
display('RM RESULTS');
for i=1:4
    display(['Case n2 = ' num2str(i)]);
    if i==1
        display(['Actions = ' mat2str([mean(y_avg) mean(z_avg)],prec)]);
        display(['Payoffs = ' mat2str([mean(u1) mean(u2)],prec)]);
    else
        display(['Actions = ' mat2str(mean(actions_avg_cell{i-1}),prec)]);
        display(['Payoffs = ' mat2str(mean(u_avg_cell{i-1}),prec)]);
    end
end
%% Values analysis
% We analyze both the mean values, as well as the values in a realization
display('RESULTS');
for i=1:4
    figure();
    display(['Case n2 = ' num2str(i)]);
    npl=i+1;
    h=[];
    if i==1
        actions=[mean(y_avg) mean(z_avg)];
        actions_r=[y_avg(end) z_avg(end)];
        payoff_emp=[mean(u1) mean(u2)];
        h=[y_avg z_avg];
        hist(h,linspace(0.05,0.25,3));
    else
        actions=mean(actions_avg_cell{i-1});
        actions_r=actions_avg_cell{i-1}(end,:);
        payoff_emp=mean(u_avg_cell{i-1});
        h=actions_avg_cell{i-1};
        hist(h,5);
    end
    grid on;
    legend('Server','Client 1','Client 2','Client 3','Client 4');
    u=obtain_u(npl,S_1,S_2,N-i,ks,kc,kd);
    payoff=obtain_payoff(npl,u,actions);
    payoff_r=obtain_payoff(npl,u,actions_r);
    display(['Actions = ' mat2str(actions,prec)]);
    display(['Payoffs emp = ' mat2str(payoff_emp,prec)]);
    display(['Payoff th = ' mat2str(payoff,prec)]);
    display(['Actions r = ' mat2str(actions_r,prec)]);
    display(['Payoff r = ' mat2str(payoff_r,prec)]);
end
%% Save
if sa==1
    save('Values_RM_learning_paper');
end