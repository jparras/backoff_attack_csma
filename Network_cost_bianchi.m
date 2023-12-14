%% BIACHI MODEL SIMULAIONS - LONG PAYLOAD
% Juan Parras, GAPS - UPM, May 2017
clear all; clc; close all

%% Input parameters
PAYLOAD=8184; % Payload bits
MACh=272; %MAC header bits
PHYh=128; %PHY header
ACK=112+PHYh; %ACK packet length in bits
RTS=160+PHYh; %RTS packet length in bits
CTS=112+PHYh; %CTS packet length in bits
Rb=1; % Bit rate (in Mbps)
pd=1; %Propagation delay in us
Ts=50; %Slot time in us
SIFS=28; %SIFS time in us
DIFS=128; % DIFS time in us
ACK_timeout=300; %Timeout in us
CTS_timeout=300; %TImeout in us
Tp=PAYLOAD/Rb;

CW_min=31;
CW_max=1023;
m=log2((CW_max+1)/(CW_min+1));

% Times in basic access (ONLY)
Tt=MACh/Rb+PHYh/Rb+Tp+SIFS+pd+ACK/Rb+DIFS+pd;
Tc=MACh/Rb+PHYh/Rb+Tp+DIFS+pd;

sa=1; %Paramter to save data

N=5; %Total number of players
W=8; % W2 of malicious node

% Output variables
nc=5; %Number of cases:
% 1 - n1=5, n2=0
% 2 - n1=4, n2=1
% 3 - n1=3, n2=2
% 4 - n1=2, n2=3
% 5 - n1=1, n2=4
S_1=zeros(1,nc); %Throughput normal station
S_2=zeros(1,nc); %Throughput malicious station
S=zeros(1,nc); %Throughput aggregated

n1_v=[5 4 3 2 1];
n2_v=[0 1 2 3 4];

%% Case 1: only 1 bad player

func_p=@(pi,n) 1-(1-pi)^(n-1);
func_pi=@(p,m,Wm) 2/(1+Wm+p.*Wm.*sum((2*p).^(0:m-1)));
%func_pi=@(p,m,Wm) 2*(1-2*p)/((1-2*p)*(Wm+1)+(p*Wm*(1-(2*p)^m)));

syms ps pis

display('Case 1');
sol=vpasolve([func_p(pis,N)==ps,func_pi(ps,m,CW_min)==pis],[ps,pis],[0,1;0,1]);
p=vpa(sol.ps);
pi=vpa(sol.pis);
Ptr=1-(1-pi)^N;
Ps=N*pi*(1-pi)^(N-1)/Ptr;
S(1)=(Ps*Ptr*Tp)/(Ps*Ptr*Tt+Ptr*(1-Ps)*Tc+(1-Ptr)*Ts);
S_1(1)=S(1)/N;
S_2(1)=S(1)/N;

%% Numerical solving cases 2-5

func_pi1=@(p,m,Wm) 2/(1+Wm+p.*Wm.*sum((2*p).^(0:m-1)));
func_pi2=@(W) 2/(W+1);
func_p1=@(pi1,pi2,n1,n2) 1-(1-pi1)^(n1-1)*(1-pi2)^n2;
func_p2=@(pi1,pi2,n1,n2) 1-(1-pi2)^(n2-1)*(1-pi1)^n1;
syms p1s pi1s p2s pi2s

for iw=2:nc
    n1=n1_v(iw);
    n2=n2_v(iw);
    sol=vpasolve([func_pi1(p1s,m,CW_min)==pi1s,func_pi2(W)==pi2s,func_p1(pi1s,pi2s,n1,n2)==p1s,func_p2(pi1s,pi2s,n1,n2)==p2s],[p1s pi1s p2s pi2s],[0,1;0,1;0,1;0,1]);
    p1=vpa(sol.p1s);
    pi1=vpa(sol.pi1s);
    p2=vpa(sol.p2s);
    pi2=vpa(sol.pi2s);
    if (isempty(p1) || isempty(pi1) || isempty(p2) || isempty(p2s))==0
        Ptr=1-((1-pi1)^n1)*((1-pi2)^n2);
        Ps1=pi1*(1-pi1)^(n1-1)*(1-pi2)^n2;
        Ps2=pi2*(1-pi2)^(n2-1)*(1-pi1)^n1;
        Es=((n1*Ps1+n2*Ps2)*Tt+(Ptr-(n1*Ps1+n2*Ps2))*Tc+(1-Ptr)*Ts);
        S_1(iw)=Ps1*Tp/Es;
        S_2(iw)=Ps2*Tp/Es;
        S(iw)=n1*S_1(iw)+n2*S_2(iw);
    end
end
%% Plot
figure();
plot(1:nc,S,'bo-',1:nc,S_1,'kx-',1:nc,S_2,'rs-');
grid on;
legend('S','S_1','S_2');


%% Save
if sa
    save('Data_network_cost_simulations','S','S_1','S_2','n1_v','n2_v')
end