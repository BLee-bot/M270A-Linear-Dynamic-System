%% Setting
clear all
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;
ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;
% figure(1)
% plot(t,u1); hold on;
% plot(t,u2);

%% Task 4-1

Mean_u1=sum(u1)/length(u1); % Average Input of u1 is almost zero
Mean_u2=sum(u2)/length(u2); % Average Input of u2 is almost zero

%% Task 4-2
% Estimate
u=[u1;u2];
lag=-5:0.025:5;
P=10000;
K=200;
Ruu11=[];
Ruu21=[];
Ruu12=[];
Ruu22=[];

for k=-K:1:K
    Ruu=0;
for q=-P:1:P 
    Ruu=Ruu+(1/(2*P))*(u(1:2,12000+k+q)*u(1:2,12000+q)');
end
    Ruu11=[Ruu11 Ruu(1,1)];
    Ruu21=[Ruu21 Ruu(2,1)];
    Ruu12=[Ruu12 Ruu(1,2)];
    Ruu22=[Ruu22 Ruu(2,2)];
end

% figure(1)
% subplot(221)
% plot(lag,Ruu11,'r'); hold on;
% axis([-5 5 -1 5])
% ylabel('$Ruu$ (1,1)','FontSize',14,'Interpreter','Latex');
% xlabel('lag factor','FontSize',14,'Interpreter','Latex');
% grid on;
% 
% subplot(222)
% plot(lag,Ruu21,'g');
% axis([-5 5 -1 5])
% ylabel('$Ruu$ (1,2)','FontSize',14,'Interpreter','Latex');
% xlabel('lag factor','FontSize',14,'Interpreter','Latex');
% grid on;
% 
% subplot(223)
% plot(lag,Ruu12,'b');
% axis([-5 5 -1 5])
% ylabel('$Ruu$ (2,1)','FontSize',14,'Interpreter','Latex');
% xlabel('lag factor','FontSize',14,'Interpreter','Latex');
% grid on;
% 
% subplot(224)
% plot(lag,Ruu22,'k'); hold off;
% axis([-5 5 -1 5])
% ylabel('$Ruu$ (2,2)','FontSize',14,'Interpreter','Latex');
% xlabel('lag factor','FontSize',14,'Interpreter','Latex');
% grid on;
% sgtitle('Ruu at P=10000')
%% Task 4-3

Ruu0=[Ruu11(1,201) Ruu12(1,201);Ruu21(1,201) Ruu22(1,201)];


%% Task 4-4

%Estimate
u=[u1;u2];
y=[y1;y2];
lag2=-0.2:0.025:2;
P2=10000;
K2=80;

Ryu11=[];
Ryu21=[];
Ryu12=[];
Ryu22=[];


for k=-8:1:K2
    Ryu=0;
for q=-P2:1:P2 
    Ryu=Ryu+(1/(2*P))*(y(1:2,12000+k+q)*u(1:2,12000+q)');
end
    Ryu11=[Ryu11 Ryu(1,1)];
    Ryu21=[Ryu21 Ryu(2,1)];
    Ryu12=[Ryu12 Ryu(1,2)];
    Ryu22=[Ryu22 Ryu(2,2)];
end
%%
load u1_impulse.mat
y11 = u1_impulse.Y(3).Data;
y21 = u1_impulse.Y(4).Data;
u1 = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi] = max(u1>0); %%% find index where pulse occurs
load u2_impulse.mat
y12 = u2_impulse.Y(3).Data;
y22 = u2_impulse.Y(4).Data;
u2 = u2_impulse.Y(2).Data;
%%% remove any offsets in output data using data prior to pulse application
y11 = y11 - mean(y11([1:mi-1]));
y12 = y12 - mean(y12([1:mi-1]));
y21 = y21 - mean(y21([1:mi-1]));
y22 = y22 - mean(y22([1:mi-1]));
%%% rescale IO data so that impulse input has magnitude 1
y11 = y11/max(u1);
y12 = y12/max(u2);
y21 = y21/max(u1);
y22 = y22/max(u2);
u1 = u1/max(u1);
u2 = u2/max(u2);
ts = 1/40; %%%% sample period
N = length(y11); %%%% length of data sets
t = [0:N-1]*ts - 1;

%%
tt=t(1,33:121);
y11=y11(1,33:121);
y21=y21(1,33:121);
y12=y12(1,33:121);
y22=y22(1,33:121);

figure(1)
plot(lag2,Ryu11/Ruu0(1,1),'ro'); hold on;
plot(tt,y11,'rx');
title('Ryu(1,1) vs Experimental response of (1,1)')
legend('Ryu(1,1)','Experimental response of (1,1)')
axis([-0.2 2 -0.1 0.06])
grid on;


figure(2)
plot(lag2,Ryu21/Ruu0(1,1),'go'); hold on;
plot(tt,y21,'gx'); 
title('Ryu(2,1) vs Experimental response of (2,1)')
legend('Ryu(2,1)','Experimental response of (2,1)')
axis([-0.2 2 -0.1 0.06])
grid on;

figure(3)
plot(lag2,Ryu12/Ruu0(2,2),'bo'); hold on;
plot(tt,y12,'bx'); 
title('Ryu(1,2) vs Experimental response of (1,2)')
legend('Ryu(1,2)','Experimental response of (1,2)')
axis([-0.2 2 -0.1 0.06])
grid on;

figure(4)
plot(lag2,Ryu22/Ruu0(2,2),'ko'); hold on;
plot(tt,y22,'kx'); 
title('Ryu(2,2) vs Experimental response of (2,2)')
legend('Ryu(2,2)','Experimental response of (2,2)')
axis([-0.2 2 -0.1 0.06])
grid on;



