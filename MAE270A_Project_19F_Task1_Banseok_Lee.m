
%% Setting
clear
close all
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

% Original Graph
%{
figure(1);
subplot(311)
plot(t,u1,'b*','LineWidth',2)
ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])
subplot(312)
plot(t,y11,'r*','LineWidth',2)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])
subplot(313)
plot(t,y21,'r*','LineWidth',2)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])
figure(2);
subplot(311)
plot(t,u2,'b*','LineWidth',2)
ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 1.1])
subplot(312)
plot(t,y12,'r*','LineWidth',2)
ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
axis([-0.2 2 -0.1 0.1])
subplot(313)
plot(t,y22,'r*','LineWidth',2)
ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
grid on
xlabel('second','FontSize',14)
set(gca,'FontSize',14)
axis([-0.2 2 -0.1 0.1])
%} 

% Task 1

%% Task 1-1
%%% H100 Matrix construct
H100=[]; Hk1=[]; n=42; k=0;
while n<142
    while k<100
    Hk1=[Hk1;y11(n+k) y12(n+k);y21(n+k) y22(n+k)];
    k=k+1;
    end
H100=[H100, Hk1];
Hk1=[];
n=n+1;
k=0;
end
%%% H101 Matrix construct
H101=[]; Hk2=[];n=43;k=0;
while n<143
    while k<100
    Hk2=[Hk2;y11(n+k) y12(n+k);y21(n+k) y22(n+k)];
    k=k+1;
    end
H101=[H101, Hk2];
Hk2=[];
n=n+1;
k=0;
end
%%% Singular Value Decomposition of H100
[U,S,V] = svd(H100);              
U106=U(:,[1,2,3,4,5,6]); V106=V(:,[1,2,3,4,5,6]); S106=S([1,2,3,4,5,6],[1,2,3,4,5,6]);
U107=U(:,[1,2,3,4,5,6,7]); V107=V(:,[1,2,3,4,5,6,7]); S107=S([1,2,3,4,5,6,7],[1,2,3,4,5,6,7]);
U110=U(:,[1,2,3,4,5,6,7,8,9,10]); V110=V(:,[1,2,3,4,5,6,7,8,9,10]); S110=S([1,2,3,4,5,6,7,8,9,10],[1,2,3,4,5,6,7,8,9,10]);
U120=U(:,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]); V120=V(:,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]); S120=S([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]);

%%% Figure of Singular values
% figure('Name', 'Singular Values')
% X=1:20;
% X=X';
% title('Singular Values')
% figure1=semilogy(X,[diag(S106);zeros(14,1)],'ro',X,[diag(S107);zeros(13,1)],'go',X,[diag(S110);zeros(10,1)],'bo',X,diag(S120),'ko')
% figure1(1).MarkerSize=10;
% figure1(2).MarkerSize=8;
% figure1(3).MarkerSize=6;
% figure1(4).MarkerSize=4;
% legend('ns=06','ns=07','ns=10','ns=20')
% axis([0 20 0 1])
% xlabel('Singular value Index')
% ylabel('Habkel singular value')
% grid on

A06=U106'*H101*V106*inv(S106);                  %%% Construct A06
A07=U107'*H101*V107*inv(S107);                  %%% Construct A07
A10=U110'*H101*V110*inv(S110);                  %%% Construct A10
A20=U120'*H101*V120*inv(S120);                  %%% Construct A20
max(abs(eig(A06)))                              %%% Max Singular of A06
max(abs(eig(A07)))                              %%% Max Singular of A07
max(abs(eig(A10)))                              %%% Max Singular of A10
max(abs(eig(A20)))                              %%% Max Singular of A20

%% Taks 1-2
%%% Impulse Response of A06                                            
Y061=[0;0]; Y062=[0;0];
N06=U106*(S106)^0.5; C06=N06([1,2],:);
M06=((S106)^0.5)*V106'; B06=M06(:,[1,2]);
Xk=zeros(1,6);
Xk=Xk';
k=1;
while k<401
Xk=A06*Xk+B06*[u1(1,k);0];    
Yk=C06*Xk;
Y061=[Y061, Yk];
k=k+1;
end
k=1;
while k<401
Xk=A06*Xk+B06*[0;u2(1,k)];    
Yk=C06*Xk;
Y062=[Y062, Yk];
k=k+1;
end
Y0611=Y061(1,:);
Y0612=Y061(2,:);
Y0621=Y062(1,:);
Y0622=Y062(2,:);

%%% Graph of A06

% figure('Name','Impulse response of A06 for input u1');
% 
% subplot(511);
% plot(t,u1,'b*','LineWidth',2);
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u1','Interpreter','Latex');
% 
% subplot(512);
% plot(t,y11,'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y0611,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A06','Interpreter','Latex')
% 
% subplot(514);
% plot(t,y21,'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y0612,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A06','Interpreter','Latex')
% xlabel('second','FontSize',14)
% sgtitle('Impulse response of A06 for input u1');
% 
% 
% figure('Name','Impulse response of A06 for input u2');
% 
% subplot(511)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u2','Interpreter','Latex');
% 
% subplot(512)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y0621,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A06','Interpreter','Latex')
% 
% subplot(514)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y0622,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% xlabel('second','FontSize',14)
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A06','Interpreter','Latex')
% sgtitle('Impulse response of A06 for input u2');
% 

%%% Impulse Response of A07
Y071=[0;0]; Y072=[0;0];
N07=U107*(S107)^0.5; C07=N07([1,2],:);
M07=((S107)^0.5)*V107'; B07=M07(:,[1,2]);
Xk=zeros(1,7); Xk=Xk'; k=1;
while k<401
Xk=A07*Xk+B07*[u1(1,k);0];    
Yk=C07*Xk;
Y071=[Y071, Yk];
k=k+1;
end
k=1;
while k<401
Xk=A07*Xk+B07*[0;u2(1,k)];    
Yk=C07*Xk;
Y072=[Y072, Yk];
k=k+1;
end
Y0711=Y071(1,:);
Y0712=Y071(2,:);
Y0721=Y072(1,:);
Y0722=Y072(2,:);

%%% Graph of A07


% figure('Name','Impulse response of A07 for input u1');
% 
% subplot(511);
% plot(t,u1,'b*','LineWidth',2);
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u1','Interpreter','Latex');
% 
% subplot(512);
% plot(t,y11,'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y0711,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A07','Interpreter','Latex')
% 
% subplot(514);
% plot(t,y21,'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y0712,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A07','Interpreter','Latex')
% xlabel('second','FontSize',14)
% sgtitle('Impulse response of A07 for input u1');
% 
% 
% figure('Name','Impulse response of A07 for input u2');
% 
% subplot(511)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u2','Interpreter','Latex');
% 
% subplot(512)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y0721,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A07','Interpreter','Latex')
% 
% subplot(514)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y0722,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% xlabel('second','FontSize',14)
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A07','Interpreter','Latex')
% sgtitle('Impulse response of A07 for input u2');



%%% Impulse Response of A10
Y101=[0;0]; Y102=[0;0];
N10=U110*(S110)^0.5; C10=N10([1,2],:);
M10=((S110)^0.5)*V110'; B10=M10(:,[1,2]);
Xk=zeros(1,10); Xk=Xk'; k=1;
while k<401
Xk=A10*Xk+B10*[u1(1,k);0];    
Yk=C10*Xk;
Y101=[Y101, Yk];
k=k+1;
end
k=1;
while k<401
Xk=A10*Xk+B10*[0;u2(1,k)];    
Yk=C10*Xk;
Y102=[Y102, Yk];
k=k+1;
end
Y1011=Y101(1,:);
Y1012=Y101(2,:);
Y1021=Y102(1,:);
Y1022=Y102(2,:);

%%% Graph of A10

% figure('Name','Impulse response of A10 for input u1');
% 
% subplot(511);
% plot(t,u1,'b*','LineWidth',2);
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u1','Interpreter','Latex');
% 
% subplot(512);
% plot(t,y11,'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y1011,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A10','Interpreter','Latex')
% 
% subplot(514);
% plot(t,y21,'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y1012,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A10','Interpreter','Latex')
% xlabel('second','FontSize',14)
% sgtitle('Impulse response of A10 for input u1');
% 
% 
% figure('Name','Impulse response of A10 for input u2');
% 
% subplot(511)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u2','Interpreter','Latex');
% 
% subplot(512)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y1021,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A10','Interpreter','Latex')
% 
% subplot(514)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y1022,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% xlabel('second','FontSize',14)
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A10','Interpreter','Latex')
% sgtitle('Impulse response of A10 for input u2');

%%% Impulse Response of A20
Y201=[0;0]; Y202=[0;0];
N20=U120*(S120)^0.5; C20=N20([1,2],:);
M20=((S120)^0.5)*V120'; B20=M20(:,[1,2]);
Xk=zeros(1,20); Xk=Xk'; k=1;
while k<401
Xk=A20*Xk+B20*[u1(1,k);0];    
Yk=C20*Xk;
Y201=[Y201, Yk];
k=k+1;
end
k=1;
while k<401
Xk=A20*Xk+B20*[0;u2(1,k)];    
Yk=C20*Xk;
Y202=[Y202, Yk];
k=k+1;
end
Y2011=Y201(1,:);
Y2012=Y201(2,:);
Y2021=Y202(1,:);
Y2022=Y202(2,:);

%%% Graph of A20

% figure('Name','Impulse response of A20 for input u1');
% 
% subplot(511);
% plot(t,u1,'b*','LineWidth',2);
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u1','Interpreter','Latex');
% 
% subplot(512);
% plot(t,y11,'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y2011,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A20','Interpreter','Latex')
% 
% subplot(514);
% plot(t,y21,'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y2012,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A20','Interpreter','Latex')
% xlabel('second','FontSize',14)
% sgtitle('Impulse response of A20 for input u1');
% 
% 
% figure('Name','Impulse response of A20 for input u2');
% 
% subplot(511)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u2','Interpreter','Latex');
% 
% subplot(512)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(513)
% plot(t,Y2021,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A20','Interpreter','Latex')
% 
% subplot(514)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(515)
% plot(t,Y2022,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% xlabel('second','FontSize',14)
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A20','Interpreter','Latex')
% sgtitle('Impulse response of A20 for input u2');


% figure('Name','Errors');
% subplot(4,4,1);
% temp=y11-Y0611;
% plot(t,abs(temp),'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A06','Interpreter','Latex')
% 
% subplot(4,4,5);
% plot(t,abs(y21-Y0612),'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A06','Interpreter','Latex')
% 
% subplot(4,4,9)
% plot(t,abs(y12-Y0621),'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A06','Interpreter','Latex')
% 
% subplot(4,4,13)
% plot(t,abs(y22-Y0622),'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A06','Interpreter','Latex')
% 
% subplot(4,4,2);
% plot(t,abs(y11-Y0711),'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A07','Interpreter','Latex')
% 
% subplot(4,4,6);
% plot(t,abs(y21-Y0712),'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A07','Interpreter','Latex')
% 
% subplot(4,4,10)
% plot(t,abs(y12-Y0721),'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A07','Interpreter','Latex')
% 
% subplot(4,4,14)
% plot(t,abs(y22-Y0722),'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A07','Interpreter','Latex')
% 
% subplot(4,4,3);
% plot(t,abs(y11-Y1011),'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A10','Interpreter','Latex')
% 
% subplot(4,4,7);
% plot(t,abs(y21-Y1012),'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A10','Interpreter','Latex')
% 
% subplot(4,4,11)
% plot(t,abs(y12-Y1021),'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A10','Interpreter','Latex')
% 
% subplot(4,4,15)
% plot(t,abs(y22-Y1022),'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A10','Interpreter','Latex')
% 
% subplot(4,4,4);
% plot(t,abs(y11-Y2011),'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A20','Interpreter','Latex')
% 
% subplot(4,4,8);
% plot(t,abs(y21-Y2012),'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_1$ (volts) for A20','Interpreter','Latex')
% 
% subplot(4,4,12)
% plot(t,abs(y12-Y2021),'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A20','Interpreter','Latex')
% 
% subplot(4,4,16)
% plot(t,abs(y22-Y2022),'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 1 -0.1 0.1])
% title('Error of $y_2$ (volts) for A20','Interpreter','Latex')
% sgtitle('Impulse response Error of A06,A07,A10,A20');


%% Taks 1-3
wnyq=20;
ws=40;
ts=1/40;
D=zeros(2,2);
W={0,wnyq};
sys06=ss(A06,B06,C06,D,ts);
sys07=ss(A07,B07,C07,D,ts);
sys10=ss(A10,B10,C10,D,ts);
sys20=ss(A20,B20,C20,D,ts);
[mag06,phase06,wout06]=bode(sys06,W);     
[mag07,phase07,wout07]=bode(sys07,W);    
[mag10,phase10,wout10]=bode(sys10,W);    
[mag20,phase20,wout20]=bode(sys20,W);    


% figure('Name','Response of (1,1) element')
% opts1 = bodeoptions;
% opts1.Title.FontSize = 14;
% opts1.Title.Interpreter = 'Latex';
% opts1.Title.String = 'Response of (1,1) element';
% opts1.FreqUnits = 'Hz';
% bodeplot(sys06(1,1),'r',sys07(1,1),'y',sys10(1,1),'g',sys20(1,1),'b',opts1)
% legend('ns=06','ns=07','ns=10','ns=20')
% grid on
% 
% figure('Name','Response of (1,2) element')
% opts2 = bodeoptions;
% opts2.Title.FontSize = 14;
% opts2.Title.Interpreter = 'Latex';
% opts2.Title.String = 'Response of (1,2) element';
% opts2.FreqUnits = 'Hz';
% bodeplot(sys06(1,2),'r',sys07(1,2),'y',sys10(1,2),'g',sys20(1,2),'b',opts2)
% legend('ns=06','ns=07','ns=10','ns=20')
% grid on
% 
% figure('Name','Response of (2,1) element')
% opts3 = bodeoptions;
% opts3.Title.FontSize = 14;
% opts3.Title.Interpreter = 'Latex';
% opts3.Title.String = 'Response of (2,1) element';
% opts3.FreqUnits = 'Hz';
% bodeplot(sys06(2,1),'r',sys07(2,1),'y',sys10(2,1),'g',sys20(2,1),'b',opts3)
% legend('ns=06','ns=07','ns=10','ns=20')
% grid on
% 
% figure('Name','Response of (2,2) element')
% opts4 = bodeoptions;
% opts4.Title.FontSize = 14;
% opts4.Title.Interpreter = 'Latex';
% opts4.Title.String = 'Response of (2,2) element';
% opts4.FreqUnits = 'Hz';
% bodeplot(sys06(2,2),'r',sys07(2,2),'y',sys10(2,2),'g',sys20(2,2),'b',opts4)
% legend('ns=06','ns=07','ns=10','ns=20')
% grid on

%% Taks 1-4


y11f = fft(y11)./fft(u1) ;
y21f = fft(y21)./fft(u1) ;
y12f = fft(y12)./fft(u2) ;
y22f = fft(y22)./fft(u2) ;

N = length(y11f);
om = [0:N-1]/(ts*N); %%%% frequency vector in hertz

sys11 = frd(y11f,om,ts);
sys11.FrequencyUnit='Hz';
sys12 = frd(y12f,om,ts);
sys12.FrequencyUnit='Hz';
sys21 = frd(y21f,om,ts);
sys21.FrequencyUnit='Hz';
sys22 = frd(y22f,om,ts);
sys22.FrequencyUnit='Hz';

figure('Name','Response of (1,1) element')
opts1 = bodeoptions;
opts1.Title.FontSize = 14;
opts1.Title.Interpreter = 'Latex';
opts1.Title.String = 'Response of (1,1) element';
opts1.FreqUnits = 'Hz';
bodeplot(sys06(1,1),'r',sys07(1,1),'y',sys10(1,1),'g',sys20(1,1),'b',sys11,'k',opts1)
legend('ns=06','ns=07','ns=10','ns=20','direct model')
grid on

figure('Name','Response of (1,2) element')
opts2 = bodeoptions;
opts2.Title.FontSize = 14;
opts2.Title.Interpreter = 'Latex';
opts2.Title.String = 'Response of (1,2) element';
opts2.FreqUnits = 'Hz';
bodeplot(sys06(1,2),'r',sys07(1,2),'y',sys10(1,2),'g',sys20(1,2),'b',sys12,'k',opts2)
legend('ns=06','ns=07','ns=10','ns=20','direct model')
grid on

figure('Name','Response of (2,1) element')
opts3 = bodeoptions;
opts3.Title.FontSize = 14;
opts3.Title.Interpreter = 'Latex';
opts3.Title.String = 'Response of (2,1) element';
opts3.FreqUnits = 'Hz';
bodeplot(sys06(2,1),'r',sys07(2,1),'y',sys10(2,1),'g',sys20(2,1),'b',sys21,'k',opts3)
legend('ns=06','ns=07','ns=10','ns=20','direct model')
grid on

figure('Name','Response of (2,2) element')
opts4 = bodeoptions;
opts4.Title.FontSize = 14;
opts4.Title.Interpreter = 'Latex';
opts4.Title.String = 'Response of (2,2) element';
opts4.FreqUnits = 'Hz';
bodeplot(sys06(2,2),'r',sys07(2,2),'y',sys10(2,2),'g',sys20(2,2),'b',sys22,'k',opts4)
legend('ns=06','ns=07','ns=10','ns=20','direct model')
grid on





