
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
U107=U(:,[1,2,3,4,5,6,7]); V107=V(:,[1,2,3,4,5,6,7]); S107=S([1,2,3,4,5,6,7],[1,2,3,4,5,6,7]);
A07=U107'*H101*V107*inv(S107);                  %%% Construct A07
max(abs(eig(A07)))    ;                          %%% Max Singular of A07

Y071=[]; Y072=[];
N07=U107*(S107)^0.5; C07=N07([1,2],:);
M07=((S107)^0.5)*V107'; B07=M07(:,[1,2]);
Xk=zeros(1,7); Xk=Xk'; k=1;
while k<402
Xk=A07*Xk+B07*[u1(1,k);0];    
Yk=C07*Xk;
Y071=[Y071, Yk];
k=k+1;
end
k=1;
while k<402
Xk=A07*Xk+B07*[0;u2(1,k)];    
Yk=C07*Xk;
Y072=[Y072, Yk];
k=k+1;
end
Y0711=Y071(1,:);
Y0712=Y071(2,:);
Y0721=Y072(1,:);
Y0722=Y072(2,:);

%%% Impulse Response of A07 Good


%%% Graph of A07

% figure('Name','Impulse response of A07');
% 
% subplot(521);
% plot(t,u1,'b*','LineWidth',2);
% ylabel('$u_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u1','Interpreter','Latex');
% 
% subplot(523);
% plot(t,y11,'r*','LineWidth',2);
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(525)
% plot(t,Y0711,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A07','Interpreter','Latex')
% 
% subplot(527);
% plot(t,y21,'r*','LineWidth',2);
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(529)
% plot(t,Y0712,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A07','Interpreter','Latex')
% xlabel('second','FontSize',14)
% 
% subplot(522)
% plot(t,u2,'b*','LineWidth',2)
% ylabel('$u_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 1.1])
% title('Input u2','Interpreter','Latex');
% 
% subplot(524)
% plot(t,y12,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_1$ (volts)','Interpreter','Latex')
% 
% subplot(526)
% plot(t,Y0621,'r*','LineWidth',2)
% ylabel('$y_1$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_1$ (volts) for A07','Interpreter','Latex')
% 
% subplot(528)
% plot(t,y22,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% grid on
% set(gca,'FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% title('Real data of $y_2$ (volts)','Interpreter','Latex');
% 
% subplot(5,2,10)
% plot(t,Y0722,'r*','LineWidth',2)
% ylabel('$y_2$ (volts)','FontSize',14,'Interpreter','Latex');
% set(gca,'FontSize',14)
% grid on
% xlabel('second','FontSize',14)
% axis([-0.2 2 -0.1 0.1])
% title('Response data of $y_2$ (volts) for A07','Interpreter','Latex')
% sgtitle('Impulse response of A07');



% Task 2

%% Task 2-1
D=zeros(2,2);
Z=tzero(A07,B07,C07,D); % Transmission zeros of the MIMO model

Z1=Z(1,1); Z2=Z(2,1); Z3=Z(3,1); Z4=Z(4,1); Z5=Z(5,1); % Label Zeros from Largest to Smallest Magnitude

for i = [1,2,3,4,5]
    if abs(Z(i,1)) > 1
        disp(['Z',i,' is unstable'])
    else 
        disp(['Z',i,' is aymptotically stable'])
    end
end
Z11=abs(Z(1,1)); Z22=abs(Z(2,1)); Z33=abs(Z(3,1)); Z44=abs(Z(4,1)); Z55=abs(Z(5,1));
%% Task 2-2

eigvals=eig(A07);

figure(1)
plot(eigvals,'x','MarkerSize',10); hold on;
plot(Z,'o','MarkerSize',10)
theta = linspace(0, 2*pi);
x = cos(theta);
y = sin(theta);
plot(x,y); %// Unit circle
axis([-3 3 -3 3])
axis square
xlabel('Real','FontSize',14,'Interpreter','Latex');
ylabel('Imaginary','FontSize',14,'Interpreter','Latex');
grid on
title('Transmission zeros and eigenvalues')
hold off;

%% Task 2-3

Z1=eigvals(1,1);
eigc1=(log(abs(Z1))+j*atan2(imag(Z1),real(Z1)))/ts
Z2=eigvals(2,1);
eigc2=(log(abs(Z2))+j*atan2(imag(Z2),real(Z2)))/ts
Z3=eigvals(3,1);
eigc3=(log(abs(Z3))+j*atan2(imag(Z3),real(Z3)))/ts
Z4=eigvals(4,1);
eigc4=(log(abs(Z4))+j*atan2(imag(Z4),real(Z4)))/ts
Z5=eigvals(5,1);
eigc5=(log(abs(Z5))+j*atan2(imag(Z5),real(Z5)))/ts
Z6=eigvals(6,1);
eigc6=(log(abs(Z6))+j*atan2(imag(Z6),real(Z6)))/ts
Z7=eigvals(7,1);
eigc7=(log(abs(Z7))+j*atan2(imag(Z7),real(Z7)))/ts

figure(2)
plot([eigc1,eigc2,eigc3,eigc4,eigc5,eigc6,eigc7],'x','MarkerSize',10); hold on;
axis([-11 3 -80 80]);
grid on; hold off;

% 2 damped oscillator (by 2 pair of complex eigenvalues)
dnfreq1=imag(eigc1)/(2*pi)
dnfreq2=imag(eigc3)/(2*pi)

%% Task 2-4

sys=ss(A07,B07,C07,D);
[wn,zeta,p] = damp(sys)

% figure('Name','Pole-Zero Map of discrete-time system')
% 
% subplot(221)
% theta = linspace(0, 2*pi);
% x = cos(theta);
% y = sin(theta);
% plot(x,y); hold on;%// Unit circle
% pzmap(sys(1,1))
% grid on;
% title('Channel (1,1)')
% 
% subplot(222)
% theta = linspace(0, 2*pi);
% x = cos(theta);
% y = sin(theta);
% plot(x,y); hold on;%// Unit circle
% pzmap(sys(1,2))
% grid on;
% title('Channel (1,2)')
% 
% subplot(223)
% theta = linspace(0, 2*pi);
% x = cos(theta);
% y = sin(theta);
% plot(x,y); hold on;%// Unit circle
% pzmap(sys(2,1))
% grid on;
% title('Channel (2,1)')
% 
% subplot(224)
% theta = linspace(0, 2*pi);
% x = cos(theta);
% y = sin(theta);
% plot(x,y); hold on;%// Unit circle
% pzmap(sys(2,2))
% grid on;
% title('Channel (2,2)')
% 
% 
% sgtitle('Pole-Zero Map of discrete-time system')
%% Task 2-5




























%% Taks 1-3
% wnyq=20;
% ws=40;
% ts=1/40;
% D=zeros(2,2);
% W={0,wnyq};
% sys06=ss(A06,B06,C06,D,ts);
% sys07=ss(A07,B07,C07,D,ts);
% sys10=ss(A10,B10,C10,D,ts);
% sys20=ss(A20,B20,C20,D,ts);
% [mag06,phase06,wout06]=bode(sys06,W);     
% [mag07,phase07,wout07]=bode(sys07,W);    
% [mag10,phase10,wout10]=bode(sys10,W);    
% [mag20,phase20,wout20]=bode(sys20,W);    


% figure('Name','Response of (1,1) element')
% opts1 = bodeoptions;
% opts1.Title.FontSize = 14;
% opts1.Title.Interpreter = 'Latex';
% opts1.Title.String = 'Response of (1,1) element';
% opts1.FreqUnits = 'Hz';
% bodeplot(sys06(1,1),'r',sys07(1,1),'y',sys10(1,1),'g',sys20(1,1),'k',opts1)
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

% 
% y11f = fft(y11)./fft(u1) ;
% y21f = fft(y21)./fft(u1) ;
% y12f = fft(y12)./fft(u2) ;
% y22f = fft(y22)./fft(u2) ;
% 
% N = length(y11f);
% om = [0:N-1]/(ts*N); %%%% frequency vector in hertz
% 
% sys11 = frd(y11f,om,ts);
% sys11.FrequencyUnit='Hz';
% sys12 = frd(y12f,om,ts);
% sys12.FrequencyUnit='Hz';
% sys21 = frd(y21f,om,ts);
% sys21.FrequencyUnit='Hz';
% sys22 = frd(y22f,om,ts);
% sys22.FrequencyUnit='Hz';

% figure('Name','Response of (1,1) element')
% opts1 = bodeoptions;
% opts1.Title.FontSize = 14;
% opts1.Title.Interpreter = 'Latex';
% opts1.Title.String = 'Response of (1,1) element';
% opts1.FreqUnits = 'Hz';
% bodeplot(sys06(1,1),'r',sys07(1,1),'y',sys10(1,1),'g',sys20(1,1),'b',sys11,'k',opts1)
% legend('ns=06','ns=07','ns=10','ns=20','direct model')
% grid on
% 
% figure('Name','Response of (1,2) element')
% opts2 = bodeoptions;
% opts2.Title.FontSize = 14;
% opts2.Title.Interpreter = 'Latex';
% opts2.Title.String = 'Response of (1,2) element';
% opts2.FreqUnits = 'Hz';
% bodeplot(sys06(1,2),'r',sys07(1,2),'y',sys10(1,2),'g',sys20(1,2),'b',sys12,'k',opts2)
% legend('ns=06','ns=07','ns=10','ns=20','direct model')
% grid on
% 
% figure('Name','Response of (2,1) element')
% opts3 = bodeoptions;
% opts3.Title.FontSize = 14;
% opts3.Title.Interpreter = 'Latex';
% opts3.Title.String = 'Response of (2,1) element';
% opts3.FreqUnits = 'Hz';
% bodeplot(sys06(2,1),'r',sys07(2,1),'y',sys10(2,1),'g',sys20(2,1),'b',sys21,'k',opts3)
% legend('ns=06','ns=07','ns=10','ns=20','direct model')
% grid on
% 
% figure('Name','Response of (2,2) element')
% opts4 = bodeoptions;
% opts4.Title.FontSize = 14;
% opts4.Title.Interpreter = 'Latex';
% opts4.Title.String = 'Response of (2,2) element';
% opts4.FreqUnits = 'Hz';
% bodeplot(sys06(2,2),'r',sys07(2,2),'y',sys10(2,2),'g',sys20(2,2),'b',sys22,'k',opts4)
% legend('ns=06','ns=07','ns=10','ns=20','direct model')
% grid on


% A07=round(A07,3);
% C07=round(C07,3);
% B07=round(B07,3);
% D07=zeros(2,2);
% Ob = obsv(A07,C07);
% unob=length(A07)-rank(Ob)
% 
% 
% 
% [num1 den1] = ss2tf(A07,B07,C07,D07,1); % iu = 1
% [num2 den2] = ss2tf(A07,B07,C07,D07,2); % iu = 2
% sys = ss(A07,B07,C07,D07);
% s=tf(sys);
% zz=tzero(A07,B07,C07,D07)
% z1=zero(s(1,1))
% z2=zero(s(2,1))
% z3=zero(s(1,2))
% z4=zero(s(2,2))
% s1=s(1,1)
% s2=s(2,1)
% s3=s(1,2)
% s4=s(2,2)
%                          
%                                    
%                                    
%                                    
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   




