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
% figure(1)
% plot(t,u1); hold on;
% plot(t,u2);
%% Task 5-1
u=[u1;u2];
y=[y1;y2]/2;
P=11999;
Yrms=0;
for k=-P:1:P
    Yrms=Yrms+(1/(2*P))*trace(y(1:2,12000+k)*y(1:2,12000+k)');
end
Yrms=sqrt(Yrms);

%% Task 5-2

Goinf=dlyap(A07',C07'*C07);
Gcinf=dlyap(A07,B07*B07');

PH2norm1=trace(B07'*Goinf*B07); % Use (9)
PH2norm1=sqrt(PH2norm1);

PH2norm2=trace(C07*Gcinf*C07'); % Use (10)
PH2norm2=sqrt(PH2norm2);

%% Task 5-3

hk=[];
sum=0;

for i =1:401
    hk=[hk [y11(1,i) y12(1,i);y21(1,i) y22(1,i)]];
end

for i=1:1:401
    hkt=hk(:,2*i-1:2*i);
    sum=sum+trace(hkt*hkt');
end
sum=sqrt(sum);


disp('RMS value scaled by y is')
Yrms
disp('P H2 norm is')
PH2norm1
PH2norm2
disp('Approxiamte P H2 norm by impulse response data is')
sum





