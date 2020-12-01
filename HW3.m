%%HW3
%Problem1
    A1=[1 1; 2 1];
    A2=[1 1; 1 2];
    A3=[-1 1; -1 -2];
    A4=[1 1; 1 1]; %Hermitian
    [U1,S1,V1]=svd(A1)
    [U2,S2,V2]=svd(A2)
    [U3,S3,V3]=svd(A3)
    [U4,S4,V4]=svd(A4)
    t=linspace(0,2*pi,360); %점을 일정하게 나누기 위해 linspace를 사용하는 것이 좋습니다.
    figure(1)
    axis square;
    hold on;
    plot(cos(t),sin(t),'b',cos(t)+sin(t),2*cos(t)+sin(t),'g')
    plot([0,U1(1,1)],[0,U1(1,2)],'r')
    plot([0,U1(2,1)],[0,U1(2,2)],'y')
    plot([0,V1(1,1)],[0,V1(2,1)],'-r')
    plot([0,V1(1,2)],[0,V1(2,2)],'-y')
    title('Unit Circle C and Ellipsoid E')
    xlabel('x1')
    ylabel('x2')
    legend('C','E')
    grid ON
    xlim([-3,3])
    ylim([-3,3])
    hold off
    



%eigshow(_)

%Problem2
    %{ 
    M=[3 2 1 0;1 1 1 1];
    [U,S,V]=svd(M,0);
    rank(M);
    s=diag(S);
    U(:,logical(s));
    V2=null(M,'r')
    u=[0.1;0.2;0.3;0.4];
    V3=0*V2(:,1)+0*V2(:,2);
    unew=u+V3;
    M*unew;
    norm(unew)^2; 
    %}
%Problem3
    %{
    C=[1 1];
    A=[1 1;0 1];
    M=[C;C*A;C*A^2;C*A^3]
    y=[0 ;2 ;4 ;6]
    x0=inv(M'*M)*M'*y
    %}
%Problem4
    %{
    y=[0;2;4;6;8;10;12;14];
    H=[y(1:4) y(2:5) y(3:6) y(4:7)];
     rank(H);
    C=[1 1];
    A=[1 1;0 1];
   	M=[C;C*A;C*A^2;C*A^3];
    x01=[-2 ;2];
    N=[x0 A*x0 (A^2)*x0 (A^3)*x0];
        M*N;
    H2=[y(2:5) y(3:6) y(4:7) y(5:8)];
    % variable T
        t11=sym('t11');
        t12=sym('t12');
        t21=sym('t21');
        t22=sym('t22');
    T=[t11 t12;t21 t22];
    t=[t11; t12; t21; t22];
    A2=simplify(inv(T)*A*T)
    C2=simplify(C*T)
    X02=simplify(inv(T)*x01)
    %}
        function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
hold off
    end
    
    
    