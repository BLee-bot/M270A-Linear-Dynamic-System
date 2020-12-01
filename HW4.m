%%HW4
%Problem1
A=[1 1 1; 1 1 1 ; 1 2 3 ; 1 2 4];
det(A'*A);
a1=A(:, 1);
a2=A(:, 2);
a3=A(:, 3);

u1=a1;
u2=a2-((a1'*a2)/(a1'*a1))*a1;
u3=a3-((a1'*a3)/(a1'*a1))*a1-((u2'*a3)/(u2'*u2))*u2;

e1=u1/norm(u1);
e2=u2/norm(u2);
e3=u3/norm(u3);

Q=[e1 e2 e3];
R=[2 3 4.5; 0 1 2.5; 0 0 0.5/0.7071];
Q*R;
y=[2; 1; 1; 2];
xls=inv(A'*A)*A'*y;
inv(R)*Q'*y;
%problem2
%problem3

