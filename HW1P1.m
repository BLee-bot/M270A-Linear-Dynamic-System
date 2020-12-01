A = [6 7 15;-1 -1 -3;-1 -1 -2];
[V,E] = eig(A);
v1 = V(:,1); %% can use one eigenvector of A to start the orthogonal basis
w2 = [0; 1; 1]; %% play with w2 and w3 to generate different orthog basis
w3 = [0; 0; 1];
v2tmp = w2 - (w2'*v1)*v1; %%% Gram-Schmidt
v2 = v2tmp/norm(v2tmp);
v3tmp = w3 - (w3'*v1)*v1 - (w3'*v2)*v2; %%% Gram-Schmidt
v3 = v3tmp/norm(v3tmp);
T1 = [v1 v2 v3]; %%%% unitary
Anew = T1'*(A*T1);
A2 = Anew(2:3,2:3);
[V,E] = eig(A2); %% can use one eigenvector of A2
v1 = V(:,1);
w4 = [1; 1]; %% change this for different orthogonal basis
v2tmp = w4 - (w4'*v1)*v1; %%% Gram-Schmidt
v2 = v2tmp/norm(v2tmp);
T2 = [v1 v2]; %%% unitary
%%% final T matrix
T = T1*[1, zeros(1,2); zeros(2,1), T2];
disp(T'*(A*T)) %% show transformed matrix in upper triangular form