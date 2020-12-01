

function gamma = HinfnormofPd(A,B,C,D,g_upper,g_lower,tol)
    Ac=-inv(eye(length(A))+A)*(eye(length(A))-A);
    Bc=sqrt(2)*inv(eye(length(A))+A)*B;
    Cc=sqrt(2)*C*inv(eye(length(A))+A);
    Dc=D-C*inv(eye(length(A))+A)*B;
    g_lower=round(g_lower);
    g_upper=round(g_upper);
    tol=round(tol);

    array=[0.4:0.0001:0.5];
    N=length(array);
    d=length(Dc'*Dc);
    Aclp=[];
    gamma=0;
    for i=1:N;
        Dg=[(array(1,i)^2)*eye(d)-Dc'*Dc];
        Aclp=[Ac+Bc*inv(Dg)*Dc'*Cc Bc*inv(Dg)*Bc';-Cc'*Cc-Cc'*Dc*inv(Dg)*Dc'*Cc -Ac'-Cc'*Dc*inv(Dg)*Bc'];
        evals=eig(Aclp);
        for k=1:length(Aclp);
            v=evals(k,1)
            if j*imag(v)==v;
                gamma=i
                break   
            end
        end
        if gamma==i
            break
        end
    end
    gamma = array(1,i);
end
