

function gamma = HinfnormofP(Ac, Bc, Cc, Dc,g_upper,g_lower,tol)
    Gamma=[g_lower:tol:g_upper]
    N=length(Gamma)
    d=length(Dc'*Dc)
    Aclp=[]
    gamma=0
    for i=1:N
        Dg=[(Gamma(1,i)^2)*eye(d)-Dc'*Dc]
        Aclp=[Ac+Bc*inv(Dg)*Dc'*Cc -Bc*inv(Dg)*Bc';Cc'*Cc+Cc'*Dc*inv(Dg)*Dc'*Dc -Ac'-Cc'*Dc*inv(Dg)*Bc']
        evals=eig(Aclp)
%         for k=1:dim(Aclp);
%             v=evals(k,1);
%             if j*imag(v)==v;
%                 gamma=i;
%                 break   
%             end
%         end
%         if gamma==i;
%             break
%         end
    end
    gamma = Gamma(1,i);
end
