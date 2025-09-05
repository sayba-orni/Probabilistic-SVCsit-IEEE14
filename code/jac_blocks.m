function [J1,J2,J3,J4] = jac_blocks(Vmag,Vang,Y)
N=numel(Vmag); n=N-1;
J1=zeros(n,n); J2=zeros(n,n); J3=zeros(n,n); J4=zeros(n,n);
[P,Q] = pq_injections(Vmag,Vang,Y);
for i=2:N
    for k=2:N
        if i==k
            Yii = Y(i,i); gii=abs(Yii); thi=angle(Yii);
            J1(i-1,k-1) = -Q(i) - Vmag(i)^2*gii*sin(thi);
            J2(i-1,k-1) =  P(i) + Vmag(i)^2*gii*cos(thi);
            J3(i-1,k-1) =  P(i) - Vmag(i)^2*gii*cos(thi);
            J4(i-1,k-1) =  Q(i) - Vmag(i)^2*gii*sin(thi);
        else
            th = angle(Y(i,k)) + Vang(k) - Vang(i);
            g  = abs(Y(i,k));
            J1(i-1,k-1) = -Vmag(i)*Vmag(k)*g*sin(th);
            J2(i-1,k-1) =  Vmag(i)*Vmag(k)*g*cos(th);
            J3(i-1,k-1) = -Vmag(i)*Vmag(k)*g*cos(th);
            J4(i-1,k-1) = -Vmag(i)*Vmag(k)*g*sin(th);
        end
    end
end
end
