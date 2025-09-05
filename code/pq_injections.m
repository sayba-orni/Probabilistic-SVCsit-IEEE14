function [P,Q] = pq_injections(Vmag, Vang, Y)
N=numel(Vmag); P=zeros(1,N); Q=zeros(1,N);
for i=1:N
    for k=1:N
        th = angle(Y(i,k)) + Vang(k) - Vang(i);
        g  = abs(Y(i,k));
        P(i) = P(i) + Vmag(i)*Vmag(k)*g*cos(th);
        Q(i) = Q(i) - Vmag(i)*Vmag(k)*g*sin(th);
    end
end
end
