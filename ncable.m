function [w] = ncable(n,m,v)
%n cable of braid with m strands denoted by v
x=length(v);
w=zeros(x*n^2,1);
for i=1:x
    e=v(i)/abs(v(i));
    for j=1:n^2
        w((i-1)*n^2+j)=e*((abs(v(i))-1)*n+(1+n)*(floor((j-1)/n)+1)-j);
    end
end
end

