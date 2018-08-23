function [w] = resolution(m,v,r)
%given tangle v on m strands, and a resolution vector r, ths outputs the
%resolved tangle
x=length(v);
w=zeros(size(v));
for i=1:x
    e=v(i)/abs(v(i));
    d=2*r(i)-1;
    if d*e<0
        w(i)=-1;
    else
        w(i)=abs(v(i));
    end
end

