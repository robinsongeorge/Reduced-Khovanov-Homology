function [perm] = quickcobordism(m,v,C,c,D,y)
%given a tangle v on m strands, a resolution r, and an integer y st r(y)=-1, this outputs a map from the unaffected colours
%of the first tangle to the unaffected colours of the second
x=length(v);
I=mod(abs(v(y)),m)+1;
J=mod(y,x)+1;
pre=union(C(y,abs(v(y))),C(J,I));
perm=zeros(1,c);
for k=1:c
    if ismember(k,pre)==1
        perm(k)=0;
    else 
        perm(k)=D(find(C==k,1));
    end
end
end