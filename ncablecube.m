function [h] = ncablecube(n,m,v)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
w=ncable(n,m,v);
M=n*m;
x=length(w);
K=zeros(x,M,2^x); %cube
C=zeros(2^x,1); %number of circles
T=zeros(2^x,1); %colours used in ignored tangle strands
R=zeros(2^x,x); %resolutions
U=zeros(2^x,1); %homological grading of diagrams
maxq=-1;
minq=x+2;
for i=1:2^x
    for j=1:x
        R(i,j)=bitget(i-1,j);
    end
    U(i)=sum(R(i,:));
    [K(:,:,i),C(i),T(i)]=circles(n,M,resolution(M,w,R(i,:)));
    maxq=max(maxq,U(i)+C(i)-T(i));
    minq=min(minq,U(i)-C(i)+T(i));
end
P=ones(maxq-minq+1,x+1);
Z=zeros(2^x,maxq-minq+1); %start dimension for each resolution
r=zeros(maxq-minq+1,x+1); %ranks of total complex
for i=1:2^x
    if C(i)~=0
        Z(i,:)=P(:,U(i)+1);
        for j=0:C(i)-T(i)
            r(U(i)-C(i)+T(i)+2*j-minq+1,U(i)+1)=r(U(i)-C(i)+T(i)+2*j-minq+1,U(i)+1)+nchoosek(C(i)-T(i),j);
            P(U(i)-C(i)+T(i)+2*j-minq+1,U(i)+1)=P(U(i)-C(i)+T(i)+2*j-minq+1,U(i)+1)+nchoosek(C(i)-T(i),j);
        end
    else
        Z(i,:)=P(:,U(i)+1);
    end
end
D=cell(maxq-minq+1,x); %differential
for i=1:x
    for j=1:maxq-minq+1
        D{j,i}=sparse(r(j,i+1),r(j,i));
    end
end
L=cell(2^x,1);
quant=cell(2^x,1);
for i=1:2^x
%    if i/1000-floor(i/1000)==0
%        i
%    end
   if C(i)==0
       continue
   end
   quant{i}=zeros(2^(C(i)-T(i)),1);
   L{i}=zeros(2^(C(i)-T(i)),C(i)-T(i)); %generators of domain
   counter=zeros(maxq-minq+1,1);
   for l=1:2^(C(i)-T(i))
       for j=1:C(i)-T(i)
            L{i}(l,j)=bitget(l-1,j);
       end
       quant{i}(l)=counter(sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq);
       counter(sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq)=counter(sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq)+1;
   end
end
for i=1:2^x
   P=find(R(i,:)==0);
   p=length(P);
   for j=1:p
       s=i+2^(P(j)-1);
       if C(s)==0
           continue
       end
       [perm]=quickcobordism(M,w,K(:,:,i),C(i),K(:,:,s),P(j));
       q=setdiff(T(i)+1:C(s),perm);
       if length(find(perm==0))==2 && isempty(setdiff(1:T(i),perm))==1
           for l=1:2^(C(i)-T(i))
               qgr=sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq;
               g=1;
               for u=1:C(i)-T(i)
                    if perm(u+T(i))==0 || L{i}(l,u)==0
                        continue
                    end
                    g=g+2^(perm(u+T(i))-T(s)-1);
               end
               B=L{i}(l,find(perm==0)-T(i));
               if isequal(B,[0,0])==1
                   D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g),Z(i,qgr)+quant{i}(l))=1;
               elseif isequal(B,[1,0]) || isequal(B,[0,1])
                   D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g+2^(q-T(s)-1)),Z(i,qgr)+quant{i}(l))=1;
               end
           end
       elseif length(find(perm==0))==2
           for l=1:2^(C(i)-T(i))
               B=L{i}(l,find(perm==0,1,'last')-T(i));
               if B==0
                   qgr=sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq;
                   g=1;
                   for u=1:C(i)-T(i)
                       if perm(u+T(i))==0 || L{i}(l,u)==0
                           continue
                       end
                       g=g+2^(perm(u+T(i))-T(s)-1);
                   end
                   D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g),Z(i,qgr)+quant{i}(l))=1;
               end
           end
       elseif isempty(setdiff(1:T(i),perm))==1
           for l=1:2^(C(i)-T(i))
               B=L{i}(l,find(perm==0)-T(i));
               qgr=sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq;
               g=1;
               J=setdiff(1:C(s),perm);
               for u=1:C(i)-T(i)
                   if perm(u+T(i))==0 || L{i}(l,u)==0
                       continue
                   end
                   g=g+2^(perm(u+T(i))-T(s)-1);
               end
               if B==0
                   D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g+2^(J(1)-T(s)-1)),Z(i,qgr)+quant{i}(l))=1;
                   D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g+2^(J(2)-T(s)-1)),Z(i,qgr)+quant{i}(l))=1;
               else
                   D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g+2^(J(1)-T(s)-1)+2^(J(2)-T(s)-1)),Z(i,qgr)+quant{i}(l))=1;
               end
           end
       else
           for l=1:2^(C(i)-T(i))
               qgr=sum(L{i}(l,:)==0)-sum(L{i}(l,:))+U(i)+1-minq;
               g=1;
               J=setdiff(T(s)+1:C(s),perm);
               for u=1:C(i)-T(i)
                   if perm(u+T(i))==0 || L{i}(l,u)==0
                       continue
                   end
                   g=g+2^(perm(u+T(i))-T(s)-1);
               end
               D{qgr,U(i)+1}(Z(s,qgr)+quant{s}(g+2^(J-T(s)-1)),Z(i,qgr)+quant{i}(l))=1;
           end
       end
   end 
end
%Now calculate homology
s=zeros(maxq-minq+1,x+1);
h=zeros(maxq-minq+1,x+1);
for j=1:maxq-minq+1
    if mod(j+minq-n,2)==0
        for i=1:x+1
            h(i,j)=0;
        end
    else
        for i=1:x
            E=refmod(D{j,i},n);
            [F,~,~]=find(E);
            if isempty(F)==1
                s(j,i)=0;
            else
                s(j,i)=max(F);
            end
            if i==1
                h(j,i)=r(j,i)-s(j,i);
            else
                h(j,i)=r(j,i)-s(j,i)-s(j,i-1);
            end
        end
        h(j,x+1)=r(j,x+1)-s(j,x);
    end
end
end