function A = refmod(A, p)
[m,n] = size(A);
i=1;
j=1;
A=mod(A,p);
while (i<=m && j<=n) 
    [~,maxi]=max(A(i:end,j));
    maxi=maxi+i-1;
    if(A(maxi,j)~=0)
        rows=1:m;
        rows(i)=maxi;
        rows(maxi)=i;
        A=A(rows,:);
        A(i,:)=mod(A(i,:)*invMod(A(i,j),p),p);
        A(i+1:end,:)=A(i+1:end,:)-A(i+1:end,j)*A(i,:);
        A(i+1:end,:)=mod(A(i+1:end,:),p);
        i=i+1;
    end
    j=j+1;
  end
end
  
  
  
function y =invMod(x,p)
[~,c,~]=gcd(x,p); 
y=mod(c,p);
end