function [C,c,y] = circles(n,m,v)
%For a resolved tangle v on m strands, this creates a matrix to name the
%circles, n is used to check if any of the first n strands backtrack
x=length(v);
C=zeros(x,m);
c=1;
b=1;
for i=1:x
for j=1:m
    if C(i,j)==0
        C(i,j)=c;
        k=1; %back to start=0
        d=0; %direction
        I=i;
        J=j;
        while k~=0
            if d==0
                if v(I)==J
                    J=J+1;
                    C(I,J)=c;
                    d=1;
                elseif v(I)==J-1
                    J=J-1;
                    C(I,J)=c;
                    d=1;
                else
                    I=mod(I,x)+1;
                    C(I,J)=c;
                end
            elseif d==1
                if (I>1 && v(I-1)==J) || (I==1 && v(x)==J)
                    J=J+1;
                    C(I,J)=c;
                    d=0;
                elseif (I>1 && v(I-1)==J-1) || (I==1 && v(x)==J-1)
                    J=J-1;
                    C(I,J)=c;
                    d=0;
                elseif I==1
                    I=x;
                    C(I,J)=c;
                else 
                    I=I-1;
                    C(I,J)=c;
                end
            end
            if i==1 && j<=n && I==1 && J<=n && d==1
                b=0;
            end
            if I==i && J==j
                k=0;
                c=c+1;
            end
        end
    end
end
end
if b==0
    c=0;
else
c=c-1;
end
y=max(C(1,1:n));
end

