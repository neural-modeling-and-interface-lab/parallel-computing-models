function y=lagu_re(k,n,p)
%k is sorder of the basis
%p is alpha
%n is the lenght of the sequence
nn=sqrt(1-p^2);n1=[0,nn];p1=[1,-p];
[a1,b1,c1,d1]=tf2ss(n1,p1);
n0=[-p,1];p0=[1,-p];
[a0,b0,c0,d0]=tf2ss(n0,p0);
a=a1;b=b1;c=c1;d=d1;

for i=2:k; 
    [a,b,c,d]=series(a0,b0,c0,d0,a,b,c,d);
end

y1=dimpulse(a,b,c,d,1,n+1);
y=y1(2:n+1);




