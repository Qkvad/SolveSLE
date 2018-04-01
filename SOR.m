function [x,res]=SOR(A,b,x0,tol,kmax,omega)
  n=max(size(A));
  normab=norm(b);
  
  if normab==0
    x0=zeros(n,1);
    res(1)=0;
    break
  end
  
  res(1)=norm(A*x0-b,2)/normab;
  if(res(1)<tol) break, end
  
  x=x0;
  for k=2:kmax
    for i=1:n
      x(i)=(1-omega)*x(i)+omega*(b(i)-(A(i,1:n)*x(1:n)-A(i,i)*x(i)))/A(i,i);
    end
    res(k)=norm(A*x-b,2)/normab;
    if(res(k)<=tol) break, end;
  end
  