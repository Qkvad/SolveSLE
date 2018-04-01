function [x,res]=GaussSeidel(A,b,x0,tol,kmax)
  n=max(size(A));
  normab=norm(b);
  normaA=norm(A,'fro');
  flag = 0;
  
  if normab==0
    x=zeros(n,1);
    k=1;
    res(1)=0;
    return
  end
  
  k=1;
  x=x0;
  res(1)=norm(b-A*x0)/normab;
  
  if(res(1)<tol) return, end
  
  for k=2:kmax
    for i=1:n
      x(i)=(b(i)-(A(i,1:n)*x(1:n)-A(i,i)*x(i)))/A(i,i);
    end
    res(k)=norm(b-A*x)/normab;
    if(res(k)<=tol), break, end
  end
  if(res(k)>tol) flag=1;end;
