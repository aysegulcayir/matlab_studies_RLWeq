clc; clear all;close all;

za(1)=100;za(2)=50;za(3)=25;za(4)=20;za(5)=10;za(6)=5;za(7)=2.5;za(8)=2;za(9)=1;
for ij=1:5
    ti=za(ij);h=za(ij);
    alfa=0.5;mu=0; a=0; b=9000; tson=9600;ro=264; x0=2000;
%h=20; ti=20; 
x=a:h:b;
n=(b-a)/h;
%t=0 aný
u=zeros(n+1,1);g=zeros(n+2,1);yh=zeros(n+1,1);er=zeros(n+1,1);u1=zeros(n+1,1);g1=zeros(n+2,1);
for i=1:n+1
    u(i)=10*exp(-((x(i)-x0)^2)/(2*ro^2));
end
uxa=0;
 ea1=1;ea2=1;
eb1=-2/h;eb2=2/h;

g(1)=(eb2*u(1)-ea2*uxa)/(eb2*ea1-ea2*eb1);
for i=1:n+1
    g(i+1)=(u(i)-ea1*g(i))/ea2;    
end

for i=1:n+1
   yh(i)=ea1*g(i)+ea2*g(i+1);    
   er(i)=abs(yh(i)-u(i));    
   
end
%disp(max(er));

syms k;
% quadratic B-spline þekil fonksiyonlarý
 ex(1)=((h-k)/h)^2;
 ex(2)=((2*h-k)^2-3*(h-k)^2)/h^2;
 ex(3)=((3*h-k)^2-3*(2*h-k)^2+3*(h-k)^2)/h^2;
a=zeros(3,3);b=zeros(3,3);c=zeros(3,3);


 for i=1:3
    for j=1:3
       a(i,j)=int(ex(i)*ex(j),k,0,h);
       b(i,j)=int(ex(i)*diff(ex(j),k),k,0,h);
       c(i,j)=int(diff(ex(i),k)*diff(ex(j),k),k,0,h);
    end
 end
 sol=zeros(n+2,n+2); sag=zeros(n+2,n+2);f=zeros(n+2,1);
 for t=ti:ti:tson
  
%  sol(1,1)=a(1,1)+ti*alfa*(b(1,1))/2+ti*mu*(c(1,1))/2;
%  sol(1,2)=a(1,2)+ti*alfa*(b(1,2))/2+ti*mu*(c(1,2))/2; 
%  sol(1,3)=a(1,3)+ti*alfa*(b(1,3))/2+ti*mu*(c(1,3))/2;
%  
%  sag(1,1)=a(1,1)-ti*alfa*(b(1,1))/2-ti*mu*(c(1,1))/2;
%  sag(1,2)=a(1,2)-ti*alfa*(b(1,2))/2-ti*mu*(c(1,2))/2; 
%  sag(1,3)=a(1,3)-ti*alfa*(b(1,3))/2-ti*mu*(c(1,3))/2;
%  
%  f(1)=g(1)*sag(1,1)+g(2)*sag(1,2)+g(3)*sag(1,3);

sol(1,1)=ea1;sol(1,2)=ea2; f(1)=0;
 
 sol(2,1)=a(2,1)+ti*alfa*(b(2,1))/2+ti*mu*(c(2,1))/2;
 sol(2,2)=a(2,2)+a(1,1)+ti*alfa*(b(2,2)+b(1,1))/2+ti*mu*(c(2,2)+c(1,1))/2; 
 sol(2,3)=a(2,3)+a(1,2)+ti*alfa*(b(2,3)+b(1,2))/2+ti*mu*(c(2,3)+c(1,2))/2;
 sol(2,4)=a(1,3)+ti*alfa*(b(1,3))/2+ti*mu*(c(1,3))/2;
 
 sag(2,1)=a(2,1)-ti*alfa*(b(2,1))/2-ti*mu*(c(2,1))/2;
 sag(2,2)=a(2,2)+a(1,1)-ti*alfa*(b(2,2)+b(1,1))/2-ti*mu*(c(2,2)+c(1,1))/2; 
 sag(2,3)=a(2,3)+a(1,2)-ti*alfa*(b(2,3)+b(1,2))/2-ti*mu*(c(2,3)+c(1,2))/2;
 sag(2,4)=a(1,3)-ti*alfa*(b(1,3))/2-ti*mu*(c(1,3))/2;
 
 f(2)=g(1)*sag(2,1)+g(2)*sag(2,2)+g(3)*sag(2,3)+g(4)*sag(2,4);
 for i=3:n
     sol(i,i-2)=a(3,1)+ti*alfa*(b(3,1))/2+ti*mu*(c(3,1))/2;
     sol(i,i-1)=a(3,2)+a(2,1)+ti*alfa*(b(3,2)+b(2,1))/2+ti*mu*(c(3,2)+c(2,1))/2;
     sol(i,i)=a(3,3)+a(2,2)+a(1,1)+ti*alfa*(b(3,3)+b(2,2)+b(1,1))/2+ti*mu*(c(3,3)+c(2,2)+c(1,1))/2;
     sol(i,i+1)=a(2,3)+a(1,2)+ti*alfa*(b(2,3)+b(1,2))/2+ti*mu*(c(2,3)+c(1,2))/2;
     sol(i,i+2)=a(1,3)+ti*alfa*(b(1,3))/2+ti*mu*(c(1,3))/2;
     
     sag(i,i-2)=a(3,1)-ti*alfa*(b(3,1))/2-ti*mu*(c(3,1))/2;
     sag(i,i-1)=a(3,2)+a(2,1)-ti*alfa*(b(3,2)+b(2,1))/2-ti*mu*(c(3,2)+c(2,1))/2;
     sag(i,i)=a(3,3)+a(2,2)+a(1,1)-ti*alfa*(b(3,3)+b(2,2)+b(1,1))/2-ti*mu*(c(3,3)+c(2,2)+c(1,1))/2;
     sag(i,i+1)=a(2,3)+a(1,2)-ti*alfa*(b(2,3)+b(1,2))/2-ti*mu*(c(2,3)+c(1,2))/2;
     sag(i,i+2)=a(1,3)-ti*alfa*(b(1,3))/2-ti*mu*(c(1,3))/2;
     f(i)=g(i-2)*sag(i,i-2)+g(i-1)*sag(i,i-1)+g(i)*sag(i,i)+g(i+1)*sag(i,i+1)+g(i+2)*sag(i,i+2);
     
 end 
 
     sol(n+1,n-1)=a(3,1)+ti*alfa*(b(3,1))/2+ti*mu*(c(3,1))/2;
     sol(n+1,n)=a(3,2)+a(2,1)+ti*alfa*(b(3,2)+b(2,1))/2+ti*mu*(c(3,2)+c(2,1))/2;
     sol(n+1,n+1)=a(3,3)+a(2,2)+ti*alfa*(b(3,3)+b(2,2))/2+ti*mu*(c(3,3)+c(2,2))/2;
     sol(n+1,n+2)=a(2,3)+ti*alfa*(b(2,3))/2+ti*mu*(c(2,3))/2;
     
     
     sag(n+1,n-1)=a(3,1)-ti*alfa*(b(3,1))/2-ti*mu*(c(3,1))/2;
     sag(n+1,n)=a(3,2)+a(2,1)-ti*alfa*(b(3,2)+b(2,1))/2-ti*mu*(c(3,2)+c(2,1))/2;
     sag(n+1,n+1)=a(3,3)+a(2,2)-ti*alfa*(b(3,3)+b(2,2))/2-ti*mu*(c(3,3)+c(2,2))/2;
     sag(n+1,n+2)=a(2,3)-ti*alfa*(b(2,3))/2-ti*mu*(c(2,3))/2;
     f(n+1)=g(n-1)*sag(n+1,n-1)+g(n)*sag(n+1,n)+g(n+1)*sag(n+1,n+1)+g(n+2)*sag(n+1,n+2);    
     
%      sol(n+2,n)=a(3,1)+ti*alfa*(b(3,1))/2+ti*mu*(c(3,1))/2;
%      sol(n+2,n+1)=a(3,2)+ti*alfa*(b(3,2))/2+ti*mu*(c(3,2))/2;
%      sol(n+2,n+2)=a(3,3)+ti*alfa*(b(3,3))/2+ti*mu*(c(3,3))/2;
%      
%      
%      
%     sag(n+2,n)=a(3,1)-ti*alfa*(b(3,1))/2-ti*mu*(c(3,1))/2;
%      sag(n+2,n+1)=a(3,2)-ti*alfa*(b(3,2))/2-ti*mu*(c(3,2))/2;
%      sag(n+2,n+2)=a(3,3)-ti*alfa*(b(3,3))/2-ti*mu*(c(3,3))/2;
%      
%      
%      f(n+2)=g(n)*sag(n+2,n)+g(n+1)*sag(n+2,n+1)+g(n+2)*sag(n+2,n+2);   

sol(n+2,n+1)=ea1;sol(n+2,n+2)=ea2; f(n+2)=0;
 A=sparse(sol);
   g1=A\f;
   g=g1;
   
  
 end
 for i=1:n+1
   yh(i)=ea1*g1(i)+ea1*g1(i+1); 
   u1(i)=10*exp(-((x(i)-x0-alfa*t)^2)/(2*ro^2));
   er(i)=abs(yh(i)-u1(i));   
 end
   ymax(ij)=max(er);
fprintf('t=%.5f,h=%.5f,ti=%.5f,maxhata=%.10f\n',t,h,ti,max(er));
end
 figure(1); plot(x,u1);
 figure(2); plot(x,yh);
 figure(3); plot(x,er);
 
for i=2:ij
    rate=log(ymax(i)/ymax(i-1))/log(za(i)/za(i-1));
    fprintf('order:%.15f\n',rate);
end 
 
 
 