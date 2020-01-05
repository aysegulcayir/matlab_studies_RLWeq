%cubic b-spline collocation method for ew equation
%iç iterasyon
clc; clear all;close all;
zaman=cputime;
mu=1; h=0.03; ti=0.05; a=0; b=30; tson=80; x0=10;
c=0.1;
x=a:h:b;
n=(b-a)/h;
%t=0 aný
k=sqrt(1/(4*mu));
u=zeros(1,n+1);f=zeros(1,n+1);g=zeros(1,n+3);g1=zeros(1,n+3);
b=zeros(1,n+1);
for i=1:n+1%matlabda 0 indisi yok o yüzden bir fazla
    u(i)=3*c*sech(k*(x(i)-x0))^2;
    f(i)=u(i);
end
ea1=1;ea2=4;ea3=1; eb1=-3/h; eb2=0; eb3=3/h; 
ec1=6/(h^2);ec2=-12/(h^2);ec3=6/(h^2);
A=zeros(n+1,n+1);
A(1,1)=ea2;A(1,2)=ea3-eb3/eb1;
for i=2:n
    A(i,i-1)=ea1;A(i,i)=ea2;A(i,i+1)=ea3;
end
A(n+1,n)=ea1-eb1/eb3;A(n+1,n+1)=ea2;
AA=sparse(A);
gc=AA\f';
for i=1:n+1
    g(i+1)=gc(i);
end

g(1)=-eb3*g(3)/eb1;
g(n+3)=-eb1*g(n+1)/eb3;
for i=1:n+1
    uy(i)=ea1*g(i)+ea2*g(i+1)+ea3*g(i+2);
    hata(i)=abs(uy(i)-u(i));
end
% figure(1); 
%plot(x,hata);
% figure(2);
% plot(x,u);


gg=g;%g=t=0 anýndaki deltalar
for t=ti:ti:tson
    for iterasyon=1:3
  A=zeros(n+1,n+1);
    uu=ea1*gg(1)+ea2*gg(2)+ea3*gg(3);
    uk=ea1*g(1)+ea2*g(2)+ea3*g(3);
    b1=ea1+ti*uu*eb1/2-mu*ec1;
    b2=ea2+ti*uu*eb2/2-mu*ec2;
    b3=ea3+ti*uu*eb3/2-mu*ec3;
    
    b4=ea1-ti*uk*eb1/2-mu*ec1;     b5=ea2-ti*uk*eb2/2-mu*ec2; 
    b6=ea3-ti*uk*eb3/2-mu*ec3;
    A(1,1)=b2-ea2*b1/ea1;A(1,2)=b3-ea3*b1/ea1;
    b(1)=b4*g(1)+b5*g(2)+b6*g(3);
    
    for i=2:n
    uu=ea1*gg(i)+ea2*gg(i+1)+ea3*gg(i+2);
    uk=ea1*g(i)+ea2*g(i+1)+ea3*g(i+2);
    b1=ea1+ti*uu*eb1/2-mu*ec1;
    b2=ea2+ti*uu*eb2/2-mu*ec2;
    b3=ea3+ti*uu*eb3/2-mu*ec3;
    
    b4=ea1-ti*uk*eb1/2-mu*ec1;     b5=ea2-ti*uk*eb2/2-mu*ec2; 
    b6=ea3-ti*uk*eb3/2-mu*ec3;
        A(i,i-1)=b1; A(i,i)=b2;A(i,i+1)=b3;
        b(i)=b4*g(i)+b5*g(i+1)+b6*g(i+2);
    end
   uu=ea1*gg(n+1)+ea2*gg(n+2)+ea3*gg(n+3);
    uk=ea1*g(n+1)+ea2*g(n+2)+ea3*g(n+3);
    b1=ea1+ti*uu*eb1/2-mu*ec1;
    b2=ea2+ti*uu*eb2/2-mu*ec2;
    b3=ea3+ti*uu*eb3/2-mu*ec3;
    
    b4=ea1-ti*uk*eb1/2-mu*ec1;     b5=ea2-ti*uk*eb2/2-mu*ec2; 
    b6=ea3-ti*uk*eb3/2-mu*ec3;
    A(n+1,n)=b1-ea1*b3/ea3;A(n+1,n+1)=b2-ea2*b3/ea3;    
    b(n+1)=b4*g(n+1)+b5*g(n+2)+b6*g(n+3);    
    AA=sparse(A);
    gc=AA\b';
    for i=1:n+1
        g1(i+1)=gc(i);
    end
    %sýnýr þartlarýndan
    g1(1)=-(ea2*g1(2)+ea3*g1(3))/ea1;
    g1(n+3)=-(ea1*g1(n+1)+ea2*g(n+2))/ea3;
    gg=g1;
    end
   g=g1;  %n+1. zamandakileri bulmak için n.zamana ihtiyaç var

for i=1:n+1
        uu(i)=3*c*sech(k*(x(i)-x0-c*t))^2;
        uy(i)=ea1*g1(i)+ea2*g1(i+1)+ea3*g1(i+2);
        hata(i)=abs(uu(i)-uy(i)); 
     %  fprintf('x=%.5f,utam=%.10f,uyak=%.10f,hata=%.10f\n',x(i),uu(i),uy(i),hata(i));
 end
    fprintf('t=%.5f,maxhata=%.10f\n',t,max(hata));
end
    
 for i=1:n+1
        uu(i)=3*c*sech(k*(x(i)-x0-c*t))^2;
        uy(i)=ea1*g1(i)+ea2*g1(i+1)+ea3*g1(i+2);
        hata(i)=abs(uu(i)-uy(i)); 
      % fprintf('x=%.5f,utam=%.10f,uyak=%.10f,hata=%.10f\n',x(i),uu(i),uy(i),hata(i));
 end
     fprintf('t=%.5f,maxhata=%.10f\n',t,max(hata));
     plot(x,hata);   
    disp(cputime-zaman);
    
    


