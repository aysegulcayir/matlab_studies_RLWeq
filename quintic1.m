%quitic b-spline collocation method
clc; clear all;close all;
alfa=0.5;mu=0; h=20; ti=20; a=0; b=9000; tson=9600;ro=264; x0=2000;
x=a:h:b;
n=(b-a)/h;
%t=0 aný
for i=1:n+1%matlabda 0 indisi yok o yüzden bir fazla
    u(i)=10*exp(-((x(i)-x0)^2)/(2*ro^2));
    f(i)=u(i);
end
ea1=1;ea2=26;ea3=66;ea4=26;ea5=1; eb1=5/h;eb2=50/h; eb3=0;eb4=-50/h;eb5=-5/h; 
ec1=20/(h^2);ec2=40/(h^2);ec3=-120/(h^2);ec4=40/(h^2);ec5=20/(h^2);
A=zeros(n+1,n+1);
A(1,1)=ea2;A(1,2)=ea3-eb3/eb1;
for i=2:n
    A(i,i-2)=ea1;A(i,i-1)=ea2;A(i,i)=ea3;A(i,i+1)=ea4;A(i,i+2)=ea5;
end
A(n+1,n)=ea1-eb1/eb3;A(n+1,n+1)=ea2;
gc=A\f';
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


A=zeros(n+1,n+1);
for t=ti:ti:tson
  %theta=0;
   % theta=1;
    theta=0.5;
    b1=ea1+theta*ti*alfa*eb1-theta*ti*mu*ec1;
    b2=ea2+theta*ti*alfa*eb2-theta*mu*ti*ec2;
    b3=ea3+theta*ti*alfa*eb3-theta*ti*mu*ec3;
    b4=ea1-(1-theta)*ti*alfa*eb1+(1-theta)*ti*mu*ec1;
    b5=ea2-(1-theta)*ti*alfa*eb2+(1-theta)*mu*ti*ec2;
    b6=ea3-(1-theta)*ti*alfa*eb3+(1-theta)*ti*mu*ec3;
    A(1,1)=b2-ea2*b1/ea1;A(1,2)=b3-ea3*b1/ea1;
    b(1)=b4*g(1)+b5*g(2)+b6*g(3);
    
    for i=2:n
        A(i,i-1)=b1; A(i,i)=b2;A(i,i+1)=b3;
        b(i)=b4*g(i)+b5*g(i+1)+b6*g(i+2);
    end
    A(n+1,n)=b1-ea1*b3/ea3;A(n+1,n+1)=b2-ea2*b3/ea3;    
    b(n+1)=b4*g(n+1)+b5*g(n+2)+b6*g(n+3);    
    
    gc=A\b';
    for i=1:n+1
        g1(i+1)=gc(i);
    end
    %sýnýr þartlarýndan
    g1(1)=-(ea2*g1(2)+ea3*g1(3))/ea1;
    g1(n+3)=-(ea1*g1(n+1)+ea2*g(n+2))/ea3;
    
   g=g1;  %n+1. zamandakileri bulmak için n.zamana ihtiyaç var

% for i=1:n+1
%         uu(i)=10*exp(-((x(i)-x0-alfa*t)^2)/(2*ro^2));
%         uy(i)=ea1*g1(i)+ea2*g1(i+1)+ea3*g1(i+2);
%         hata(i)=abs(uu(i)-uy(i)); 
%      %  fprintf('x=%.5f,utam=%.10f,uyak=%.10f,hata=%.10f\n',x(i),uu(i),uy(i),hata(i));
%  end
%     fprintf('t=%.5f,maxhata=%.10f\n',t,max(hata));
end
    
 for i=1:n+1
        uu(i)=10*exp(-((x(i)-x0-alfa*tson)^2)/(2*ro^2));
        uy(i)=ea1*g1(i)+ea2*g1(i+1)+ea3*g1(i+2);
        hata(i)=abs(uu(i)-uy(i)); 
      % fprintf('x=%.5f,utam=%.10f,uyak=%.10f,hata=%.10f\n',x(i),uu(i),uy(i),hata(i));
 end
     fprintf('t=%.5f,maxhata=%.10f\n',t,max(hata));
     plot(x,hata);   