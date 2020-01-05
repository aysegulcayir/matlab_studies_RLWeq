%quinitc b-spline collocation method

clc; clear all;close all;
zaman=cputime;
mu=1; h=0.03; ti=0.05; a=-20; b=40; tson=80; x0=10;
c=0.1;
x=a:h:b;
n=(b-a)/h;
%t=0 aný
k=sqrt(1/(4*mu));
for i=1:n+1%matlabda 0 indisi yok o yüzden bir fazla
    u(i)=3*c*sech(k*(x(i)-x0))^2;
    f(i)=u(i);
end
uxa=0;uxxa=0;uxb=0;uxxb=0;
ea1=1;ea2=26;ea3=66;ea4=26;ea5=1;
eb1=-5/h; eb2=-50/h; eb3=0; eb4=50/h;eb5=5/h;
ec1=20/(h^2);ec2=40/(h^2);ec3=-120/(h^2);ec4=40/(h^2);ec5=20/(h^2);
A=zeros(n+5,n+5);ff=zeros(n+5,1);

%1. denklem ux(a,0)=0 dan
A(1,1)=eb1;A(1,2)=eb2;A(1,3)=eb3;A(1,4)=eb4;A(1,5)=eb5;ff(1)=0;
%2. denklem uxx(a,0)=0 dan
A(2,1)=ec1;A(2,2)=ec2;A(2,3)=ec3;A(2,4)=ec4;A(2,5)=ec5;ff(2)=0;
for i=3:n+3
    A(i,i-2)=ea1;A(i,i-1)=ea2;A(i,i)=ea3;A(i,i+1)=ea4;A(i,i+2)=ea5;
    ff(i)=u(i-2);
end

%n+4. denklem ux(b,0)=0 
A(n+4,n+1)=eb1;A(n+4,n+2)=eb2;A(n+4,n+3)=eb3;A(n+4,n+4)=eb4;
A(n+4,n+5)=eb5;ff(n+4)=0;
%n+5. denklem uxx(b,0)=0 
A(n+5,n+1)=ec1;A(n+5,n+2)=ec2;A(n+5,n+3)=ec3;A(n+5,n+4)=ec4;
A(n+5,n+5)=ec5;ff(n+5)=0;

AA= sparse(A);
 g=AA\ff;
    
% b3=(-eb5*ec2+eb2*ec5)/(eb1*ec2-ec1*eb2);b2=(-eb4*ec2+eb2*ec4)/(eb1*ec2-ec1*eb2);
% b1=(-eb3*ec2+eb2*ec3)/(eb1*ec2-ec1*eb2);bd1=(-eb2*uxxa+uxa*ec2)/(eb1*ec2-ec1*eb2);
% b6=-(-ec1*eb5+eb1*ec5)/(eb1*ec2-ec1*eb2);b5=-(eb1*ec4-ec1*eb4)/(eb1*ec2-ec1*eb2);
% b4=-(eb1*ec3-ec1*eb3)/(eb1*ec2-ec1*eb2);bd2=-(-eb1*uxxa+ec1*uxa)/(eb1*ec2-ec1*eb2);
% A(1,1)=ea3+ea1*b1+ea2*b4;
% A(1,2)=ea4+ea1*b2+ea2*b5;A(1,3)=ea5+ea1*b3+ea2*b6;
% f(1)=f(1)-ea1*bd1-ea2*bd2;
% A(2,1)=ea2+ea1*b4;A(2,2)=ea3+ea1*b5;A(2,3)=ea4+ea1*b6;
% A(2,4)=ea5;
% f(2)=f(2)-ea1*bd2;    
% for i=3:n-1
%     A(i,i-2)=ea1;A(i,i-1)=ea2;A(i,i)=ea3;A(i,i+1)=ea4;A(i,i+2)=ea5;
% end
% s1=-(-ec1*eb5+eb1*ec5)/(eb4*ec5-ec4*eb5);s2=-(-eb5*ec2+eb2*ec5)/(eb4*ec5-ec4*eb5);s3=-(-eb5*ec3+eb3*ec5)/(eb4*ec5-ec4*eb5);
% sd1=-(-uxb*ec5+eb5*uxxb)/(eb4*ec5-ec4*eb5);
% s4=(eb1*ec4-ec1*eb4)/(eb4*ec5-ec4*eb5);s5=(-eb4*ec2+eb2*ec4)/(eb4*ec5-ec4*eb5);s6=(-eb4*ec3+ec4*eb3)/(eb4*ec5-ec4*eb5);
% sd2=(-ec4*uxb+eb4*uxxb)/(eb4*ec5-ec4*eb5);
% A(n,n-2)=ea1;A(n,n-1)=ea2+ea5*s1;
% A(n,n)=ea3+ea5*s2;A(n,n+1)=ea4+ea5*s3;
% A(n+1,n-1)=ea1+s1*ea4+s4*ea5;A(n+1,n)=ea2+s2*ea4+s5*ea5;
% A(n+1,n+1)=ea3+s3*ea4+s6*ea5;
% f(n)=f(n)-ea5*sd1;
% f(n+1)=f(n+1)-ea4*sd1-ea5*sd2;
%  AA= sparse(A);
%  gc=AA\f';
% for i=1:n+1
%     g(i+2)=gc(i);
% end
%  
%  g(1)=b1*g(3)+b2*g(4)+b3*g(5);
%  g(2)=b4*g(3)+b5*g(4)+b6*g(5);
%  g(n+4)=s1*g(n+1)+s2*g(n+2)+s3*g(n+3);
%  g(n+5)=s4*g(n+1)+s5*g(n+2)+s6*g(n+3);

 for i=1:n+1
    uy(i)=ea1*g(i)+ea2*g(i+1)+ea3*g(i+2)+ea4*g(i+3)+ea5*g(i+4);
    hata(i)=abs(uy(i)-u(i));
end
% figure(1); 
% plot(x,hata);
% figure(2);
% plot(x,u);
% 
% 
A=zeros(n+5,n+5);
for t=ti:ti:tson
%   uu=ea1*g(1)+ea2*g(2)+ea3*g(3)+ea4*g(4)+ea5*g(5);
%   ux=eb1*g(1)+eb2*g(2)+eb3*g(3)+eb4*g(4)+eb5*g(5);
%     kat1=(1+ti*ux/2)*ea1+ti*uu*eb1/2-mu*ec1;
%     kat2=(1+ti*ux/2)*ea2+ti*uu*eb2/2-mu*ec2;
%     kat3=(1+ti*ux/2)*ea3+ti*uu*eb3/2-mu*ec3;
%     kat4=(1+ti*ux/2)*ea4+ti*uu*eb4/2-mu*ec4;
%     kat5=(1+ti*ux/2)*ea5+ti*uu*eb5/2-mu*ec5;
%     
%     kat6=ea1-mu*ec1;
%     kat7=ea2-mu*ec2;
%     kat8=ea3-mu*ec3;
%     kat9=ea4-mu*ec4;
%     kat10=ea5-mu*ec5;
%     
% b3=(-eb5*ec2+eb2*ec5)/(eb1*ec2-ec1*eb2);b2=(-eb4*ec2+eb2*ec4)/(eb1*ec2-ec1*eb2);
% b1=(-eb3*ec2+eb2*ec3)/(eb1*ec2-ec1*eb2);bd1=(-eb2*uxxa+uxa*ec2)/(eb1*ec2-ec1*eb2);
% b6=-(-ec1*eb5+eb1*ec5)/(eb1*ec2-ec1*eb2);b5=-(eb1*ec4-ec1*eb4)/(eb1*ec2-ec1*eb2);
% b4=-(eb1*ec3-ec1*eb3)/(eb1*ec2-ec1*eb2);bd2=-(-eb1*uxxa+ec1*uxa)/(eb1*ec2-ec1*eb2);
% 
% A(1,1)=kat3+ea1*b1+ea2*b4;
% A(1,2)=kat4+ea1*b2+ea2*b5;A(1,3)=kat5+ea1*b3+ea2*b6;
% for i=1:1
%     f(i)=kat6*g(i)+kat7*g(i+1)+kat8*g(i+2)+kat9*g(i+3)+kat10*g(i+4);
% end
%   uu=ea1*g(2)+ea2*g(3)+ea3*g(4)+ea4*g(5)+ea5*g(6);
%   ux=eb1*g(2)+eb2*g(3)+eb3*g(4)+eb4*g(5)+eb5*g(6);
%     kat1=(1+ti*ux/2)*ea1+ti*uu*eb1/2-mu*ec1;
%     kat2=(1+ti*ux/2)*ea2+ti*uu*eb2/2-mu*ec2;
%     kat3=(1+ti*ux/2)*ea3+ti*uu*eb3/2-mu*ec3;
%     kat4=(1+ti*ux/2)*ea4+ti*uu*eb4/2-mu*ec4;
%     kat5=(1+ti*ux/2)*ea5+ti*uu*eb5/2-mu*ec5;
%     
%     kat6=ea1-mu*ec1;
%     kat7=ea2-mu*ec2;
%     kat8=ea3-mu*ec3;
%     kat9=ea4-mu*ec4;
%     kat10=ea5-mu*ec5;
% A(2,1)=kat2+ea1*b4;A(2,2)=kat3+ea1*b5;A(2,3)=kat4+ea1*b6;
% A(2,4)=kat5;
% for i=2:2
%     f(i)=kat6*g(i)+kat7*g(i+1)+kat8*g(i+2)+kat9*g(i+3)+kat10*g(i+4);
% end
%1. denklem u(a,0)=0 dan
A(1,1)=ea1;A(1,2)=ea2;A(1,3)=ea3;A(1,4)=ea4;A(1,5)=ea5;f(1)=0;
%2. denklem ux(a,0)=0 dan
A(2,1)=eb1;A(2,2)=eb2;A(2,3)=eb3;A(2,4)=eb4;A(2,5)=eb5;f(2)=0;
for i=3:n+3
  uu=ea1*g(i-2)+ea2*g(i-1)+ea3*g(i)+ea4*g(i+1)+ea5*g(i+2);
  ux=eb1*g(i-2)+eb2*g(i-1)+eb3*g(i)+eb4*g(i+1)+eb5*g(i+2);
  
    kat1=(1+ti*ux/2)*ea1+ti*uu*eb1/2-mu*ec1;
    kat2=(1+ti*ux/2)*ea2+ti*uu*eb2/2-mu*ec2;
    kat3=(1+ti*ux/2)*ea3+ti*uu*eb3/2-mu*ec3;
    kat4=(1+ti*ux/2)*ea4+ti*uu*eb4/2-mu*ec4;
    kat5=(1+ti*ux/2)*ea5+ti*uu*eb5/2-mu*ec5;
    
    kat6=ea1-mu*ec1;
    kat7=ea2-mu*ec2;
    kat8=ea3-mu*ec3;
    kat9=ea4-mu*ec4;
    kat10=ea5-mu*ec5;
    A(i,i-2)=kat1;A(i,i-1)=kat2;A(i,i)=kat3;
    A(i,i+1)=kat4;A(i,i+2)=kat5; 
    f(i)=kat6*g(i-2)+kat7*g(i-1)+kat8*g(i)+kat9*g(i+1)+kat10*g(i+2);
end

%n+4. denklem u(b,0)=0 
A(n+4,n+1)=ea1;A(n+4,n+2)=ea2;A(n+4,n+3)=ea3;A(n+4,n+4)=ea4;
A(n+4,n+5)=ea5;f(n+4)=0;
%n+5. denklem ux(b,0)=0 
A(n+5,n+1)=eb1;A(n+5,n+2)=eb2;A(n+5,n+3)=eb3;A(n+5,n+4)=eb4;
A(n+5,n+5)=eb5;f(n+5)=0;

AA= sparse(A);
 g1=AA\f';

% s1=-(-ec1*eb5+eb1*ec5)/(eb4*ec5-ec4*eb5);s2=-(-eb5*ec2+eb2*ec5)/(eb4*ec5-ec4*eb5);s3=-(-eb5*ec3+eb3*ec5)/(eb4*ec5-ec4*eb5);
% sd1=-(-uxb*ec5+eb5*uxxb)/(eb4*ec5-ec4*eb5);
% s4=(eb1*ec4-ec1*eb4)/(eb4*ec5-ec4*eb5);s5=(-eb4*ec2+eb2*ec4)/(eb4*ec5-ec4*eb5);s6=(-eb4*ec3+ec4*eb3)/(eb4*ec5-ec4*eb5);
% sd2=(-ec4*uxb+eb4*uxxb)/(eb4*ec5-ec4*eb5);
% 
% uu=ea1*g(n)+ea2*g(n+1)+ea3*g(n+2)+ea4*g(n+3)+ea5*g(n+4);
%   ux=eb1*g(n)+eb2*g(n+1)+eb3*g(n+2)+eb4*g(n+3)+eb5*g(n+4);
%     kat1=(1+ti*ux/2)*ea1+ti*uu*eb1/2-mu*ec1;
%     kat2=(1+ti*ux/2)*ea2+ti*uu*eb2/2-mu*ec2;
%     kat3=(1+ti*ux/2)*ea3+ti*uu*eb3/2-mu*ec3;
%     kat4=(1+ti*ux/2)*ea4+ti*uu*eb4/2-mu*ec4;
%     kat5=(1+ti*ux/2)*ea5+ti*uu*eb5/2-mu*ec5;
%     
%     kat6=ea1-mu*ec1;
%     kat7=ea2-mu*ec2;
%     kat8=ea3-mu*ec3;
%     kat9=ea4-mu*ec4;
%     kat10=ea5-mu*ec5;
% A(n,n-2)=kat1;A(n,n-1)=kat2+ea5*s1;
% A(n,n)=kat3+ea5*s2;A(n,n+1)=kat4+ea5*s3;
% for i=n:n
%     f(i)=kat6*g(i)+kat7*g(i+1)+kat8*g(i+2)+kat9*g(i+3)+kat10*g(i+4);
% end
% uu=ea1*g(n+1)+ea2*g(n+2)+ea3*g(n+3)+ea4*g(n+4)+ea5*g(n+5);
%   ux=eb1*g(n+1)+eb2*g(n+2)+eb3*g(n+3)+eb4*g(n+4)+eb5*g(n+5);
%     kat1=(1+ti*ux/2)*ea1+ti*uu*eb1/2-mu*ec1;
%     kat2=(1+ti*ux/2)*ea2+ti*uu*eb2/2-mu*ec2;
%     kat3=(1+ti*ux/2)*ea3+ti*uu*eb3/2-mu*ec3;
%     kat4=(1+ti*ux/2)*ea4+ti*uu*eb4/2-mu*ec4;
%     kat5=(1+ti*ux/2)*ea5+ti*uu*eb5/2-mu*ec5;
%     
%     kat6=ea1-mu*ec1;
%     kat7=ea2-mu*ec2;
%     kat8=ea3-mu*ec3;
%     kat9=ea4-mu*ec4;
%     kat10=ea5-mu*ec5;
% A(n+1,n-1)=kat1+s1*ea4+s4*ea5;A(n+1,n)=kat2+s2*ea4+s5*ea5;
% A(n+1,n+1)=kat3+s3*ea4+s6*ea5;
% for i=n+1:n+1
%     f(i)=kat6*g(i)+kat7*g(i+1)+kat8*g(i+2)+kat9*g(i+3)+kat10*g(i+4);
% end
% AA= sparse(A);
%  gc=AA\f';
% for i=1:n+1
%     g1(i+2)=gc(i);
% end
%  
%  g1(1)=b1*g1(3)+b2*g1(4)+b3*g1(5);
%  g1(2)=b4*g1(3)+b5*g1(4)+b6*g1(5);
%  g1(n+4)=s1*g1(n+1)+s2*g1(n+2)+s3*g1(n+3);
%  g1(n+5)=s4*g1(n+1)+s5*g1(n+2)+s6*g1(n+3);

    
   g=g1;  %n+1. zamandakileri bulmak için n.zamana ihtiyaç var

for i=1:n+1
        uu(i)=3*c*sech(k*(x(i)-x0-c*t))^2;
        uy(i)=ea1*g1(i)+ea2*g1(i+1)+ea3*g1(i+2)+ea4*g1(i+3)+ea5*g1(i+4);
        hata(i)=abs(uu(i)-uy(i)); 
     %  fprintf('x=%.5f,utam=%.10f,uyak=%.10f,hata=%.10f\n',x(i),uu(i),uy(i),hata(i));
 end
    fprintf('t=%.5f,maxhata=%.10f\n',t,max(hata));
end
    
 for i=1:n+1
        uu(i)=3*c*sech(k*(x(i)-x0-c*t))^2;
        uy(i)=ea1*g1(i)+ea2*g1(i+1)+ea3*g1(i+2)+ea4*g1(i+3)+ea5*g1(i+4);
        hata(i)=abs(uu(i)-uy(i)); 
      fprintf('x=%.5f,utam=%.10f,uyak=%.10f,hata=%.10f\n',x(i),uu(i),uy(i),hata(i));
 end
     fprintf('t=%.5f,maxhata=%.10f\n',t,max(hata));
     plot(x,hata);   
     disp(abs(cputime-zaman));
    
    
    


