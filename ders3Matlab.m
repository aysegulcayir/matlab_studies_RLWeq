% clc; clear all;
% f=@(x)x+2;%f(x)=x+2;
% disp(f(2));
% g=@(x,y) x*y+2;
% disp(g(2,3));
%ödev sorusu
% clc;clear all;
% f=@(x)x*exp(-x)+2; a=-5;b=5;h=0.001;
% x=a:h:b;%n+1 nokta
% n=(b-a)/h;% n aralýk
% for i=1:n
%     if f(x(i))*f(x(i+1))<0
%         a=x(i); b=x(i+1); break;
%     end
% end
% fprintf('kokun olduðu aralýk [%.12f,%.12f]\n',a,b);
% y(1)=(a+b)/2;%newton icin baslangýc deðeri
% %ikiye bolme
% hata=1;tol=10^(-10);sayac=0;
% while(hata>tol)
%     c=(a+b)/2;
%     if f(a)*f(c)<0 
%         b=c;
%     else
%         a=c;
%     end
%     hata=b-a;     sayac=sayac+1;
% end
% fprintf('ikiye bolme: %g adým sonunda aralýk [%.12f,%.12f]\n',sayac,a,b); 
% fx=@(x)exp(-x)-x*exp(-x);%newton rtaphson icin f nin türevi gerekliydi 
% hata=1;i=1;
% while(hata>tol)
%     y(i+1)=y(i)-f(y(i))/fx(y(i));
%     hata=abs(y(i+1)-y(i));
%    % fprintf('newton: %g adým sonunda kok=%.12f]\n',i,y(i+1)); 
%     i=i+1;
% end
% fprintf('newton: %g adým sonunda kok=%.12f]\n',i-1,y(i));

%*************************yamuklar
% clc; clear all;
% n=input('n='); h=(1-0)/n; x=0:h:1; f=@(x)x^2; sonuc=0;
% for i=1:n
%     sonuc=sonuc+h*(f(x(i))+f(x(i+1)))/2;
% end
% fprintf('n=%g icin sonuc=%.15f\n',n,sonuc);
% clc; clear all;
% for j=0:5

% n=10^j; h=(1-0)/n; x=0:h:1; f=@(x)x^2; sonuc=0;
% for i=1:n
%     sonuc=sonuc+h*(f(x(i))+f(x(i+1)))/2;
% end
% fprintf('n=%g icin sonuc=%.15f\n',n,sonuc);
% end
%%n=100,200,300,...,1000 için bir önceki soruyu simpsonla yapalým

% clc; clear all;
% for n=100:100:1000
% h=(1-0)/n; x=0:h:1; f=@(x)x^2; sonuc=0;
% for i=1:2:n-1
%     sonuc=sonuc+h*(f(x(i))+4*f(x(i+1))+f(x(i+2)))/3;
% end
% fprintf('n=%g icin sonuc=%.15f\n',n,sonuc);
% end
%Soru: 

clc; clear all;
for n=50:50:1000
h=(2-1)/n; x=1:h:2; f=@(x)exp(x); yamuk=0;simpson=0;tam=exp(2)-exp(1);
for i=1:n
     yamuk=yamuk+h*(f(x(i))+f(x(i+1)))/2;
 end
for i=1:2:n-1
    simpson=simpson+h*(f(x(i))+4*f(x(i+1))+f(x(i+2)))/3;
end
fprintf('n=%g icin yamuk=%.10f,yamuk hata=%.10f,simpson=%.10f,simpson hata=%.10f\n'....
    ,n,yamuk,abs(yamuk-tam),simpson,abs(simpson-tam));
end



















