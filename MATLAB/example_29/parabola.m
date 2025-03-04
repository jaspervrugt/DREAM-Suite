function c = parabola(x,y)

A = [x(1)^2 x(1) 1;x(2)^2 x(2) 1;x(3)^2 x(3) 1];
x=linspace(-10,10,100);
b=y;
c=A/b;
f=c(1,1)*x.^2+c(2,1)*x+c(3,1);
plot(x,f)
xlabel('x')
ylabel('y')

% A = [   x(1)^2 x(1) 1; 
%         x(2)^2 x(2) 1; 
%         x(3)^2 x(3) 1];
% b = y(:);
% c = A\b;
% b = y;
% c2 = A/b
% plot(x,y,'ro'); hold on; 
% x = linspace(-10,10,100);
% y = c(1)*x.^2+c(2)*x+c(3);
% y2 = c2(1)*x.^2+c2(2)*x+c2(3);
% plot(x,y,'k');
% plot(x,y2,'b');

end