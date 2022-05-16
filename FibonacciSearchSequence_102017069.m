clear all
clc
% defining the objective function

f=@(x) x.*(x-2);
l=0;
r=1.5;
moe=0.25;
n=1/moe;
lo=r-l;

% plotting the function
t = linspace(l,r,100);
plot(t,f(t),'k','LineWidth',2);

% compute the fibonacci series
fib = ones(1,n); %set all 1 initially
for i=3:n+1
    fib(i)=fib(i-1)+fib(i-2);
end

for k=1:n
    ratio=(fib(n+1-k)./fib(n+2-k));
    x2=l+ratio.*lo;
    x1=l+r-x2;
    fx1=f(x1);
    fx2=f(x2);
    rsl(k,:)=[l r x1 x2 fx1 fx2]; % for printing purpose
    
    if fx1<fx2
        r=x2;
    elseif fx1>fx2
        l=x1;
    elseif fx1==fx2
        if min(abs(x1),abs(l))==abs(l)
            r=x2;
        else
            l=x1;
        end
    end
end

variables={'l','r','x1','x2','fx1','fx2'};
resl=array2table(rsl);
resl.Properties.VariableNames(1:6)=variables;
disp(resl)
xopt=(l+r)/2
fopt=f(xopt)
fprintf('Optimal value of x = %f \n',xopt);
fprintf('Optimal value of f(x) = %f \n',fopt);
