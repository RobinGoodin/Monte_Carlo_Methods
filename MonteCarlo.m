clc
clear
close
rand('seed',666);

%%
x0 = 1;
K = @(x,y) x*y^3;
F = 4;
N = 100;
n = 1000;
W = zeros(N,1);
W(1)=1;
g = rand(N,1);
g(1) = x0;
for i=2:N
    W(i)=W(i-1)*K(g(i-1), g(i));
end
theta = zeros(N,1);
for i=1:N
    for j=1:n
        theta(i) = theta(i) + W(i)*F;
    end
    theta(i) = theta(i)/n;
end
res = 0;
for i=1:N
    res = res + theta(i);
end
res
