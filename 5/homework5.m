clc
clear
close
rand('seed',666);
% G - ошибка погрешности, вообще тут еще коэффициент должен быть
% посмотреть первые лекции про вероятную ошибку.
%% 5.25
x0 = 1;
K = @(x,y) x*y^3;
F = @(x) 4;
p = @(x,y) 4*y^3;
S = @(x) x/4;

%% Основные параметры
N = 100000;
W = 1;
U = 0;
it = 1;
max = it;
%% Главный цикл
theta = zeros(N,1);
for s = 1:N
    flag = 1;
    ksi1 = x0;
    theta(s) = W * F(ksi1);
    while (flag == 1)
        
        gamma = rand();
        if (gamma > S(ksi1))
            flag = 0;
            if(it > max)
                max = it;
            end
            it = 1;
        else
            ksi1 = (4*gamma/ksi1)^(1/4);
            theta(s) = theta(s) + W * F(ksi1);
            it = it + 1;
        end
    end
end
%% Результат
for s = 1:N
    U = U + theta(s);
end
U = U / N

T = 0;
for s = 1:N
    T = T + theta(s)*theta(s);
end
T = T / N;

D = (T - U*U)/(N-1)

G = sqrt(D)

max