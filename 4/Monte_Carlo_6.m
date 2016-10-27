clc
clear
close
%% Задание основных параметров
rand('seed',666);
N = 50000;
N_start  = 1000;
N_step = 1000;
C = 16.1975;
Ni = N_start:N_step:N;
res_old = zeros(length(Ni),1);
res_new = zeros(length(Ni),1);

func = @(x,y,z,u) exp(x*y*z*u)/sqrt(x*y*z*u);
f = @(x,y,z,u) 16*exp(x*y*z*u);
p = @(x,y,z,u) 1/(16*sqrt(x*y*z*u));

I = 0;
for k=0:1000
   I = I + 16*1/(factorial(k)*(2*k+1)^4);
end

y = rand(N,4);

%% Вычисление последовательности Холтона
p2 = zeros(1, N);
p3 = zeros(1, N);
p5 = zeros(1, N);
for k = 1:N
    M = k + 2;
    a2 = zeros(1, M);
    a3 = zeros(1, M);
    a5 = zeros(1, M);
        
    a2(1) = k;
    a3(1) = k;    
    a5(1) = k;
    
    for j = 2:M
        a2(j) = fix(a2(j - 1)/2);
        a3(j) = fix(a3(j - 1)/3);
        a5(j) = fix(a5(j - 1)/5);
        a2(j - 1) = mod(a2(j - 1), 2);
        a3(j - 1) = mod(a3(j - 1), 3);    
        a5(j - 1) = mod(a5(j - 1), 5);
        p2(k) = p2(k) + a2(j - 1)*2^(-j + 1);
        p3(k) = p3(k) + a3(j - 1)*3^(-j + 1);
        p5(k) = p5(k) + a5(j - 1)*5^(-j + 1);
    end
    if(mod(k,1000) == 0)
        k
    end
end

%% Основной цикл
for k = 1:length(Ni)
    
    % Старый алгоритм
    M_old = 0;
    for i=1:Ni(k)
        M_old = M_old + C + 16*(exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ...
                - 1 - y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ;
    end
    M_old = M_old / Ni(k);
    res_old(k) = log(abs(M_old - I))/log(10);
    
    % Новый алгоритм
    M_new = 0;
    for i=1:Ni(k)
        M_new = M_new + C + 16*(exp((i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2) ...
                - 1 - (i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2) ;
    end
    M_new = M_new / Ni(k);
    res_new(k) = log(abs(M_new - I))/log(10);
    
end

%% Построени графиков
% МНК
e1 = polyfit(log(Ni)/log(10), res_old',1);
e2 = polyfit(log(Ni)/log(10), res_new',1);

grid on;
hold on
    plot(log(Ni)/log(10), res_old, 'DisplayName', 'Псевдослучайные числа');
    plot(log(Ni)/log(10), res_new, 'DisplayName', 'Последовательность Холтона');
    plot(log(Ni)/log(10), e1(2)+e1(1)*log(Ni)/log(10), 'DisplayName', 'Псевдослучайные числа МНК');
    plot(log(Ni)/log(10), e2(2)+e2(1)*log(Ni)/log(10), 'DisplayName', 'Последовательность Холтона МНК');
hold off
xlabel('lg N');
ylabel('lg |teta - I|');
legend('show');
title('Сходимость алгоритма с выделением главной части');
