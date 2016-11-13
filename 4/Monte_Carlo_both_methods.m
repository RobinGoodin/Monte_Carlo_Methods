clc
clear
close
%% Задание основных параметров
rand('seed',666);
N = 50000;
N_start  = 1000;
N_step = 1000;
C = 16.1975308642;
Ni = N_start:N_step:N;
res_old1 = zeros(length(Ni),1);
res_new1 = zeros(length(Ni),1);
res_old2 = zeros(length(Ni),1);
res_new2 = zeros(length(Ni),1);

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
    M_old1 = 0;
    for i=1:Ni(k)
        M_old1 = M_old1 + 16*exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2);
    end
    M_old1 = M_old1 / Ni(k);
    res_old1(k) = log(abs(M_old1 - I))/log(10);
    
    % Новый алгоритм
    M_new1 = 0;
    for i=1:Ni(k)
        M_new1 = M_new1 + 16*exp((i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2);
    end
    M_new1 = M_new1 / Ni(k);
    res_new1(k) = log(abs(M_new1 - I))/log(10);
    
    % Старый алгоритм
    M_old2 = 0;
    for i=1:Ni(k)
        M_old2 = M_old2 + C + 16*(exp(y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ...
                - 1 - y(i,1)^2*y(i,2)^2*y(i,3)^2*y(i,4)^2) ;
    end
    M_old2 = M_old2 / Ni(k);
    res_old2(k) = log(abs(M_old2 - I))/log(10);
    
    % Новый алгоритм
    M_new2 = 0;
    for i=1:Ni(k)
        M_new2 = M_new2 + C + 16*(exp((i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2) ...
                - 1 - (i/Ni(k))^2*p2(i)^2*p3(i)^2*p5(i)^2) ;
    end
    M_new2 = M_new2 / Ni(k);
    res_new2(k) = log(abs(M_new2 - I))/log(10);
    
end

%% Построени графиков
% МНК
e11 = polyfit(log(Ni)/log(10), res_old1',1);
e12 = polyfit(log(Ni)/log(10), res_new1',1);
e21 = polyfit(log(Ni)/log(10), res_old2',1);
e22 = polyfit(log(Ni)/log(10), res_new2',1);

grid on;
hold on
    plot(log(Ni)/log(10), res_old1, 'DisplayName', 'Стандартный алгоритм (ПСЧ)');
    plot(log(Ni)/log(10), res_new1, 'DisplayName', 'Стандартный алгоритм (Холтон)');
    plot(log(Ni)/log(10), e11(2)+e11(1)*log(Ni)/log(10), 'DisplayName', 'Стандартный алгоритм (ПСЧ) МНК');
    plot(log(Ni)/log(10), e12(2)+e12(1)*log(Ni)/log(10), 'DisplayName', 'Стандартный алгоритм (Холтон) МНК');

    plot(log(Ni)/log(10), res_old2, 'DisplayName', 'Выделение главной части (ПСЧ)');
    plot(log(Ni)/log(10), res_new2, 'DisplayName', 'Выделение главной части (Холтон)');
    plot(log(Ni)/log(10), e21(2)+e21(1)*log(Ni)/log(10), 'DisplayName', 'Выделение главной части (ПСЧ) МНК');
    plot(log(Ni)/log(10), e22(2)+e22(1)*log(Ni)/log(10), 'DisplayName', 'Выделение главной части (Холтон) МНК');
hold off
xlabel('lg N');
ylabel('lg |teta - I|');
legend('show');
title('Сходимость алгоритма с выделением главной части');
