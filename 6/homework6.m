clc
clear
close
rand('seed',666);

%% Основные параметры
z0 = 0; % Начальная точка траектории
zN = 1; % Конечная точка траектории
u_x_z = @(x,z) 5*x*(exp(-4/5*z) - exp(-z)) + 4*exp(-z); % Аналитическое решение
u_x = @(x) u_x_z(x,zN); % Аналитическое решение в точке zN
x_step = 0.05; % Шаг сетки
x_grid = 0:x_step:1; % Сетка по координате x
res = u_x(x_grid+x_step/2); % Ответ
f = 4; % Неоднородность в уравнении, вообще говоря, f = f(x);
p_x = 1; % Плотность вероятности по x. Общая плотность p(x,z) = A*f(x)*delta(z)
sigma_t = 1; % Полное сечение рассеяния (коэффициент при u(x) в уравнении)
g = @(z) exp(-sigma_t*(1-z));
k = 4; % Отношение плотности источников s(x,y) = f(x)*delta(z) к общей плотности p(x,y) = A*f(x)*delta(z), A = 4 из условия номировки на 1
%K = @(x,y) x*z^3; % Ядро интеграла
sigma_s = @(x,z) x*z^3; % Сечение рассеяния
omega = 1; % Направление движение. Формально единичный трехмерный вектор, но в нашем случае это (0, 0, 1)
N = 100000; % Число экспериментов
w0 = 1; % Начальный вес
U_x = zeros(1, length(x_grid)); % Результаты
U_x_n = zeros(N, length(x_grid)); % Промежуточные результаты
len = 0; % Длина пробега

%% Главный цикл
for s = 1:N
    x0 = rand();
    U_x_n(s, theta(x0, x_grid)) = U_x_n(s, theta(x0, x_grid)) + 4*w0*g(z0);
    len = -log(rand());
    z = z0 + len*omega;
    i = 1;
    x_pr = x0;
    while(z < zN)
        x = sqrt(rand());
        U_x_n(s, theta(x, x_grid)) = U_x_n(s, theta(x, x_grid)) + 4*weight(i,x_pr)*g(z);
        len = -log(rand());
        i = i + 1;
        z = z + len*omega;
        x_pr = x;
    end
end

%% Результат
for s = 1:N
    for k = 1:length(x_grid)
        U_x(k) = U_x(k) + U_x_n(s, k)/(N*x_step);
    end
end
for k = 1:length(x_grid)-1
    fprintf('%f %f %f\n',x_grid(k)+x_step/2,res(k),U_x(k));
end
