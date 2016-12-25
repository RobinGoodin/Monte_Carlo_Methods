clc
clear
close
rand('seed',666);

%% �������� ���������
z0 = 0; % ��������� ����� ����������
zN = 1; % �������� ����� ����������
u_x_z = @(x,z) 5*x*(exp(-4/5*z) - exp(-z)) + 4*exp(-z); % ������������� �������
u_x = @(x) u_x_z(x,zN); % ������������� ������� � ����� zN
x_step = 0.05; % ��� �����
x_grid = 0:x_step:1; % ����� �� ���������� x
res = u_x(x_grid+x_step/2); % �����
f = 4; % �������������� � ���������, ������ ������, f = f(x);
p_x = 1; % ��������� ����������� �� x. ����� ��������� p(x,z) = A*f(x)*delta(z)
sigma_t = 1; % ������ ������� ��������� (����������� ��� u(x) � ���������)
g = @(z) exp(-sigma_t*(1-z));
k = 4; % ��������� ��������� ���������� s(x,y) = f(x)*delta(z) � ����� ��������� p(x,y) = A*f(x)*delta(z), A = 4 �� ������� ��������� �� 1
%K = @(x,y) x*z^3; % ���� ���������
sigma_s = @(x,z) x*z^3; % ������� ���������
omega = 1; % ����������� ��������. ��������� ��������� ���������� ������, �� � ����� ������ ��� (0, 0, 1)
N = 100000; % ����� �������������
w0 = 1; % ��������� ���
U_x = zeros(1, length(x_grid)); % ����������
U_x_n = zeros(N, length(x_grid)); % ������������� ����������
len = 0; % ����� �������

%% ������� ����
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

%% ���������
for s = 1:N
    for k = 1:length(x_grid)
        U_x(k) = U_x(k) + U_x_n(s, k)/(N*x_step);
    end
end
for k = 1:length(x_grid)-1
    fprintf('%f %f %f\n',x_grid(k)+x_step/2,res(k),U_x(k));
end
