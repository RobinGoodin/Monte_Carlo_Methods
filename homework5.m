clc
clear
close
rand('seed',666);
% G - ������ �����������, ������ ��� ��� ����������� ������ ����
% ���������� ������ ������ ��� ��������� ������.
%% 5.25
x0 = 1;
K = @(x,y) x*y^3;
F = @(x) 4;
p = @(x,y) 4*y^3;
S = @(x) x/4;

%% �������� ���������
N = 1000000;
W = 1;
U = 0;
%% ������� ����
theta = zeros(N,1);
for s = 1:N
    flag = 1;
    ksi1 = x0;
    theta(s) = W * F(ksi1);
    while (flag == 1)
        gamma = rand();
        if (gamma > S(ksi1))
            flag = 0;
        else
            ksi1 = (4*gamma/ksi1)^(1/4);
            theta(s) = theta(s) + W * F(ksi1);
        end
    end
end
%% ���������
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

