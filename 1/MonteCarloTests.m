%Система тестов для проверки генератора псевдослучайных чисел
clc;
clear;
close all;
N = 600;
seed = 10;
rand('seed', seed);
y = rand([1 N]);
y = sort(y); %Последовательность упорядоченных псевдослучайных чисел

display('Тест 1');
%Проверка выборочной функции распределения по критерию Колмогорова
d = zeros(N,2);
for k=1:N
    d(k,1) = abs(y(k)-(k-1)/N);
    d(k,2) = abs(k/N-y(k));
end
D = max(max(d));
Kappa = D*sqrt(N)
if(Kappa < 1.36)
    display('Тест пройден');
else
    display('Тест НЕ пройден');
end;

display('Тест 2');
%Проверка расположения одномерных точек
CHI = zeros(10,1);
for s = 1:10
    N = 600*2^s;
    rand('seed', seed);
    y = rand([1 N]);
    %y = sort(y);
    v = zeros(16,1);
    v1 = 0;
    v2 = 1/16;
    next = 1/16;
    for j=1:16
        for i=1:N
            if((y(i) > v1) && (y(i) <= v2))
                v(j)=v(j)+1;
            end
        end
        v1=v1+next;
        v2=v2+next;
    end
    chi = 0;
    for j=1:16
       chi = chi + (v(j)-N/16)^2;
    end
    CHI(s) = 16/N*chi;
end;
chi1 = max(CHI)
%r = 16;
%K = @(x) chi2pdf(x,r-1);
%p = int(K,chi,Inf)*100
if(chi1 < 25)%22.3
    display('Тест пройден');
else
    display('Тест НЕ пройден');
end;
    
display('Тест 3');
%Проверка расположения двумерных точек
CHI = zeros(10,1);
for s = 2:10
    N = 600*2^s;
    rand('seed', seed);
    y = rand([1 N]);
    %y = sort(y);
    v = zeros(8,8);
    next = 1/8;
    x1 = 0;
    x2 = 1/8;
    for i=1:8
        y1 = 0;
        y2 = 1/8;
        for j=1:8
            for k=1:2:N
                if((y(k) > x1) && (y(k) <= x2) ...
                  && (y(k+1) > y1) && (y(k+1) <= y2))
                    v(i,j)=v(i,j)+1;
                end
            end
            y1=y1+next;
            y2=y2+next;
        end
        x1=x1+next;
        x2=x2+next;
    end
    chi = 0;
    for i=1:8
        for j=1:8
            chi = chi + (v(i,j)-N/(2*64))^2;
        end
    end
    CHI(s) = 64*2/N*chi;
end;
chi2 = max(CHI)
%r = 64;
%K = @(x) chi2pdf(x,r-1);
%p = integral(K,chi,Inf)*100
if(chi2 < 82.5)%77.7
    display('Тест пройден');
else
    display('Тест НЕ пройден');
end;

display('Тест 4');
%Проверка расположения трехмерных точек
CHI = zeros(10,1);
for s = 4:10
    N = 600*2^s;
    rand('seed', seed);
    y = rand([1 N]);
    %y = sort(y);
    v = zeros(5,5,5);
    next = 1/5;
    x1 = 0;
    x2 = 1/5;
    for i=1:5
        y1 = 0;
        y2 = 1/5;
        for j=1:5
            z1 = 0;
            z2 = 1/5;
            for k=1:5
                for t=1:3:N
                    if((y(t) > x1) && (y(t) <= x2) ...
                      && (y(t+1) > y1) && (y(t+1) <= y2) ...
                      && (y(t+2) > z1) && (y(t+2) <= z2))
                        v(i,j,k)=v(i,j,k)+1;
                    end
                end
                z1=z1+next;
                z2=z2+next;
            end
            y1=y1+next;
            y2=y2+next;
        end
        x1=x1+next;
        x2=x2+next;
    end
    chi = 0;
    for i=1:5
        for j=1:5
            for k=1:5
                chi = chi + (v(i,j,k)-N/(3*125))^2;
            end
        end
    end
    CHI(s) = 125*3/N*chi;
end;
chi3 = max(CHI)
%r = 125;
%K = @(x) chi2pdf(x,r-1);
%p = integral(K,chi,Inf)*100
if(chi3 < 151)%145
    display('Тест пройден');
else
    display('Тест НЕ пройден');
end;

display('Тест 5');
%Проверка расположения четырехмерных точек
CHI = zeros(10,1);
for s = 8:10
    N = 600*2^s;
    rand('seed', seed);
    y = rand([1 N]);
    %y = sort(y);
    v = zeros(4,4,4,4);
    next = 1/4;
    x1 = 0;
    x2 = 1/4;
    for i=1:4
        y1 = 0;
        y2 = 1/4;
        for j=1:4
            z1 = 0;
            z2 = 1/4;
            for k=1:4
                t1 = 0;
                t2 = 1/4;
                for m=1:4
                    for t=1:4:N
                        if((y(t) > x1) && (y(t) <= x2) ...
                          && (y(t+1) > y1) && (y(t+1) <= y2) ...
                          && (y(t+2) > z1) && (y(t+2) <= z2) ...
                          && (y(t+3) > t1) && (y(t+3) <= t2))
                            v(i,j,k,m)=v(i,j,k,m)+1;
                        end
                    end
                    t1=t1+next;
                    t2=t2+next;
                end
                z1=z1+next;
                z2=z2+next;
            end
            y1=y1+next;
            y2=y2+next;
        end
        x1=x1+next;
        x2=x2+next;
    end
    chi = 0;
    for i=1:4
        for j=1:4
            for k=1:4
                for m=1:4
                    chi = chi + (v(i,j,k,m)-N/(4*256))^2;
                end
            end
        end
    end
    CHI(s) = 256*4/N*chi;
end;
chi4 = max(CHI)
if(chi4 < 293)%284
    display('Тест пройден');
else
    display('Тест НЕ пройден');
end;

display('Тест 6');
%Проверка расположения четырехмерных точек
CHI = zeros(10,1);
for s = 8:10
    N = 600*2^s;
    rand('seed', seed);
    y = rand([1 N]);
    n1 = zeros(5,1);
    buf = floor(y(1)*10);
    m = 0;
    for i=2:N
        if(floor(y(i)*10) == buf)
            m = m + 1;
        elseif(m <= 4)
            n1(m + 1) = n1(m + 1) + 1;
            if(i <= N-1)
                buf = floor(y(i)*10);
            end
            m = 0;
        elseif(m > 5)
            n1(5) = n1(5) + 1;
            if(i <= N-1)
                buf = floor(y(i)*10);
            end
            m = 0;
        end
    end
    e = 0;
    n = sum(n1);
    pj = [0.9,0.09,0.009,0.0009,10^(-4)];
    for j=1:4
        e = e + (n1(j) - n*pj(j))^2/(n*pj(j));
    end
    e = e + (n1(5) - n*pj(5))^2/(n*pj(5));
    CHI(s) = e;
end
chi5 = max(CHI)
if(chi5 < 9.49)%7.78
    display('Тест пройден');
else
    display('Тест НЕ пройден');
end;
%{
Kappa
chi1
chi2
chi3
chi4
chi5
%}
