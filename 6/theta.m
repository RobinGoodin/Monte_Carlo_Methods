function k = theta(x, x_grid)
%% Описание
% Функция возвращает номер ячейки x_grid, в которую попал x 

%% Реализация
N = length(x_grid);

if(x < x_grid(1) || x > x_grid(N))
    error('Недопустимое значение x');
end
for i = 2:N
    if(x < x_grid(i))
        k = i-1;
        return;
    end
end