function k = theta(x, x_grid)
%% ��������
% ������� ���������� ����� ������ x_grid, � ������� ����� x 

%% ����������
N = length(x_grid);

if(x < x_grid(1) || x > x_grid(N))
    error('������������ �������� x');
end
for i = 2:N
    if(x < x_grid(i))
        k = i-1;
        return;
    end
end