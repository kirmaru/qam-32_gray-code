clc; clear all; close all;
 
%% параметры
q     = 32;            
k     = log2(q);        
N     = 1000;            
Vmod  = 300;             
T     = 1/Vmod;          
dt    = T/N;             
t     = 0:dt:T-dt;        
 
%% базисные функции
f0   = 9000;
phi1 = sqrt(2/T)*cos(2*pi*f0*t);
phi2 = sqrt(2/T)*sin(2*pi*f0*t);
 
%% генерация созвездия и шахматная маска
levels     = -7:2:7;
[I_all,Q_all] = meshgrid(levels, levels);
%% созвездие до прореживания
figure;
scatter(I_all(:), Q_all(:), 60, 'b', 'filled');
xlabel('I'); ylabel('Q');
title('Созвездие q=64');
axis equal; grid on;
 
%%
mask_chess   = mod((0:7)' + (0:7), 2)==0;
I_all = I_all(:);  Q_all = Q_all(:);
base_points = [I_all(mask_chess), Q_all(mask_chess)];  % 32 точки
 
%% созвездие после прореживания
figure;
scatter(base_points(:,1), base_points(:,2), 60, 'b', 'filled');
xlabel('I'); ylabel('Q');
title('Созвездие q=32 после прореживания');
axis equal; grid on;
 
%% поворот на 45 градусов
theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
rotated_points = (R * base_points.').';   
 
%% перестановка по коду Грея
bin2gray = zeros(q,1);
for b = 0:q-1
    bin2gray(b+1) = bitxor(b, bitshift(b,-1));
end
signal_points = rotated_points(bin2gray+1, :);
 
%% средняя энергия
Es = mean(sum(signal_points.^2, 2));
 
%% Формирование сигналов
signals = zeros(q, length(t));
for idx = 1:q
    signals(idx,:) = signal_points(idx,1)*phi1 + signal_points(idx,2)*phi2;
end
 
%% созвездие после поворота
figure; hold on;
scatter(signal_points(:,1), signal_points(:,2), 60, 'b', 'filled');
xlabel('I'); ylabel('Q');
title('Созвездие q=32 после прореживания и поворота');
axis equal; grid on;
hold off;
 
%% созвездие с решающими областями
margin = 3;
x_min = min(signal_points(:,1)) - margin;
x_max = max(signal_points(:,1)) + margin;
y_min = min(signal_points(:,2)) - margin;
y_max = max(signal_points(:,2)) + margin;
 
[xq, yq] = meshgrid(linspace(x_min, x_max, 1000), linspace(y_min, y_max, 1000));
D = zeros([size(xq) q]);
for idx = 1:q
    D(:,:,idx) = (xq - signal_points(idx,1)).^2 + (yq - signal_points(idx,2)).^2;
end
[~, regions] = min(D, [], 3);
 
figure; hold on;
contour(xq, yq, regions, q-1, 'k', 'LineWidth', 1.5);
scatter(signal_points(:,1), signal_points(:,2), 60, 'b', 'filled');
xlabel('I'); ylabel('Q');
title('Решающие области');
axis equal; grid on;
xlim([x_min x_max]);
ylim([y_min y_max]);
hold off;
 
%% построение LUT таблицы
 
projI = unique(signal_points(:,1))';
projQ = unique(signal_points(:,2))';
Ni = numel(projI);
Nq = numel(projQ);
 
bI = [-inf, (projI(1:end-1)+projI(2:end))/2, +inf];
bQ = [-inf, (projQ(1:end-1)+projQ(2:end))/2, +inf]; 
 
lut = zeros(Nq, Ni);
for qi = 1:Nq
    for ii = 1:Ni
        pt = [projI(ii), projQ(qi)];
        [~, lut(qi,ii)] = min( sum((signal_points - pt).^2, 2) );
    end
end
 
%% моделирование приема
EbN0_dB   = 0:2:18;
EbN0      = 10.^(EbN0_dB/10);
EsN0      = EbN0 * k;
theoryBER = (1/k)*(1 - (1 - 2*qfunc(sqrt(3*k*EbN0/(q-1)))).^2);
theorySER = 1 - (1 - 2*qfunc(sqrt(3*EsN0/(q-1)))).^2;
bit_errors = zeros(size(EbN0));
sym_errors = zeros(size(EbN0));
 
fprintf('Теоретические значения:\n');
for ei = 1:length(EbN0_dB)
    fprintf('Eb/N0 = %2d dB: theory BER = %.3e, theory SER = %.3e\n', ...
        EbN0_dB(ei), theoryBER(ei), theorySER(ei));
end
 
parfor ei = 1:length(EbN0)
    N0    = Es/(k * EbN0(ei));       
    sigma = sqrt(N0/(2*dt));
    errs      = 0; total_bits = 0;
    sym_errs  = 0; total_sym  = 0;
    target_err = 1e4;
    max_bits   = 10e6 * k;
 
    while errs < target_err && total_bits < max_bits
        bits = randi([0 1],1,k);
        db   = bi2de(bits,'left-msb');
        gz   = bitxor(db, bitshift(db,-1));
        idx  = gz + 1;
        s    = signals(idx, :);
 
        r = s + sigma*randn(size(s));
 
        r1 = sum(r .* phi1)*dt;
        r2 = sum(r .* phi2)*dt;
 
        i_idx   = discretize(r1, bI);   
        q_idx   = discretize(r2, bQ);  
        idx_hat = lut(q_idx, i_idx);    
 
        g_rx    = idx_hat - 1;
        b_rx    = g_rx;
        for sh = 1:k-1
            b_rx = bitxor(b_rx, bitshift(g_rx, -sh));
        end
        bits_rx = de2bi(b_rx, k, 'left-msb');
 
        errs      = errs + sum(bits ~= bits_rx);
        total_bits= total_bits + k;
        sym_errs  = sym_errs + double(g_rx ~= gz);
        total_sym = total_sym + 1;
        
    end
 
    bit_errors(ei) = errs / total_bits;
    sym_errors(ei) = sym_errs / total_sym;
    fprintf('Eb/N0 = %2d dB: BER = %.2e, SER = %.2e\n', ...
            EbN0_dB(ei), bit_errors(ei), sym_errors(ei));
 
end
 
%% графики вероятностей ошибок
figure;
semilogy(EbN0_dB, bit_errors, 'ro-', EbN0_dB, theoryBER, 'b*-');
xlabel('Eb/N0, dB'); ylabel('BER');
legend('симуляция','теория'); grid on;
title('Вероятность ошибки на бит');
 
figure;
semilogy(EbN0_dB, sym_errors, 'ks--', EbN0_dB, theorySER, 'g^-');
xlabel('Eb/N0, dB'); ylabel('SER');
legend('симуляция','теория'); grid on;
title('Вероятность ошибки на символ');
