% ALGO1 FROM FATKULLIN VANDEN EIJNDEN -small step != big step || euler
% steps
% variables
clear
close all

% per full model
% N = 1000;
% dt = 1/1000;

%% Lorenz 96 parameters
F = 10;                  % forcing parameter
hx = -0.8;
hy = 1;
epsi = 2^(-7)

exp = 7;
Dt = 2^(-exp);            % timestep for SLOW var
%dt = 0;
dt = 2^(-14);           % timestep for FAST var

t_final = 2^(10);          % duration of the simulation in natural time units
NT = fix(t_final/Dt);   % no. of discrete time steps for SLOW var
nosc = 9;               % no. of oscillators [K]
fps = 8;                % no. of fast variables per slow variable [J]
nosc_fast = fps*nosc;   % no. of fast variables
%N0 = 1e0;
N0 = 1e6;
N1 = 4;
R = 1;
% Initial conditions
x = zeros([nosc NT]);
x(1:nosc,1) = 2*rand([nosc 1])-1;
z = zeros([nosc_fast NT R]);

F_m = zeros(nosc,1);
g = @(y, x, j, N, fps, hy) 1/epsi*( y(mod(j,N)+1)*(y(mod(j-2,N)+1) - y(mod(j+1,N)+1)) - y(j) + hy*x(1 + floor((j-1)/fps)));

(R*N1*dt)/Dt
NT
%% Euler integration
tic
k=1;
    for m = 2:NT
        
        if(mod(m,floor(NT/100))==0)
            disp("----------" + k + "% ----------");
            %disp("step " + m + " of " + NT);
            k=k+1;
        end
        F_m = zeros(nosc,1);
        % fast variables
        for r = 1:R
            
            % sample the fast variables (modify to assign the input if available)
            if (m==2)
                z(:,1, r) = rand([nosc_fast 1]) - 1/2;
            end
            
            % scarto i primi N1 termini
            clearvars zz
            zz(:, 1) = z(:, m-1, r);

            if (m==2)
                for n = 2:N0+1
                    if(mod(n,floor(N0/10))==0)
                        disp("step " + n + " of " + N0);
                    end
                    id_x = 1;
                    for j = 1:nosc_fast
                        id_x = 1 + floor((j-1)/fps);
                        zz(j,n) = zz(j,n-1) + dt/epsi*( zz(mod(j,nosc_fast)+1,n-1)*(zz(mod(j-2,nosc_fast)+1,n-1) - zz(mod(j+1,nosc_fast)+1,n-1)) - zz(j,n-1) + hy*x(id_x,m-1) ) ;
                        % one step euler
                        % zz(j,n) = zz(j,n-1) + dt*g(zz(:,n-1), x(:, m-1), j, nosc_fast, fps, hy);
                        % one step RK4
                        %zz(j,n) = rungekutta4(@(t, z) g([zz(1:j-1,n-1); z; zz(j+1:end,n-1)], x(:,m-1), j, nosc_fast, fps, hy), (n-2)*dt, dt, zz(j,n-1), 1);
                    end %j
                end
            else
                for n = 2:N1+1
                    
                    id_x = 1;
                    for j = 1:nosc_fast
                        id_x = 1 + floor((j-1)/fps);
                        zz(j,n) = zz(j,n-1) + dt/epsi*( zz(mod(j,nosc_fast)+1,n-1)*(zz(mod(j-2,nosc_fast)+1,n-1) - zz(mod(j+1,nosc_fast)+1,n-1)) - zz(j,n-1) + hy*x(id_x,m-1) ) ;
                        % one step euler
                        % zz(j,n) = zz(j,n-1) + dt*g(zz(:,n-1), x(:, m-1), j, nosc_fast, fps, hy);
                        % one step RK4
                        %zz(j,n) = rungekutta4(@(t, z) g([zz(1:j-1,n-1); z; zz(j+1:end,n-1)], x(:,m-1), j, nosc_fast, fps, hy), (n-2)*dt, dt, zz(j,n-1), 1);
                    end %j
                end
            end
            z(:, m, r) = zz(:, end);

           for i=1:nosc
                F_m(i) = F_m(i) + 1/R*( -x(mod(i-2, nosc)+1,m-1)*(x(mod(i-3, nosc)+1,m-1) - x(mod(i,nosc)+1,m-1)) - x(i,m-1) + F + hx/fps*sum(z((i-1)*fps+(1:fps), m, r)) ) ;
           end
        end

        % slow variables
        for i=1:nosc
            % one step euler
            x(i,m) = x(i,m-1) + Dt*F_m(i);
            % on step RK4
        end %i
        
        
    end %m
    
    ok = isempty( find( isnan(x(1,:)) | isinf(x(1,:)), 1 ) );   %Check wheter the integration have worked or not

    if not(ok) 
        disp('INTEGRATION DID NOT WORK') 
    end
toc

%% save data
writematrix(x, "data/x_comparison/x_fat_" + exp + "_R" + R + "_N" + N1 + ".csv");

%% FIGURES

%% pdf (mean of all of them)
clear p

for k = 1:nosc
    [p(k,:), a] = histcounts(x(k,:), 50, 'Normalization','pdf', 'BinLimits', [-10 15]);

end

meanpdf = mean(p);

figure()
plot(a(1:end-1), meanpdf)
legend('Fatkullin multiscale method')
title('pdf')
subtitle("\Delta t = 2^{-" + exp + "}, R = " + R + ", N1 = " + N1)
xlabel('X_k')
ylabel('probability density')

%% error
norm(meanpdf'- direct_pdf)/norm(direct_pdf)

%% kolmogorov smirnov : returns 0 if the null hp is not rejected
x_direct = readmatrix('direct_x.csv');

for i=1:nosc
    x_direct_conc((i-1)*size(x_direct,2) + 1 : i*size(x_direct,2)) = x_direct(i, :)';
    x_conc((i-1)*NT + 1 : i*NT) = x(i, :);
    end

[h,p, kstat] = kstest2(x_direct_conc, x_conc)

%% 

%% auto correlation function
 clear acf acf_z lags1 meanacf meanacf_z
x_direct = readmatrix('direct_x.csv');

% manual way
% nosc = 9;
% for k = 1:nosc-1
%     k
%     acf(k,:) = calculate_correlation(x(k,:), x(k,:));
%     [acf_2(k,:),lags1] = autocorr(x(k,:), size(x,2)-1);
% end
% 
% 
% meanacf_manual = mean(acf_2);
% meanacf = mean(acf);
%%
% dt = 2^(-9);
% plot(lags1*dt, meanacf), grid on; hold on;
% plot(lags1*dt, meanacf_manual);
% legend('buit-in', 'manual')
% xlim([0, 30])

%% built-in function
for k = 1:5
    [acf(k,:), lags1] = autocorr(x(k,:), size(x,2)-1);
    [acf_direct(k,:), lags_direct] = autocorr(x_direct(k,:), size(x_direct,2)-1);
end
meanacf = mean(acf);
meanacf_R = mean(acf_R);
direct_acf = mean(acf_direct);

% fast var
% for k = 1:nosc_fast-1
%     [acf_z(k,:)] = autocorr(z(k,:), size(x,2)-1);
% end
% meanacf_z = mean(acf_z);

%%
figure()
plot(lags1*Dt, meanacf), grid on; hold on;
%plot(lags1*Dt, meanacf_z), grid on; hold on;

plot(lags_direct*dt, direct_acf)
xlim([0,30])
title('normalized auto correlation function')
subtitle("\Delta t = 2^{-" + exp + "}, R = " + R + ", N1 = " + N1)
ylim([-0.5 1.1])
legend('multiscale method', 'direct solver Au')
legend('X', 'X direct solver')
%%

% %% cross correlation function
% clear ccf ccf_z
% for k = 1:nosc-1
%     [ccf(k,:),lags2] = crosscorr(x(k+1,:), x(k,:), size(x,2)-1);
% end
% meanccf = mean(ccf);
% 
% %fast var
% for k = 1:nosc_fast-1
%     [ccf_z(k,:)] = crosscorr(z(k+1,:,1), z(k,:,1), size(z,2)-1);
% end
% meanccf_z = mean(ccf_z);
% 
% figure()
% plot(lags2*dt, meanccf), grid on; hold on;
% plot(lags2*dt, meanccf_z)
% xlim([0,30])
% ylim([-0.5 1.1])
% title('normalized cross correlation function')
% 

