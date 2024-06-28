% lorenz 96 model multiscale approximation via averaging
% Very crude method 
%% creating the denstity function of the fast variable for the reduced models
%% by sampling the full model 
% parameter values for initial full model
epsilon = 1/128;
K = 9;
J = 8;
F = 10;
hx = -0.8;
hy = 1;
dt = 2^(-14);
T = 40;
TIME = floor(T/dt);
tt = 1;
MeanY = [K,TIME];
% initialize x and y
X(:,tt) = -1 + (2)*rand(K,1);
Y(:,tt) = -1 + (2)*rand(K*J,1);
Y(:,tt) = -epsilon + (2*epsilon)*rand(K*J,1);
for tt = 2:TIME
    for ii = 1:K 
        summ = 0;
        for jj = 1:J  % compute all Y
            M = ((ii-1)*J)+jj;  % subscript of Y, runs from 1 to J*K 
            if     ii == 1 && jj == 1 
                Y(M,tt) = Y(M,tt-1) + (1/epsilon)*dt*(Y(M+1,tt-1)*(Y(J*K,tt-1)- Y(M+2,tt-1) )-Y(M,tt-1)+hy*X(ii,tt-1));
            elseif ii == K && jj == J-1 
                Y(M,tt) = Y(M,tt-1) + (1/epsilon)*dt*(Y(M+1,tt-1)*(Y(M-1,tt-1)- Y(1  ,tt-1) )-Y(M,tt-1)+hy*X(ii,tt-1));
            elseif ii == K && jj == J 
                Y(M,tt) = Y(M,tt-1) + (1/epsilon)*dt*(Y(1  ,tt-1)*(Y(M-1,tt-1)- Y(2  ,tt-1) )-Y(M,tt-1)+hy*X(ii,tt-1));
            else
                Y(M,tt) = Y(M,tt-1) + (1/epsilon)*dt*(Y(M+1,tt-1)*(Y(M-1,tt-1)- Y(M+2,tt-1) )-Y(M,tt-1)+hy*X(ii,tt-1));
            end
            summ = summ + Y(M,tt); % sum over all y from 1 to J
        end
        MeanY(ii,tt) = (1/J)*summ;
        Bk = (hx/J)*summ; % transform into Bk
        nan = MeanY(ii,tt);
        if isnan(nan)
            break
        end
        % compute all X
        if ii == 1
            X(ii,tt) = X(ii,tt-1) + dt*(X(K   ,tt-1)*(X(ii+1,tt-1)-X(K-1 ,tt-1))-X(ii,tt-1) + F + Bk);
        elseif ii == 2
            X(ii,tt) = X(ii,tt-1) + dt*(X(ii-1,tt-1)*(X(ii+1,tt-1)-X(K   ,tt-1))-X(ii,tt-1) + F + Bk);
        elseif ii == K
            X(ii,tt) = X(ii,tt-1) + dt*(X(ii-1,tt-1)*(X(1   ,tt-1)-X(ii-2,tt-1))-X(ii,tt-1) + F + Bk);
        else
            X(ii,tt) = X(ii,tt-1) + dt*(X(ii-1,tt-1)*(X(ii+1,tt-1)-X(ii-2,tt-1))-X(ii,tt-1) + F + Bk);
        end
    end
end
toc
%% plot Xi and YiJ over time
kk=7; % pick a k
figure()
[~,B] = size(MeanY);
A=1;
%B=2676691;
xAs = A:B;
plot(xAs,X(kk,A:B))
hold on
plot(xAs,MeanY(kk,A:B))
legend('X','Y')
title('X and Y over time')
%% Histogram Xk
figure()
BinN = 50;  % number of bins
XLorenz4 = X(:,400000:TIME);
hh = histogram(XLorenz4 ,BinN, 'Normalization','pdf');
[values, edges] = histcounts(XLorenz4,BinN, 'Normalization', 'pdf');
centers = (edges(1:end-1)+edges(2:end))/2;
hold on
figure()
plot(centers, values, 'k-')
%% Make distributions of the mean value of Y_k for every X_k
BinDistr = 16; % bins as in Crommelin
T1 = 400000; % relaxation time
BkXval = zeros(BinDistr,TIME-T1);
TotalBIN = zeros(BinDistr,1);
Binteller = zeros(BinDistr,TIME-T1); % Counts number of X in the same bin 
for timeN = T1+1:TIME
    timeN
    for kk = 1:K
        BIN = floor(X(kk,timeN)+0.5);  % determine the bin for every Xkk
        meanJ = MeanY(kk,timeN); %(1/4)*sum(Y((8*(kk-1)+3):(8*(kk-1)+6),timeN)); % compute corresponding mean Y
        if BIN <= -5               % all bins smaller than -4.5 are one bin
            BkXval(1,timeN-T1) = BkXval(1,timeN-T1) + meanJ;
            Binteller(1,timeN-T1) = Binteller(1,timeN-T1) + 1;
        elseif BIN > 10            % all bins smaller than 9.5 are one bin
            BkXval(16,timeN-T1) = BkXval(16,timeN-T1) + meanJ;
            Binteller(16,timeN-T1) = Binteller(16,timeN-T1) + 1;
        else                       % all other bins
            BkXval(BIN+6,timeN-T1) = BkXval(BIN+6,timeN-T1) + meanJ;
            Binteller(BIN+6,timeN-T1) = Binteller(BIN+6,timeN-T1) + 1;
        end
    end
    % all values in BkXval are divided by their number of entries to
    % compute the mean of Y_k
    for kk2 = 1:K
        if Binteller(kk2,timeN-T1) ~= 0
            BkXval(kk2,timeN-T1) = BkXval(kk2,timeN-T1)/Binteller(kk2,timeN-T1);
            TotalBIN(kk2) = TotalBIN(kk2) + Binteller(kk2);
        end
    end
end
BkXval0 = BkXval;
BkXval(BkXval == 0) = NaN;
Distr = [BinDistr,3];  % Distributions(Xk, Mean(Bk), Std(Bk))
for distr = 1:BinDistr
    Distr(distr,1) = distr-6;   % values of the bins
    Mean = mean(BkXval(distr,:),'omitnan');  % compute the mean of each bin
    Std = std(BkXval(distr,:),'omitnan');    % standard deviation of each bin
    Distr(distr,2) = Mean;  
    Distr(distr,3) = Std;
end
%% Reduced model run, using the means of Y_k
% create all matrices
E = eye(K);
E1 = circshift(E,-1,2); % matrix for Xk-1
E2 = circshift(E,1,2); % matrix for Xk+1
E3 = circshift(E,-2,2); % matrix for Xk-2
% initialize X 
Xx2(:,1) = X(:,400000);  % -5 + (10)*rand(K,1);
% F(X) = E1*X*(E2-E3)*X-X+E*F+Bx
for ii = 2:TIME
    Bx = zeros(K,1);
    for jj = 1:K
        BIN = floor(Xx2(jj)-0.5);
        if BIN <= -5
            Bx(jj) = Distr(1,2); % Distributions wordt gemaakt in Lorenz4.m
        elseif BIN > 10
            Bx(jj) = Distr(16,2); % Distributions wordt gemaakt in Lorenz4.m
        else
            Bx(jj) = Distr(BIN+6,2); % Distributions wordt gemaakt in Lorenz4.m
        end
    end
    Xx2(:,ii) = Xx2(:,ii-1) + dt*(diag((E1*Xx2(:,ii-1))*((E2-E3)*Xx2(:,ii-1))') - Xx2(:,ii-1) + 1*ones(K,1)*F + hx*Bx);
end
figure()
for kk = 1:K
    plot(Xx2(kk,1:TIME-1))
    hold on
end
title('Xx2')
%%



