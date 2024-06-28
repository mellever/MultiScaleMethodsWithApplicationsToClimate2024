%%  Make distributions Bk for every Xk coupled to Xk-2, Xk-1 and Xk+1
% ingredients; X, YMean from Lorenz96.m
% makes BkXvalues, forcing term for the reduced model
K = 9;
J = 8;
BinN = 4;
Xvalues = zeros(BinN,1);
BkXvalues4 = zeros(BinN,BinN,BinN,BinN,2);
Teller4 = zeros(BinN,BinN,BinN,BinN,1); % keeps track of the number of entries 
for timeN = 400000:TIME
    Location = zeros(9,4);
    teller = 1;
% cycle for k = 1 
    kk = 1;
    BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN1 <=-2
        BIN1 = -1;
    end
    BIN2 = min(floor((X(K,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN2 <=-2
        BIN2 = -1;
    end
    BIN3 = min(floor((X(K-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN3 <=-2
        BIN3 = -1;
    end
    BIN4 = min(floor((X(kk+1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN4 <=-2
        BIN4 = -1;
    end
    Xvalues(BIN1+2,1) = Xvalues(BIN1+2,1)+1; % Count number of Xk per bin
    bk = MeanY(kk,timeN);
    BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = bk;
    Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
    teller = teller+1;
% cycle for k = 2 
    kk = 2;
    BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN1 <=-2
        BIN1 = -1;
    end
    BIN2 = min(floor((X(kk-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN2 <=-2
        BIN2 = -1;
    end
    BIN3 = min(floor((X(K,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN3 <=-2
        BIN3 = -1;
    end
    BIN4 = min(floor((X(kk+1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN4 <=-2
        BIN4 = -1;
    end
    Xvalues(BIN1+2,1) = Xvalues(BIN1+2,1)+1; % Count number of Xk per bin
    bk = MeanY(kk,timeN);
    if BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) == 0
        BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = bk;
        Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
        teller = teller+1;
    else
        BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = (BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1)+bk)/2;
        Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
        teller = teller+1;
    end
% cycle for k = 3 ... k = K-1
    for kk = 3:(K-1)
        BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN1 <=-2
            BIN1 = -1;
        end
        BIN2 = min(floor((X(kk-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN2 <=-2
            BIN2 = -1;
        end
        BIN3 = min(floor((X(kk-2,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN3 <=-2
            BIN3 = -1;
        end
        BIN4 = min(floor((X(kk+1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN4 <=-2
            BIN4 = -1;
        end
        Xvalues(BIN1+2,1) = Xvalues(BIN1+2,1)+1; % Count number of Xk per bin
        bk = MeanY(kk,timeN);
        if BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) == 0
            BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = bk;
            Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
            teller = teller+1;
        else
            BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = (BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1)+bk)/2;
            Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
            teller = teller+1;
        end 
    end
    % cycle for k=K
    kk = 9;
    BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN1 <=-2
        BIN1 = -1;
    end
    BIN2 = min(floor((X(kk-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN2 <=-2
        BIN2 = -1;
    end
    BIN3 = min(floor((X(kk-2,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN3 <=-2
        BIN3 = -1;
    end
    BIN4 = min(floor((X(1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN4 <=-2
        BIN4 = -1;
    end
    Xvalues(BIN1+2,1) = Xvalues(BIN1+2,1)+1; % Count number of Xk per bin
    bk = MeanY(kk,timeN);
    if BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) == 0
        BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = bk;
        Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
    else
        BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1) = (BkXvalues4(BIN1+2,BIN2+2,BIN3+2,BIN4+2,1)+bk)/2;
        Location(teller,:) = [BIN1+2 BIN2+2 BIN3+2 BIN4+2];
    end
    % average the mean in each of the used bins
    for loop = 1:teller
        Abin = Location(loop,:);
        check = Teller4(Abin(1),Abin(2),Abin(3),Abin(4),1);
        if check == 0
            BkXvalues4(Abin(1),Abin(2),Abin(3),Abin(4),2) = BkXvalues4(Abin(1),Abin(2),Abin(3),Abin(4),1);
            Teller4(Abin(1),Abin(2),Abin(3),Abin(4),1) = Teller4(Abin(1),Abin(2),Abin(3),Abin(4),1) + 1;
        else
            BkXvalues4(Abin(1),Abin(2),Abin(3),Abin(4),2) = (check*BkXvalues4(Abin(1),Abin(2),Abin(3),Abin(4),2)...
                +BkXvalues4(Abin(1),Abin(2),Abin(3),Abin(4),1))/(check+1);
            Teller4(Abin(1),Abin(2),Abin(3),Abin(4),1) = check+1;
        end
    end
end
%% Run Xx3
% Function F(X), reduced model, taking X_k-2, X_k-1 and X_k+1 into account
% uses input BkXvalues4
%% parameters
epsilon = 1/128;
K = 9;
J = 8;
F = 10;
hx = -0.8;
hy = 1;
dt = 2^(-8);
T = 2^10;
TIME = floor(T/(2^(-14)));
% create all matrices
E = eye(K);
E1 = circshift(E,-1,2);
E2 = circshift(E,1,2);
E3 = circshift(E,-2,2);
% initialize X
Xx3_4_8(:,1) = X(:,400000);
BkXvalues4 = zeros(4,4,4,4,2);
for timeN = 2:TIME
    Bx = zeros(K,1);
% cycle for k = 1 
    kk = 1;
    BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN1 <=-2
        BIN1 = -1;
    end
    BIN2 = min(floor((X(K,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN2 <=-2
        BIN2 = -1;
    end
    BIN3 = min(floor((X(K-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN3 <=-2
        BIN3 = -1;
    end
    BIN4 = min(floor((X(kk+1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN4 <=-2
        BIN4 = -1;
    end
    Bx(1) = BkXvalues4(BIN1+2, BIN2+2, BIN3+2, BIN4+2, 2);
% cycle for k = 2
    kk = 2;
    BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN1 <=-2
        BIN1 = -1;
    end
    BIN2 = min(floor((X(kk-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN2 <=-2
        BIN2 = -1;
    end
    BIN3 = min(floor((X(K,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN3 <=-2
        BIN3 = -1;
    end
    BIN4 = min(floor((X(kk+1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN4 <=-2
        BIN4 = -1;
    end
    Bx(2) = BkXvalues4(BIN1+2, BIN2+2, BIN3+2, BIN4+2, 2);
% cycle for k = 3 ... k = K-1
    for kk = 3:(K-1)
        BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN1 <=-2
            BIN1 = -1;
        end
        BIN2 = min(floor((X(kk-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN2 <=-2
            BIN2 = -1;
        end
        BIN3 = min(floor((X(kk-2,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN3 <=-2
            BIN3 = -1;
        end
        BIN4 = min(floor((X(kk+1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
        if BIN4 <=-2
            BIN4 = -1;
        end
        Bx(kk) = BkXvalues4(BIN1+2, BIN2+2, BIN3+2, BIN4+2, 2);
    end
    % cycle for k=K
    kk = K;
    BIN1 = min(floor((X(kk,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN1 <=-2
        BIN1 = -1;
    end
    BIN2 = min(floor((X(kk-1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN2 <=-2
        BIN2 = -1;
    end
    BIN3 = min(floor((X(kk-2,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN3 <=-2
        BIN3 = -1;
    end
    BIN4 = min(floor((X(1,timeN)/4)+0.375),2);  % determine the bin for every Xkk
    if BIN4 <=-2
        BIN4 = -1;
    end
    Xx3_4_8(:,timeN) = Xx3_4_8(:,timeN-1) + dt*( diag((E1*Xx3_4_8(:,timeN-1))*((E2-E3)*Xx3_4_8(:,timeN-1))') - Xx3_4_8(:,timeN-1) + 1*ones(K,1)*F + hx*Bx);
end
%%
figure()
for kk = 1:K
    plot(Xx3_4_8(kk,1:TIME-1))
    hold on
end
title('Xx3')
