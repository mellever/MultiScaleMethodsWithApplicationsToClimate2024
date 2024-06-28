%% Lorenz 96 monoscale model
%% parameters
epsilon = 1/128;
K = 9;
F = 10;
dt = 2^(-14);
T = 2^10;
TIME = floor(T/(2^-14));
%% create all matrices
E = eye(K);
E1 = circshift(E,-1,2);
E2 = circshift(E,1,2);
E3 = circshift(E,-2,2);
%% initialize X
Monoscale(:,1) = X(:,400000); % Xinitial; %-5 + (10)*rand(K,1);
%% F(X) = E1*X*(E2-E3)*X-X+E*F
for ii = 2:TIME
    Monoscale(:,ii) = Monoscale(:,ii-1) + dt*( diag((E1*Monoscale(:,ii-1))*((E2-E3)*Monoscale(:,ii-1))') - Monoscale(:,ii-1) + 1*ones(K,1)*F);
end

