function VecFld = LLTA(X, Y, conf)
%% Initialization
gamma = conf.gamma;
beta = conf.beta;
lambda = conf.lambda;
theta = conf.theta;
a = conf.a;
MaxIter = conf.MaxIter;
ecr = conf.ecr;
minP = conf.minP;
M = conf.M;
Kn = conf.Kn;

%% Initiaization R = I, t = 0, P = P = I_N*N
[N, D]=size(X);
A = eye(D);
t = zeros(D,1);
s = 1;
P = eye(N);
iter=1;  tecr=1; E=1;     %ntol = tol+10, L = 0
sigma2=sum(sum((Y-X).^2))/(N*D);
%%
% Search the k nearest neighbors for each point X
X2= sum(X.^2,2);   %N*1
distance = repmat(X2,1,N)+repmat(X2',N,1)-2*X*X';
[~, index] = sort(distance);
neighborhood = index(2:(1+Kn),:); %K*N

%%
% Compute W by minimizeing the cost function E_LLE(W) = sum(square(xi - sum(Wij*xj)))
if(Kn>D) 
  tol=1e-5; % regularlizer in case constrained fits are ill conditioned
else
  tol=0;
end

W = zeros(Kn,N);%wij represprent xi  is the set of neighbors of xj and the wight is wij
WW = sparse(1:N,1:N,zeros(1,N),N,N,4*Kn*N);%N*N
IsubW  = sparse(1:N,1:N,ones(1,N),N,N,4*Kn*N);%N*N
for i = 1:N
    z = X(neighborhood(:,i),:) - repmat(X(i,:),Kn,1); % shift ith pt to origin  K*D
    G = z*z';                                         % local covariance  K*K
    G = G + eye(Kn,Kn)* tol * trace(G);                 % regularlization
    W(:,i) = G\ones(Kn,1);                             % solve Gw = 1
    W(:,i) = W(:,i)/sum(W(:,i));                     % normalize
    w = W(:,i);
    j = neighborhood(:,i);
    WW(i,j) = w';
    IsubW(i,j) = IsubW(i,j) - w';
end


%%
% EM iteration
PS = [];
while (iter<MaxIter) && (abs(tecr) > ecr) && (sigma2 > 1e-8) 
    %% E-step.
    % Update P
    
    E_old = E;
    T = repmat(t, [1,N]);
    TX = A*X'+ T; %D*N
    [P1, E] = get_P(Y, TX', sigma2 ,gamma, a);  
    PS = [PS, P1];
    P1 = max(P1, minP);
    P = sparse(1:N,1:N,P1,N,N,N);%N*N
    Np = sum(P1);
    muX = X'* P1 / Np;
    muY = Y'* P1 / Np;

    YtPX =  Y'* P * X - muY * muX' * Np;
    Q = sparse(1:N,1:N,ones(1,N),N,N,4*Kn*N);%N*N 
    Q = IsubW'*P*IsubW;
    XtQX = X'*Q*X;
    A = (Y'* P * X - muY * muX' * Np) / (2 * lambda * sigma2 * XtQX  +  X'*P*X - muX * muX' * Np );
    t = muY - A * muX;
    %update E 
    E = E + lambda * norm(A * X' * sqrt(P)* IsubW','fro');
    tecr=(E-E_old)/E;
    %% M-step.
    % update C by solving linear system
    % update sigma^2 and gamma by tr(v_TPv)/D.tr(P)  and tr(P)/N
    T = repmat(t, [1,N]);
    TX = A*X'+ T;
    V = Y - TX';
    sigma2 = sum(sum(V.^2.*repmat(P1,1,D)))/(D*Np);
    %sigma2 = trace(V'*P*V)/(D*trace(P));
    numcorr = length(find(P > theta));
    gamma=numcorr/N;
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    iter=iter+1;    
end
%%
VecFld.X = X;
VecFld.Y = Y;
VecFld.A = A;
VecFld.t = t;
VecFld.beta = beta;
VecFld.TX= TX';
VecFld.P = diag(P);
VecFld.PS = PS;
VecFld.VFCIndex = find(VecFld.P > theta);

% disp('Removing outliers succesfully completed.');


%%%%%%%%%%%%%%%%%%%%%%%%
function K=con_K(x,y,beta)
% CON_K constructs the kernel K, 
%   where K(i, j) = k(x, y) = exp(-beta*||x-y||^2).

[n, d]=size(x); [m, d]=size(y);

K=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
K=squeeze(sum(K.^2,2));
K=-beta * K;
K=exp(K);

%%%%%%%%%%%%%%%%%%%%%%%%
function [P, E]=get_P(Y, TX, sigma2 ,gamma, a)
% GET_P estimates the posterior probability and part of the energy.

D = size(Y, 2);
temp1 = exp(-sum((Y-TX).^2,2)/(2*sigma2));
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)/(gamma*a);
P=temp1./(temp1+temp2);
E=P'*sum((Y-TX).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2 - log(gamma)*sum(P) - log(1-gamma)*sum(1-P);
