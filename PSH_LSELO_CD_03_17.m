tic
format short
clc;clear;close;
%% === model parameters ===
n = 400; p = 500;
N = 200;
rho = [0.25 0.5];
mu = zeros(p,1); a = (1:p);
ama = bsxfun(@minus,a,a'); 
rho = rho(1);
sigma = rho.^(abs(ama));


%% Case A: tau1-->censoring rate=25%
% tau1 = 1.6;  tau2 = 0.75 ;   % rho = 0.25
% tau1 = 1.65;  tau2 = 0.78;   % rho = 0.5

%% Case B: tau2-->censoring rate=25%
tau1 = 1.10;  tau2 = 0.45;   % rho = 0.25
% tau1 = 1.2;  tau2 = 0.5;   % rho = 0.5



%% === 'lambda' parameter ===
lambda_selo = 0.15*(0.001).^((1:100)/100);
lambda_MCP = 0.15*(0.001).^((1:100)/100);
lambda_scad = 0.15*(0.001).^((1:100)/100);
lambda_Alasso = 0.15*(0.001).^((1:100)/100);
lambda_lasso = 0.15*(0.001).^((1:100)/100);

% lambda_selo = 0.15*(0.001).^((1:50)/100);
% lambda_MCP = 0.15*(0.001).^((1:50)/100);
% lambda_scad = 0.15*(0.001).^((1:50)/100);
% lambda_Alasso = 0.15*(0.001).^((1:50)/100);
% lambda_lasso = 0.15*(0.001).^((1:50)/100);

%% === 'theta' parameter ===
theta1 = (1.35:-0.005:0.900)*10^2; 
%theta1 = (1.35:-0.002:0.900)*10^2; 

block_size = 20; p_size = 10;

%% === True Beta ===
% Beta = 0.5*[-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];   %% Case A
Beta = [-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)];       %% Case B

% Beta = [-0.7;-0.7;0;0;0;-0.7;zeros(p-6,1)];
% Beta = [-0.4;-0.4;0;0;0;-0.4;zeros(p-6,1)];

Censorrate = zeros(N,1);
index = find(Beta~=0); 
lselo_lqa = zeros(p,N); % gselo_lqa = zeros(p,N);
se_lselo = zeros(1,N); se_gselo =zeros(1,N);
opt_theta1 =zeros(1,N); opt_theta2 =zeros(1,N);

for iter = 1:N
    iter 
    rng(iter)   % % set seeds
    [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau2,iter);
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n);
    initial_beta(:,iter) = ini_beta(Z,1e-5,Delta,Iota,R,W,theta2,block_size);  % % initial beta

    % [lselo_lqa(:,iter),opt_theta2(iter)] = lselo_LQA(n,initial_beta(:,iter),Z,...
    %     1e-5,Delta,Iota,R,W,theta2);   % % LSELO
    [lselo_cd(:,iter),opt_theta3(iter)] = lselo_CD(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,theta2,block_size,p_size);
    % [selo_ista(:,iter),opt_lambda_selo(iter)] = ista_SELO(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_selo);    % % SELO 
    % MIC_ista(:,iter) = ista_MIC(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W);    % % MIC
    % [MCP_ista(:,iter),opt_lambda_MCP(iter)] = ista_MCP(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_MCP);       % % MCP
    % [scad_ista(:,iter),opt_lambda_scad(iter)] = ista_scad(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_scad);    % % SCAD
    % [Alasso_ista(:,iter),opt_lambda_Alasso(iter)] = ista_Alasso(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_Alasso);  % % Alasso
    % [lasso_ista(:,iter),opt_lambda_lasso(iter)] = ista_lasso(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,lambda_lasso);      % % Lasso

    % se_lselo(iter) = (lselo_lqa(:,iter)-Beta)'*sigma*(lselo_lqa(:,iter)-Beta);   % % mse_lselo
    % se_gselo(iter) = (gselo_lqa(:,iter)-Beta)'*(gselo_lqa(:,iter)-Beta);   % % mse_gselo
    se_lselo_cd(iter) = (lselo_cd(:,iter)-Beta)'*(lselo_cd(:,iter)-Beta);   % % mse_lselo_cd
    % se_selo(iter) = (selo_ista(:,iter)-Beta)'*sigma*(selo_ista(:,iter)-Beta);    % % mse-selo
    % se_MIC(iter) = (MIC_ista(:,iter)-Beta)'*sigma*(MIC_ista(:,iter)-Beta);       % % mse-MIC
    % se_MCP(iter) = (MCP_ista(:,iter)-Beta)'*sigma*(MCP_ista(:,iter)-Beta);       % % mse-MCP
    % se_scad(iter) = (scad_ista(:,iter)-Beta)'*sigma*(scad_ista(:,iter)-Beta);    % % mse-scad
    % se_Alasso(iter) = (Alasso_ista(:,iter)-Beta)'*sigma*(Alasso_ista(:,iter)-Beta); % % mse-Alasso
    % se_lasso(iter) = (lasso_ista(:,iter)-Beta)'*sigma*(lasso_ista(:,iter)-Beta);    % % mse-lasso
    
end



%% ----------------------------------------
%      Assessment Criteria :
% ---------------------------------------------------------------------------------------------------------------
% corr_lselo = sum((all(lselo_lqa(index,:))).*(1-any(lselo_lqa(setdiff(1:1:p, index),:))))/N;
% MSE_lselo = mean(se_lselo);
% N_plus_lselo = sum(sum(lselo_lqa(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_lselo = sum(sum(lselo_lqa(index,:)==0))/N;
% Size_lselo = sum(sum(lselo_lqa(:,:)~=0))/N;

corr_lselo_cd = sum((all(lselo_cd(index,:))).*(1-any(lselo_cd(setdiff(1:1:p, index),:))))/N;
MSE_lselo_cd = mean(se_lselo_cd);
N_plus_lselo_cd = sum(sum(lselo_cd(setdiff(1:1:p, index),:)~=0))/N;
N_minus_lselo_cd = sum(sum(lselo_cd(index,:)==0))/N;
Size_lselo_cd = sum(sum(lselo_cd(:,:)~=0))/N;

% corr_gselo = sum((all(gselo_lqa(index,:))).*(1-any(gselo_lqa(setdiff(1:1:p, index),:))))/N;
% MSE_gselo = mean(se_gselo);
% N_plus_gselo = sum(sum(gselo_lqa(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_gselo = sum(sum(gselo_lqa(index,:)==0))/N;
% Size_gselo = sum(sum(gselo_lqa(:,:)~=0))/N;

% corr_selo = sum((all(selo_ista(index,:))).*(1-any(selo_ista(setdiff(1:1:p, index),:))))/N;
% MSE_selo = mean(se_selo);
% N_plus_selo = sum(sum(selo_ista(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_selo = sum(sum(selo_ista(index,:)==0))/N;
% Size_selo = sum(sum(selo_ista(:,:)~=0))/N;
% 
% corr_MIC = sum((all(MIC_ista(index,:))).*(1-any(MIC_ista(setdiff(1:1:p, index),:))))/N;
% MSE_MIC = mean(se_MIC);
% N_plus_MIC = sum(sum(MIC_ista(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_MIC = sum(sum(MIC_ista(index,:)==0))/N;
% Size_MIC = sum(sum(MIC_ista(:,:)~=0))/N;
% 
% corr_MCP = sum((all(MCP_ista(index,:))).*(1-any(MCP_ista(setdiff(1:1:p, index),:))))/N;
% MSE_MCP = mean(se_MCP);
% N_plus_MCP = sum(sum(MCP_ista(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_MCP = sum(sum(MCP_ista(index,:)==0))/N;
% Size_MCP = sum(sum(MCP_ista(:,:)~=0))/N;
% 
% corr_scad = sum((all(scad_ista(index,:))).*(1-any(scad_ista(setdiff(1:1:p, index),:))))/N;
% MSE_scad = mean(se_scad);
% N_plus_scad = sum(sum(scad_ista(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_scad = sum(sum(scad_ista(index,:)==0))/N;
% Size_scad = sum(sum(scad_ista(:,:)~=0))/N;
% 
% corr_Alasso = sum((all(Alasso_ista(index,:))).*(1-any(Alasso_ista(setdiff(1:1:p, index),:))))/N;
% MSE_Alasso = mean(se_Alasso);
% N_plus_Alasso = sum(sum(Alasso_ista(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_Alasso = sum(sum(Alasso_ista(index,:)==0))/N;
% Size_Alasso = sum(sum(Alasso_ista(:,:)~=0))/N;
% 
% corr_lasso = sum((all(lasso_ista(index,:))).*(1-any(lasso_ista(setdiff(1:1:p, index),:))))/N;
% MSE_lasso = mean(se_lasso);
% N_plus_lasso = sum(sum(lasso_ista(setdiff(1:1:p, index),:)~=0))/N;
% N_minus_lasso = sum(sum(lasso_ista(index,:)==0))/N;
% Size_lasso = sum(sum(lasso_ista(:,:)~=0))/N;

%% output
 %Criteria = [corr MSE N_plus N_minus Size]
 Criteria = [%corr_lselo  MSE_lselo  N_plus_lselo  N_minus_lselo  Size_lselo; 
corr_lselo_cd  MSE_lselo_cd  N_plus_lselo_cd  N_minus_lselo_cd  Size_lselo_cd;
             % corr_gselo  MSE_gselo  N_plus_gselo  N_minus_gselo  Size_gselo;
             % corr_selo   MSE_selo   N_plus_selo   N_minus_selo   Size_selo;
             % corr_MIC    MSE_MIC    N_plus_MIC    N_minus_MIC    Size_MIC;
             % corr_MCP    MSE_MCP    N_plus_MCP    N_minus_MCP    Size_MCP;
             % corr_scad   MSE_scad   N_plus_scad   N_minus_scad   Size_scad;
             % corr_Alasso MSE_Alasso N_plus_Alasso N_minus_Alasso Size_Alasso;
             % corr_lasso  MSE_lasso  N_plus_lasso  N_minus_lasso  Size_lasso;
             ];


% % ----save data---- %  %
% fid = fopen('ARM.txt', 'w');
% fid = fopen('PSH_VS.txt', 'a');
% fprintf(fid,'%9.4s\t', 'n', 'p', 'rho','censoring','Beta');
% fprintf(fid, '\n' );
% fprintf(fid, '%9.4f\t',[n;p;rho;mean(Censorrate);Beta(index)]);
% fprintf(fid, '\n\n' );
% fprintf(fid,'%9.4s\t', 'LSELO', 'GSELO', 'SELO', 'MIC','MCP','SCAD','ALASSO','LASSO');
% fprintf(fid, '\n' );
% fprintf(fid, [repmat('%9.4f\t', 1 ,size(Criteria',2)), '\n'],Criteria);
% fprintf(fid, '\n\n' );
% fclose(fid);
% open('PSH_VS.txt')
Criteria
censorrate = mean(Censorrate)

time = toc  %runtime

%% ============================================================
%                  SUBFUNCTIONS
% ===============%==%==%==%==%==%==%===========================
%                 survival_data()
% ===============%==%==%==%==%==%==%===========================
function [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau,iter)
%Z:n*p
%X:n*1
%R:n*n
%p = length(Beta);
Z = mvnrnd(mu,sigma,n);
Beta2 = -Beta;

F = unifrnd(0,1,[n,1]);
lambda = exp(Z*Beta)+exp(Z*Beta2);
T = -log(1-F)./lambda;  %T
C = unifrnd(0,tau,[n,1]);
Delta = (T<=C);
X = min(T,C);

prop = exp(Z*Beta)./lambda; 
alphabet = [1 0];
Iota = zeros(n,1);
for i=1:n
    Iota(i,1) = randsrc(1,1,[alphabet; prop(i,1) 1-prop(i,1)]);
end

%**************************************************
% At risk set：R()
%**************************************************
R = zeros(n,n);
for i=1:n                        
    for j=1:n
        if X(j,1)>=X(i,1)||(X(j,1)<=X(i,1)&&Delta(j,1)==1&&Iota(j,1)~=1)
            R(j,i)=1;
        end
    end
end

end


%% ===================================================
%                 Weight()
% ============================================================
function W = Weight(X,T,C,Delta,n)
W = zeros(n,n);
% C_Delta = 1-Delta;
[G,id] = C_KM(X,Delta);

for i=1:n
    for j=1:n
        if X(i) < X(j)
            W(j,i) = (C(i)>=min(T(i),X(j)))*G(id(j))/(G(id(i))+1e-6);
        else
            W(j,i) = C(i)>=min(T(i),X(j));
        end

    end
end

end

%% ===================================================
%                 C_KM()
% ============================================================
function [G,id]=C_KM(X,Delta)
% [f,x] = ecdf(C,'Censoring',Delta); %KM estimator of C 
% G = 1-f;
% ecdf(C,'censoring',Delta,'function','survivor'); %survival function of C

[n,~] = size(X);
cens = 1-Delta;
[~,id] = sort(X);
cens = cens(id);  % indicator--censoring

temp_G=1; G = zeros(n,1);

for i=1:n 
    if cens(i,1)==1 && (n-i~=0)
        temp_G = temp_G*(1-1/(n-i+1));
        G(i,1) = temp_G;
    elseif cens(i,1)==0 && (n-i~=0)
        G(i,1) = temp_G;
    
    elseif cens(i,1)==1 && (n-i==0)
        G(i,1) = 0;
    else
        G(i,1) = temp_G;
    end
end


end


%% ===================================================
%                 ini_beta()
% ============================================================
function initial_beta = ini_beta(Z,r,Delta,Iota,R,W,theta,block_size)
[n,p] = size(Z);
beta = zeros(p,1);
f2 = zeros(n,p);

n0 = sum(Delta);
lambda0 = log(n0);
opt_BIC = 1e+10;


%% Block coordinate descent parameter settings
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% index--blocks
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end


k = 1;err = 0; tk = 24;
while k<=1000&&err==0
    % k
    % for i=1:n
    %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
    %     %tem1 = R(:,i).*exp(Z*beta);
    %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
    % end

    % revised on 17/03/2025
    tem1 = (W.*R).*exp(Z*beta);
    sum_tem1 = sum( tem1, 1 );
    f2 = (Z'*tem1)./sum_tem1;
    f2 = f2';

    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

    %% Block coordinate descent
    if p<=10
        beta1 =  beta - L_prime/tk;
    else
        for iter = 1:50
            % Random block update sequence
            block_order = randperm(num_blocks);
            %%
            for b = block_order
                current_block = blocks{b};
                beta_block = beta(current_block);
                beta_new1 = beta_block - L_prime(current_block)/tk;

                beta1(current_block,1) = beta_new1;

            end

            delta = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);
            if delta < 1e-6
                break;
            end
        end
    end

    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;

    % initial_beta = beta;

end

for j = 1:length(theta)

k=1;err=0; tk = 16;
while k<=1000 && err==0
    %k
    W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
    u = eye(p) + 2*W1/tk;

    % for i=1:n
    %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
    %     %tem1 = R(:,i).*exp(Z*beta);
    %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
    % end

    % revised on 17/03/2025
    tem1 = (W.*R).*exp(Z*beta);
    sum_tem1 = sum( tem1, 1 );
    f2 = (Z'*tem1)./sum_tem1;
    f2 = f2';

    L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

    beta_tilde = beta - L_prime/tk;
    beta1 = u\beta_tilde;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end

beta2 = beta.*(abs(beta)>=2*1e-4);

% for i=1:n
%     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
% end

% revised on 17/03/2025
tem2 = sum( (W.*R).*exp(Z*beta),1 );

ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
sel = beta2~=0;
BIC = ell+sum(sel)*log(n)/n;

if BIC<=opt_BIC
    opt_BIC = BIC;
    opt_beta = beta2;
end

end

initial_beta = opt_beta;


end

%% ===================================================
%                 lselo_LQA()
% ============================================================
function [opt_beta,opt_theta] = lselo_LQA(n,ini_beta,Z,r,Delta,Iota,R,W,theta)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta = ini_beta;  % % initial value
% beta = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

for j = 1:length(theta)
    k=1;err=0; tk = 2; beta = ini_beta;
    while k<=1000 && err==0
        %k
        W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
        u = eye(p) + 2*W1/tk;

        % for i=1:n
        %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
        %     %tem1 = R(:,i).*exp(Z*beta);
        %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
        % end

        % revised on 17/03/2025
        tem1 = (W.*R).*exp(Z*beta);
        sum_tem1 = sum( tem1, 1 );
        f2 = (Z'*tem1)./sum_tem1;
        f2 = f2';

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
%         beta_tilde = beta - L_prime/(n*tk);
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
   
    beta2 = beta.*(abs(beta)>=2*1e-4);

    % for i=1:n
    %     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
    % end

    % revised on 17/03/2025
    tem2 = sum( (W.*R).*exp(Z*beta),1 );

    ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;

    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_theta = theta(j);
        opt_beta = beta2;
    end

end

end





%% ===================================================
%                 lselo_CD()
% ============================================================
function [opt_beta,opt_theta] = lselo_CD(n,ini_beta,Z,r,Delta,Iota,R,W,theta,block_size,p_size)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
% beta = ini_beta;  % % initial value
% beta1 = zeros(p,1);
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

% a = 0.00059;
% a = 0.00025;  % learning rate


%% Block coordinate descent parameter settings
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

%% index--blocks
for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end

for j = 1:length(theta)
    k=1; err=0; tk = 8; beta = ini_beta;
    while k<=1000 && err==0
        %k
        W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
        u = eye(p) + 2*W1/tk;

        % for i=1:n
        %     tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
        %     %tem1 = R(:,i).*exp(Z*beta);
        %     f2(i,:) = sum(tem1.*Z)/sum(tem1);
        % end

        % revised on 17/03/2025
        tem1 = (W.*R).*exp(Z*beta);
        sum_tem1 = sum( tem1, 1 );
        f2 = (Z'*tem1)./sum_tem1;
        f2 = f2';

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

        %% Block coordinate descent
        if p<=p_size
            beta_tilde = beta - L_prime/tk;
            % beta_tilde = beta - L_prime/(n*tk);
            beta1 = u\beta_tilde;
        else
            for iter = 1:50
                % Random block update sequence
                block_order = randperm(num_blocks);
                %%
                for b = block_order
                    current_block = blocks{b};
                    beta_block = beta(current_block);

                    % calculate W1
                    W12 = lambda0*theta(j)*diag(1./(1+theta(j)*beta_block.^2).^2);             
                    u1 = eye(length(current_block)) + 2*W12/tk;

                    beta_tilde = beta_block - L_prime(current_block)/tk;
                    beta_new1 = u1\beta_tilde;

                    beta1(current_block,1) = beta_new1;

                end

                delta = norm(beta - beta1,2)^2 / (norm(beta,2)^2 + 1e-8);
                if delta < 1e-6
                    break;
                end
            end
        end

        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;

        k = k+1;

    end
   
    beta2 = beta.*(abs(beta)>=2*1e-4);

    % for i=1:n
    %     tem2(i,1) = sum((W(:,i).*R(:,i)).*exp(Z*beta2));
    % end

    % revised on 17/03/2025
    tem2 = sum( (W.*R).*exp(Z*beta),1 );

    ell = -sum((Delta.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;

    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_theta = theta(j);
        opt_beta = beta2;
    end

end


end



%% ===================================================
%                 ista_lasso()
% ============================================================
function [opt_beta,opt_lambda] = ista_lasso(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
beta = ini_beta;
opt_BIC = 1e+10;

for i = 1:length(lambda)
    k = 1;err = 0; tk = 4; % beta = ini_beta;
    while k<=1000&&err==0
        %K
        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   
        L_prime = -sum((status.*Iota).*(Z-f2))'/n;
        
        u = beta - L_prime/tk;
        beta1 = max(zeros(p,1),abs(u)-lambda(i)/tk).*sign(u);
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta2;
    end
end

end

%% ===================================================
%                 ista_Alasso()
% ============================================================
function [opt_beta,opt_lambda] = ista_Alasso(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
WW = 1./(abs(ini_beta)+1e-6);
beta = ini_beta;
opt_BIC = 1e+10;
% beta3 = zeros(p,1);

for i = 1:length(lambda)
    k = 1;err = 0;tk = 4; % beta = ini_beta;
    while k<=1000&&err==0
        %       k
        % % === solving by Lasso form ===
        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   
        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        u = beta - L_prime/tk;
        beta1 = max(zeros(p,1),abs(u)-lambda(i).*WW/tk).*sign(u);      
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k=k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta1;
    end
    
end

end



%% ===================================================
%                 ista_scad()
% ============================================================
function [opt_beta,opt_lambda] = ista_scad(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
beta = ini_beta;
opt_BIC=1e+10; a = 3.7;

for j=1:length(lambda)
    k=1;err=0; tk = 2; % beta = ini_beta;
    while k<=1000 && err==0
%         k
        W1 = diag( ( lambda(j)*(abs(beta)<=lambda(j)) +...
            ((lambda(j)<abs(beta))&(abs(beta)<a*lambda(j))).*((a*lambda(j)-abs(beta))/(a-1)) )./(abs(beta)+1e-6) );
        u = eye(p) + W1/tk;

        for i=1:n
            tem1 = (W(:,i).*R(:,i)).*exp(Z*beta);
            f2(i,:) = sum(tem1.*Z)/sum(tem1);
        end   
        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(j);
        opt_beta = beta2;
    end
    
end

end


%% ===================================================
%                 ista_MCP()
% ============================================================
function [opt_beta,opt_lambda] = ista_MCP(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
beta = ini_beta;
tv = 2.7;    % according to Cao et al.(2017)
opt_BIC=1e+10;

for i=1:length(lambda)
    k=1;err=0; tk = 4; % beta = ini_beta;
    while k<=1000 && err==0
%        k
        W1=diag(max(0,(tv*lambda(i)-abs(beta))/(tv*lambda(i)) )./(abs(beta)+1e-6));
        u = eye(p) + W1/tk;

        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   

        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end
    
    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta2;
    end
end

end


%% ===================================================
%                 ista_SELO()
% ============================================================
function [opt_beta,opt_lambda] = ista_SELO(n,ini_beta,Z,r,status,Iota,R,W,lambda)
[~,p] = size(Z);
f2 = zeros(n,p);
tem2 = zeros(n,1);
opt_beta = zeros(p,1);
beta = ini_beta ;
% tau = 1.5*0.01;
tau = 0.01;    % according to Dicker et al.(2013)
opt_BIC=1e+10;

for i=1:length(lambda)
    k=1;err=0; tk = 4; % beta = ini_beta;
    while k<=1000 && err==0
        W1=diag((lambda(i)*tau/(2*log(2)))./((abs(beta)+1e-6).*(2*beta.^2+3*tau*abs(beta)+tau^2)));
        u = eye(p) + 2*W1/tk;

        for j=1:n
            tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
            f2(j,:) = sum(tem1.*Z)/sum(tem1);
        end   

        L_prime = -sum((status.*Iota).*(Z-f2))'/n;

        beta_tilde = beta - L_prime/tk;
        beta1 = u\beta_tilde;
        w = beta1-beta;
        err = norm(w,2)^2 <= r*norm(beta,2)^2;
        beta = beta1;
        k = k+1;
    end

    beta2 = beta.*(abs(beta)>=2*1e-4);

    for t=1:n
        tem2(t,1) = sum((W(:,t).*R(:,t)).*exp(Z*beta2));
    end

    ell = -sum((status.*Iota).*(Z*beta2-log(tem2)))/n;
    sel = beta2~=0;
    BIC = ell+sum(sel)*log(n)/n;
    if BIC<=opt_BIC
        opt_BIC = BIC;
        opt_lambda = lambda(i);
        opt_beta = beta2;
    end
end

end



%% ===================================================
%                 ista_MIC()
% ============================================================
function opt_beta = ista_MIC(n,ini_beta,Z,r,status,Iota,R,W)
[~,p] = size(Z);
f2 = zeros(n,p);
beta = ini_beta;
n0 = sum(status);
lambda0 = log(n0);
a = 70;
%a = n0;

k=1;err=0; tk = 2;
while k<=1000 && err==0
    %   k
    W1 = diag( (4*a*lambda0*exp(2*a*beta.^2))./((exp(2*a*beta.^2)+1).^2) );
    u = eye(p) + 2/tk*W1;
    for j=1:n
        tem1 = (W(:,j).*R(:,j)).*exp(Z*beta);
        f2(j,:) = sum(tem1.*Z)/sum(tem1);
    end  

    L_prime = -sum((status.*Iota).*(Z-f2))'/n;

    beta_tilde = beta - L_prime/tk;
    beta1 = u\beta_tilde;
    w = beta1-beta;
    err = norm(w,2)^2 <= r*norm(beta,2)^2;
    beta = beta1;
    k = k+1;
end
opt_beta = beta.*(abs(beta)>2*1e-4);


end




