tic
format short
clc;clear;close;
%% === model parameters ===
n = 400; p = 100;
N = 200; % Replication
rho = [0.25 0.5];
mu = zeros(p,1); a = (1:p);
ama = bsxfun(@minus,a,a'); 
rho = rho(1);
sigma = rho.^(abs(ama));

tau1 = 1.10;  tau2 = 0.45;  
theta = (0.60:-0.005:0.50)*1.25*10^4;
block_size = 20; p_size = 10;

%% === True Beta ===
Beta = [-1.0;1.0;0;0;0;-1.0;zeros(p-6,1)]; 


Censorrate = zeros(N,1);
index = find(Beta~=0); 
lselo_lqa = zeros(p,N); 
se_lselo = zeros(1,N); se_gselo =zeros(1,N);
opt_theta1 =zeros(1,N); opt_theta2 =zeros(1,N);

for iter = 1:N
    iter 
    rng(iter)   % random seed
    [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau2,iter);
    Censorrate(iter) = 1-mean(Delta);
    W = Weight(X,T,C,Delta,n);
    initial_beta(:,iter) = ini_beta(Z,1e-5,Delta,Iota,R,W,theta,block_size);  % % initial beta
    [lselo_cd(:,iter),opt_theta3(iter)] = lselo_CD(n,initial_beta(:,iter),Z,1e-5,Delta,Iota,R,W,theta,block_size,p_size);
    se_lselo_cd(iter) = (lselo_cd(:,iter)-Beta)'*(lselo_cd(:,iter)-Beta);   % % mse_lselo_cd
    
end



%% ----------------------------------------
%      Assessment Criteria :
% ---------------------------------------------------------------------------------------------------------------

corr_lselo_cd = sum((all(lselo_cd(index,:))).*(1-any(lselo_cd(setdiff(1:1:p, index),:))))/N;
MSE_lselo_cd = mean(se_lselo_cd);
N_plus_lselo_cd = sum(sum(lselo_cd(setdiff(1:1:p, index),:)~=0))/N;
N_minus_lselo_cd = sum(sum(lselo_cd(index,:)==0))/N;
Size_lselo_cd = sum(sum(lselo_cd(:,:)~=0))/N;

%% output

 Criteria = [corr_lselo_cd MSE_lselo_cd N_plus_lselo_cd N_minus_lselo_cd Size_lselo_cd];

Criteria

censorrate = mean(Censorrate)

time = toc  %运行时间

%% ===============%==%==%==%==%==%==%===========================
%                  SUBFUNCTIONS
% ===============%==%==%==%==%==%==%===========================



%% ===================================================
%                 survival_data()
% ============================================================
function [Z,X,T,C,Iota,Delta,R] = survival_data(n,Beta,mu,sigma,tau,iter)
%Z:n*p
%X:n*1
%R:n*n
%p = length(Beta);
Z = mvnrnd(mu,sigma,n);
Beta2 = -Beta;

F = unifrnd(0,1,[n,1]);
lambda = exp(Z*Beta)+exp(Z*Beta2);
T = -log(1-F)./lambda;  
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
% ：R()
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

[n,~] = size(X);
cens = 1-Delta;
[~,id] = sort(X);
cens = cens(id);  

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
% f2 = zeros(n,p);

n0 = sum(Delta);
lambda0 = log(n0);
opt_BIC = 1e+10;


%% Block-wise CD 
% block_size = 10;
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

for b = 1:num_blocks
    start_idx = (b-1)*block_size + 1;
    end_idx = min(b*block_size, p);
    blocks{b} = start_idx:end_idx;
end

k = 1;err = 0; tk = 24;
while k<=1000&&err==0
    % k

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
            block_order = randperm(num_blocks);
            
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
end

for j = 1:length(theta)

k=1;err=0; tk = 16;
while k<=1000 && err==0
    %k
    W1 = lambda0*theta(j)*diag(1./(1+theta(j)*beta.^2).^2);
    u = eye(p) + 2*W1/tk;

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
%                 lselo_CD()
% ============================================================
function [opt_beta,opt_theta] = lselo_CD(n,ini_beta,Z,r,Delta,Iota,R,W,theta,block_size,p_size)
[~,p] = size(Z);
% f2 = zeros(n,p);
% tem2 = zeros(n,1);
opt_beta = zeros(p,1);
% beta = ini_beta;  
opt_BIC = 1e+10;
n0 = sum(Delta);
lambda0 = log(n0);

%% Block-wise CD
num_blocks = ceil(p / block_size);
blocks = cell(num_blocks, 1);

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

        tem1 = (W.*R).*exp(Z*beta);
        sum_tem1 = sum( tem1, 1 );
        f2 = (Z'*tem1)./sum_tem1;
        f2 = f2';

        L_prime = -sum((Delta.*Iota).*(Z-f2))'/n;

        %% Block coordinate descent
        if p<=p_size
            beta_tilde = beta - L_prime/tk;
            beta1 = u\beta_tilde;
        else
            for iter = 1:50
                block_order = randperm(num_blocks);

                for b = block_order
                    current_block = blocks{b};
                    beta_block = beta(current_block);

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


