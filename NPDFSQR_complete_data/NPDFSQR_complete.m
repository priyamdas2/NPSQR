clear all
rng(1)
data = xlsread('transformed_data.xls');
data_size = max(size(data));
x_array = data(:,1);% x_array should contain values between 0 and 1 (after required transformation)
y_array = data(:,2);% y_array should contain values between 0 and 1 (after required transformation)
mp = 4; % (value of m and p for B-spline basis expansion)
mcmc_no = 2000; % number of iterartions
mean_start = 1000; % Burn-in
r = 1.05; % movement coefficient


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Warm-start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GCDVSMS parameters %%%%%

no_iter = 5000;
tol_fun = 10^(-2);
tol_fun_2 = 10^(-2);
parameter_cut_off = 10^(-3);
epsilon_cut_off = 10^(-2);
epsilon_decreasing_factor_1 = 2; % default is 2
epsilon_decreasing_factor_2 = 1.05;
max_runs = 200;


m = mp; % quantile function resolution
p = mp; % spatial resolution
no_parameters = (mp+2)*(mp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basis Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = zeros(p+2,3*p);


c(1,1) = p^2;
c(1,2) = -2*p;
c(1,3) = 1;
c(2,1) = -(3/2)*p^2;
c(2,2) = 2*p;
c(2,3) = 0;
c(2,4) = p^2/2;
c(2,5) = -2*p;
c(2,6) = 2;
c(p+2,3*p-2) = p^2;
c(p+2,3*p-1) = -2*p*(p-1);
c(p+2,3*p) = (p-1)^2;
c(p+1,3*p-5) = p^2/2;
c(p+1,3*p-4) = -p*(p-2);
c(p+1,3*p-3) = ((p-2)^2)/2;
c(p+1,3*p-2) = -(3/2)*p^2;
c(p+1,3*p-1) = p*(3*p-2);
c(p+1,3*p) = -(3*p^2-4*p)/2;

for i=0:(p-3)
    c(i+3,3*i+1) = p^2/2;
    c(i+3,3*i+2) = -i*p;
    c(i+3,3*i+3) = i^2/2;
    c(i+3,3*i+4) = -p^2;
    c(i+3,3*i+5) = (2*i+3)*p;
    c(i+3,3*i+6) = -(2*i^2+6*i+3)/2;
    c(i+3,3*i+7) = p^2/2;
    c(i+3,3*i+8) = -(i+3)*p;
    c(i+3,3*i+9) = ((i+3)^2)/2;
end

Basis = zeros(m+2,3*m);

Basis(1,1) = m^2;
Basis(1,2) = -2*m;
Basis(1,3) = 1;
Basis(2,1) = -(3/2)*m^2;
Basis(2,2) = 2*m;
Basis(2,3) = 0;
Basis(2,4) = m^2/2;
Basis(2,5) = -2*m;
Basis(2,6) = 2;
Basis(m+2,3*m-2) = m^2;
Basis(m+2,3*m-1) = -2*m*(m-1);
Basis(m+2,3*m) = (m-1)^2;
Basis(m+1,3*m-5) = m^2/2;
Basis(m+1,3*m-4) = -m*(m-2);
Basis(m+1,3*m-3) = ((m-2)^2)/2;
Basis(m+1,3*m-2) = -(3/2)*m^2;
Basis(m+1,3*m-1) = m*(3*m-2);
Basis(m+1,3*m) = -(3*m^2-4*m)/2;

for i=0:(m-3)
    Basis(i+3,3*i+1) = m^2/2;
    Basis(i+3,3*i+2) = -i*m;
    Basis(i+3,3*i+3) = i^2/2;
    Basis(i+3,3*i+4) = -m^2;
    Basis(i+3,3*i+5) = (2*i+3)*m;
    Basis(i+3,3*i+6) = -(2*i^2+6*i+3)/2;
    Basis(i+3,3*i+7) = m^2/2;
    Basis(i+3,3*i+8) = -(i+3)*m;
    Basis(i+3,3*i+9) = ((i+3)^2)/2;
end

linear_basis = zeros(m+1,2*m);

linear_basis(1,1) = -m;
linear_basis(1,2) = 1;
linear_basis(2,1) = m;
linear_basis(2,2) = 0;
linear_basis(2,3) = -m;
linear_basis(2,4) = 2;
for i=1:(m-2)
    linear_basis(i+2,2*i+1) = m;
    linear_basis(i+2,2*i+2) = -i;
    linear_basis(i+2,2*i+3) = -m;
    linear_basis(i+2,2*i+4) = i+2;
end
linear_basis(m+1,2*m-1) = m;
linear_basis(m+1,2*m) = -(m-1);


starting_point = (1/(m+1))*ones(m+1,p+2);


theta = starting_point;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



M = m+1;
n = p+2;

like_corresponding_run = zeros(max_runs,1);
S = 'Warming up ... MCMC starting soon';
disp(S)
for ii = 1:max_runs
    old_theta = theta;
    epsilon = 1;
    if(ii == 1)
        epsilon_decreasing_factor = epsilon_decreasing_factor_1;
    else
        epsilon_decreasing_factor = epsilon_decreasing_factor_2;
    end
    
    array_of_values = zeros(no_iter,1);
    
    for i=1:no_iter
        count_1 = 0;
        count_2 = 0;
        check_each_best = zeros(M,n);
        corresponding_likelihood = zeros(n,1);
        
        current_lh_1 = -log_like_dist(Basis,linear_basis, theta, x_array, y_array);
        
        for j = 1:n
            if(min(ge(theta(:,j),0)) == 0)
                stop('error')
            end
            
            total_lh_pos = zeros(1,M);
            total_lh_neg = zeros(1,M);
            
            matrix_pos_update_at_h = zeros(M,M);
            matrix_neg_update_at_h = zeros(M,M);
            
            % pos
            for positive_change_loc = 1:M
                possibility_pos = theta(:,j);
                temp_possibility_pos = theta(:,j);
                temp_possibility_pos(positive_change_loc) = 0; % To find all significant positions except the positive_change_loc-th
                significant_positions = find(gt(temp_possibility_pos, parameter_cut_off*ones(M,1)));
                if(isempty(significant_positions) == 1)
                    possibility_pos = theta(:,j);
                else
                    possibility_pos(positive_change_loc) = theta(positive_change_loc,j) + epsilon;
                    possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon/size(significant_positions,1);
                    epsilon_temp_pos = epsilon;
                    
                    if(min(ge(possibility_pos,0)) == 0 && epsilon_temp_pos > epsilon_cut_off)
                        epsilon_temp_pos = epsilon_temp_pos/epsilon_decreasing_factor;
                        possibility_pos = theta(:,j);
                        possibility_pos(positive_change_loc) = theta(positive_change_loc,j) + epsilon_temp_pos;
                        possibility_pos(significant_positions) = possibility_pos(significant_positions) - epsilon_temp_pos/size(significant_positions,1);
                    end
                end
                
                if(min(ge(possibility_pos,0)) == 0 || isequal(possibility_pos,theta(:,j)) == 1)
                    possibility_pos = theta(:,j);
                    total_lh_pos(positive_change_loc) = current_lh_1;
                else
                    proxy_coef_theta_pos = theta;
                    proxy_coef_theta_pos(:,j) = possibility_pos;
                    total_lh_pos(positive_change_loc) = -log_like_dist(Basis,linear_basis,proxy_coef_theta_pos, x_array, y_array);
                end
                
                matrix_pos_update_at_h(:,positive_change_loc) = possibility_pos;
            end
            
            % neg
            for negative_change_loc = 1:M
                possibility_neg = theta(:,j);
                temp_possibility_neg = theta(:,j);
                temp_possibility_neg(negative_change_loc) = 0;
                significant_positions = find(gt(temp_possibility_neg, parameter_cut_off*ones(M,1)));
                if(isempty(significant_positions) == 1)
                    possibility_neg = theta(:,j);
                else
                    possibility_neg(negative_change_loc) = theta(negative_change_loc,j) - epsilon;
                    possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon/size(significant_positions,1);
                    epsilon_temp_neg = epsilon;
                    
                    if(min(ge(possibility_neg,0)) == 0 && epsilon_temp_neg > epsilon_cut_off)
                        epsilon_temp_neg = epsilon_temp_neg/epsilon_decreasing_factor;
                        possibility_neg = theta(:,j);
                        possibility_neg(negative_change_loc) = theta(negative_change_loc,j) - epsilon_temp_neg;
                        possibility_neg(significant_positions) = possibility_neg(significant_positions) + epsilon_temp_neg/size(significant_positions,1);
                    end
                end
                
                if(min(ge(possibility_neg,0)) == 0 || isequal(possibility_neg,theta(:,j)) == 1)
                    possibility_neg = theta(:,j);
                    total_lh_neg(negative_change_loc) = current_lh_1;
                else
                    proxy_coef_theta_neg = theta;
                    proxy_coef_theta_neg(:,j) = possibility_neg;
                    total_lh_neg(negative_change_loc) = -log_like_dist(Basis,linear_basis,proxy_coef_theta_neg, x_array, y_array);
                end
                matrix_neg_update_at_h(:,negative_change_loc) = possibility_neg;
            end
            
            
            % general
            
            [M_pos_1,I_pos_1] = min(total_lh_pos);
            [M_neg_1,I_neg_1] = min(total_lh_neg);
            
            check_each_best(:,j) = theta(:,j);
            corresponding_likelihood(j) = current_lh_1;
            
            if(min(M_pos_1,M_neg_1) < current_lh_1)
                count_1 = count_1+1;
                if(M_pos_1 < M_neg_1)
                    check_each_best(:,j) = matrix_pos_update_at_h(:,I_pos_1);
                    corresponding_likelihood(j) = M_pos_1;
                else
                    check_each_best(:,j) = matrix_neg_update_at_h(:,I_neg_1);
                    corresponding_likelihood(j) = M_neg_1;
                end
            end
            
        end
        
        
        [num, idx] = min(corresponding_likelihood(:));
        
        
        % Best move selection + sparsity control
        
        parameter_1 = check_each_best(:,idx);
        sparsity_positions = lt(parameter_1,parameter_cut_off*ones(M,1));
        garbage = sum(parameter_1(sparsity_positions));
        if(garbage > 0)
            parameter_1(sparsity_positions) = 0;
            rich_positions = ge(parameter_1,parameter_cut_off*ones(M,1));
            parameter_1(rich_positions) = parameter_1(rich_positions)+garbage/nnz(rich_positions);
        end
        theta(:,idx) = parameter_1;
        
        
        array_of_values(i) = current_lh_1;
        
        % Epsilon change
        if(i > 1)
            if(abs(array_of_values(i) - array_of_values(i-1)) < tol_fun)
                if(epsilon > epsilon_decreasing_factor*epsilon_cut_off)
                    epsilon = epsilon/epsilon_decreasing_factor;
                else
                    break
                end
            end
        end
        
    end
    
    like_corresponding_run(ii)= current_lh_1;
    if(ii >= 2)
        if(abs(like_corresponding_run(ii) - like_corresponding_run(ii-1)) < tol_fun_2)
            break
        end
    end
end

starting_point = theta;

%%%%% END of GCDMSVS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCMC loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = starting_point;

start_log_lh = log_like_dist(Basis,linear_basis, theta,x_array,y_array);
coef_diff_old = theta; 
coef_diff = theta; 


all_coef_diff = cell(1,mcmc_no);
total_lh_old = log_like_dist(Basis,linear_basis, theta, x_array, y_array);
acceptance_no = 1;
lh_array = zeros(mcmc_no,1);

no_iter = 0;

for i=2:mcmc_no
    if(round(i/100) == i/100)
        mcmc_iteration_no = i;
        mcmc_iteration_no        
    end
    
    for j = 1:(p+2)
        no_iter = no_iter+1;
        
        for h = 1:(m+1)
            coef_diff(h,j) = max(10^(-6),(coef_diff_old(h,j)/r + rand(1)*(r * coef_diff_old(h,j) - coef_diff_old(h,j)/r)));
        end
        coef_diff(:,j) = coef_diff(:,j)/sum(coef_diff(:,j));
        
        total_lh = log_like_dist(Basis,linear_basis, coef_diff,x_array,y_array);
        
        
        first_ratio = coef_diff_old(:,j)./coef_diff(:,j);
        min_1 = min(first_ratio);
        max_1 = max(first_ratio);
        quan_1 = log((r*min_1)^(m+1)-(max_1/r)^(m + 1))-sum(log(coef_diff_old(:,j)));
        first_ratio_inv = 1./first_ratio;
        inv_min_1 = min(first_ratio_inv);
        inv_max_1 = max(first_ratio_inv);
        quan_2 = log((r*inv_min_1)^(m+1)-(inv_max_1 /r)^(m+1))-sum(log(coef_diff(:,j)));
        
        
        transition = total_lh-quan_1+quan_2-total_lh_old;
        acceptance_prob = min(1,exp(transition));
        if(acceptance_prob == 1)
            acceptance_no = acceptance_no + 1;
            coef_diff_old = coef_diff;
            total_lh_old = total_lh;
        else
            coin = rand(1);
            if(coin <= acceptance_prob)
                acceptance_no = acceptance_no + 1;
                coef_diff_old = coef_diff;
                total_lh_old = total_lh;
            end
        end
    end
    
    all_coef_diff{i} = coef_diff;
    lh_array(i) = total_lh_old;
    
    r_deciding_factor = acceptance_no/(no_iter);
    if(r_deciding_factor < 0.15)
        aaa = r-1;
        aaa = aaa/2;
        r = 1 + aaa;
    end
    if(r_deciding_factor > 0.45)
        aaa = r-1;
        aaa = aaa*2;
        r = 1 + aaa;
    end
end

TOTAL_coef_diff = zeros(m+1,p+2);

for i = mean_start:mcmc_no
    TOTAL_coef_diff = TOTAL_coef_diff + all_coef_diff{i};
end
MEAN_coef_diff = (1/(mcmc_no-mean_start+1))*TOTAL_coef_diff;


acceptance_ratio = acceptance_no/(no_iter);


theta_mats = MEAN_coef_diff;

%%%%%%%%%%%%%%%%%% END of MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PLOTS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taus = linspace(0.05,0.95,19);


x_gridding = linspace(0, 1, 100);


DIST_tau_x_grid_mat = zeros(19, 100);

for k = 1:100
    DIST_tau_x_grid_mat(:,k) = DIST_marginal_tau_given_x(Basis, theta_mats, taus, x_gridding(k), 101);
end



figure
hold on
col=winter(19);
h = zeros(19,1);
for i = 1:19
    h(i)=plot(x_gridding,DIST_tau_x_grid_mat(i,:),'color',col(i,:));
end
scatter(data(:,1),data(:,2),'filled')
axis([0 1 0 1])
legend([h(19) h(18) h(17) h(16) h(15) h(14) h(13) h(12) h(11) h(10) h(9) h(8) h(7) h(6) h(5) h(4) h(3) h(2) h(1)],{'\tau=0.95', ...
    '\tau=0.90','\tau=0.85','\tau=0.80','\tau=0.75','\tau=0.70','\tau=0.65','\tau=0.60','\tau=0.55','\tau=0.50',...
    '\tau=0.45','\tau=0.40','\tau=0.35','\tau=0.30','\tau=0.25','\tau=0.20','\tau=0.15','\tau=0.10','\tau=0.05'});
title(['Simultaneous quantiles (NPDFSQR Complete Data)'])
xlabel('X')
ylabel('Y')
hold off


acceptance_ratio
