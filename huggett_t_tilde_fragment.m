% PROGRAM NAME: ps4huggett.m
clear, clc

% PARAMETERS
beta = .9932; % discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix

PI_star = % invariant distribution



% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5; %upper bound of grid points - can try upper bound of 3
num_a = 10; %try 700 points

a = linspace(a_lo, a_hi, num_a); % asset (row) vector

% INITIAL GUESS FOR q
q_min = 0.98;
q_max = 1;
q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
aggsav = 1 ;
while abs(aggsav) >= 0.01 ;
    
    % CURRENT RETURN (UTILITY) FUNCTION
    cons = bsxfun(@minus, a', q_guess * a); % where cons is 3 dim - a row vector, a', subtract q*a
    cons = bsxfun(@plus, cons, permute(y_s, [1 3 2])); % permute - rearranges - adding the third dimension
    ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
    ret (cons<0) = -Inf;
    
    % INITIAL VALUE FUNCTION GUESS
    v_guess = zeros(2, num_a); % 2xN 
    
    % VALUE FUNCTION ITERATION
    v_tol = 1;
    while v_tol >.0001;
    
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
        
        value_mat = ret + beta * ...
            repmat(permute((PI * v_guess), [3 2 1]), [num_a 1 1]); % multiplying PI*v_guess - getting expectation
        
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        
        [vfn, pol_indx] = max(value_mat, [], 2);
        vfn = permute(vfn,[3 1 2]);
        
        v_tol = abs(max(v_guess(:) - vfn(:)));
        
        v_guess = vfn; % update value functions
        
  
    end;
    
    % KEEP DECSISION RULE
    pol_indx=permute(pol_indx,[3 1 2]);
    pol_fn = a(pol_indx);
    
    % SET UP INITITAL DISTRIBUTION
    Mu=zeros(2,num_a); %any initial distribution works, as long as they sum up to 1 - can be uniform dist or put all mass at one point
    Mu(1,4) = 1; %suppose full mass at one point
    
    %function distribution = (pol_fn, PI);
    
    % ITERATE OVER DISTRIBUTIONS
    
    
    mu_tol = 1;
    
    while mu_tol> 1e-08
        [emp_ind, a_ind, mass] = find(Mu > 0); % only looping over nonzero indices- find non-zero indices - employment and asset index
        
                  
        MuNew = zeros(size(Mu));
        for ii = 1:length(emp_ind)
            apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
            MuNew(:, apr_ind) = MuNew(:, apr_ind) + ... % which mass of households goes to which exogenous state?
                (PI(emp_ind(ii), :) * mass)';
        end
    
    mu_tol = max(abs(MuNew(:) - Mu(:)));
    
    Mu = MuNew;
    
    plot(Mu') % look at the distribution
    plot(Mu(2,:)') % look at the unemployed distribution
    sum(Mu(:)) % check if it sums to 1
    
    % AGGREGATE/ INTEGRATE AND CHECK FOR MARKET CLEARING
    
    Mu.*pol_fn % multiply MU * pol-fn, tells us how much total saving of people at any given state
    
    sum(sum(Mu.*pol_fn)) %to get agg saving, sum it up; check if it is close to 0; so now adjust bond price, and repeat until get close to 0
       
    
        
end
