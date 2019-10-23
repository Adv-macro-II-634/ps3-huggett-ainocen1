% PROGRAM NAME: ps4huggett.m
clear, clc

% PARAMETERS
beta = .9932; % discount factor 
sigma = 1.5; % coefficient of risk aversion
b = 0.5; % replacement ratio (unemployment benefits)
y_s = [1, b]; % endowment in employment states
PI = [.97 .03; .5 .5]; % transition matrix

% ASSET VECTOR
a_lo = -2; %lower bound of grid points
a_hi = 5; %upper bound of grid points - can try upper bound of 3
num_a = 900; %try 700 points

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
    while v_tol >.000001;
    
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
                (PI(emp_ind(ii), :) * Mu(emp_ind(ii),a_ind(ii)))';
        end
    
        mu_tol = max(abs(MuNew(:) - Mu(:)));

        Mu = MuNew;

    end
    
    plot(Mu') % look at the distribution
    plot(Mu(2,:)') % look at the unemployed distribution
    sum(Mu(:)) % check if it sums to 1
    
    
    
    % AGGREGATE/ INTEGRATE AND CHECK FOR MARKET CLEARING
    
    % Mu.*pol_fn; % multiply MU * pol-fn, tells us how much total saving of people at any given state
    aggsav = sum(sum(Mu.*pol_fn)); %to get agg saving, sum it up; check if it is close to 0; so now adjust bond price, and repeat until get close to 0
    if aggsav>0 %suppose saving too much, q might be too low
        q_min=(q_min + q_max) / 2;
    else
        q_max=(q_min + q_max) / 2;
    end
     
    q_guess = (q_min + q_max) / 2;
    
        
end

% resulting risk-free interest rate (q) = 0.9945


%Lorenz curve and Gini for earnings and wealth

% define wealth = asset plus income
w1 = [a+1 a+0.5].*Mu(:)';%wealth
wealth = w1(w1>=0); 
%earnings=[y_s(1) y_s(2)].*Mu(:)';
popn = Mu(:)';

G=gini(popn,wealth,true);
figure
%G2=gini(popn,earnings,true);
figure



%Gini fcn
function [g,l,a] = gini(pop,val,makeplot)
    % check arguments
    assert(nargin >= 2, 'gini expects at least two arguments.')
    if nargin < 3
        makeplot = false;
    end
    assert(numel(pop) == numel(val), ...
        'gini expects two equally long vectors (%d ~= %d).', ...
        size(pop,1),size(val,1))
    pop = [0;pop(:)]; val = [0;val(:)];     % pre-append a zero
    isok = all(~isnan([pop,val]'))';        % filter out NaNs
    if sum(isok) < 2
        warning('gini:lacking_data','not enough data');
        g = NaN; l = NaN(1,4);
        return;
    end
    pop = pop(isok); val = val(isok);
    
    assert(all(pop>=0) && all(val>=0), ...
        'gini expects nonnegative vectors (neg elements in pop = %d, in val = %d).', ...
        sum(pop<0),sum(val<0))
    
    % process input
    z = val .* pop;
    [~,ord] = sort(val);
    pop    = pop(ord);     z    = z(ord);
    pop    = cumsum(pop);  z    = cumsum(z);
    relpop = pop/pop(end); relz = z/z(end);
    
    % Gini coefficient
    g = 1 - sum((relz(1:end-1)+relz(2:end)) .* diff(relpop));
    
    % Lorenz curve
    l = [relpop,relz];
    a = [pop,z];
    if makeplot   % ... plot it?
        area(relpop,relz,'FaceColor',[0.5,0.5,1.0]);    % the Lorentz curve
        hold on
        plot([0,1],[0,1],'--k');                        % 45 degree line
        axis tight      % ranges of abscissa and ordinate are by definition exactly [0,1]
        axis square     % both axes should be equally long
        set(gca,'XTick',get(gca,'YTick'))   % ensure equal ticking
        set(gca,'Layer','top');             % grid above the shaded area
        grid on;
        title(['\bfGini coefficient = ',num2str(g)]);
        xlabel('share of population');
        ylabel('share of value');
    end
    
end
