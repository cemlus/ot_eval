# ot_eval

```matlab
%% Simplex

% Max Z = 4x1 + 3x2
% x1 + x2 <= 8
% 2x1 + x2 <= 10

cost = [4 3 0 0];
A = [1 1 1 0;
    2 1 0 1];
b = [8 ; 10];

[m,n] = size(A);

bv_index = [3 4];

Y = [A b];

for s = 1:50
    cb = cost(bv_index);
    xb = Y(:, end);
    zj = cb * Y(:, 1:n);
    zjcj = zj - cost;

    if all(zjcj >= 0)
        sol = zeros(1,n);
        sol(bv_index) = xb;
        disp("optimal solutions: ")
        disp(sol(1:2));
        Z = cb * Y(:, end);
        disp("optimal value: ")
        disp(Z);
    end

    [~, EV] = min(zjcj);
    ratio = inf(1:m);
    for j = 1:m
        if Y(j,EV) ~= 0
            ratio(j) = xb(j) / Y(j,EV);
        end
    end

    [~, LV] = min(ratio);
    pivot = Y(LV, EV);
    bv_index(LV) = EV;

    Y(LV, :) = Y(LV, :) / pivot;

    for j = 1:m
        if j ~= LV
            Y(j, :) = Y(j, :) - Y(j, EV) * Y(LV, :);
        end
    end
end
```

```matlab
%% Big M

clc;
clear;

M = 1000;

% Convert minimization → maximization
c = [-2 -8 0 0 -M -M];   % 🔥 only change

b = [150 ; 20 ; 14];

A = [5 10 0 0 1 0;
    1 0 1 0 0 0;
    0 1 0 -1 0 1];

[m,n] = size(A);

% Initial basis (artificial + slack)
bv_index = [5 3 6];

Y = [A b];

for s = 1:50

    % -------- SAME AS SIMPLEX --------
    cb = c(bv_index);
    cb = cb(:)';

    xb = Y(:,end);

    zj = cb * Y(:,1:n);      % 🔥 fixed
    zjcj = zj - c;

    Table = [zjcj 0; Y];
    disp("Table:");
    disp(Table);

    % -------- SAME OPTIMALITY --------
    if all(zjcj >= 0)
        solution = zeros(1,n);
        solution(bv_index) = xb;

        disp("Optimal solution (x1,x2):");
        disp(solution(1:2));

        Z = cb * xb;
        disp("Optimal value:");
        disp(Z);

        break;
    end

    % -------- SAME ENTERING --------
    [~, EV] = min(zjcj);

    ratio = inf(m,1);

    for j = 1:m
        if Y(j,EV) > 0
            ratio(j) = xb(j) / Y(j,EV);
        end
    end

    % -------- SAME LEAVING --------
    [~, LV] = min(ratio);

    bv_index(LV) = EV;

    % -------- SAME PIVOT --------
    pivot = Y(LV,EV);
    Y(LV,:) = Y(LV,:) / pivot;

    for j = 1:m
        if j ~= LV
            Y(j,:) = Y(j,:) - Y(j,EV)*Y(LV,:);
        end
    end
end
```

```matlab
%% Minimum Cost Method

clc;
clear;

% Cost matrix
cost = [11 13 17 14;
    16 18 14 10;
    21 24 13 10];

% Supply and Demand
supply = [10 5 9];
demand = [8 7 15 4];

[m,n] = size(cost);

S = sum(supply);
D = sum(demand);

% -------- Balance --------
if (S == D)
    disp("Problem is balanced")
elseif (S < D)
    cost(end+1,:) = zeros(1,n);
    supply(end+1) = D - S;
else
    cost(:,end+1) = zeros(m,1);
    demand(end+1) = S - D;
end

disp("Balanced Transportation Problem")

[m,n] = size(cost);

X = zeros(m,n);
Icost = cost;

% -------- Least Cost Method --------
while (any(supply ~= 0) && any(demand ~= 0))

    min_cost = min(cost(:));

    [r, c] = find(cost == min_cost, 1);

    aloc = min(supply(r), demand(c));

    X(r,c) = X(r,c) + aloc;

    supply(r) = supply(r) - aloc;
    demand(c) = demand(c) - aloc;

    % Block row/column
    if supply(r) == 0
        cost(r,:) = inf;
    end

    if demand(c) == 0
        cost(:,c) = inf;
    end
end

% -------- Total Cost --------
cost_ec = X .* Icost;
final_cost = sum(cost_ec(:));

disp("Allocation Matrix:")
disp(X)

disp("Total Cost:")
disp(final_cost)
```

```matlab
% STEEPEST DESCENT METHOD
% QUESTION IS - PERFORM 4 ITERATIONS OF STEEPEST DESCENT ALGORITHM TO MINIMIZE
% f(x1,x2) = x1 - x2 + 2*x1^2 + 2*x1*x2 + x2^2

% Starting from the point X1 = [1,1]

format short         
clear all             
clc                   

% Phase 1- Define objective function
syms x1 x2

f1 = x1 - x2 + 2*x1^2 + 2*x1*x2 + x2^2;

fx = inline(f1);   
fobj = @(x) fx(x(:,1), x(:,2));   

% Phase 2- Compute gradient of f
grad = gradient(f1);
G = inline(grad);  
gradx = @(x) G(x(:,1), x(:,2));

% Phase 3- Compute Hessian Matrix
H1 = hessian(f1);
Hx = inline(H1);

% Phase 4- Iterations
x0 = [1 1];         
maxiter = 4;        
tol = 10^(-3);     

iter = 0;           
X = [];             

while norm(gradx(x0)) > tol && iter < maxiter
    X = [X; x0];           

    S = -gradx(x0);         
    H = Hx(x0);             

    lambda = (S'*S) / (S'*H*S);  

    Xnew = x0 + lambda.*S'; 
    x0 = Xnew;             

    iter = iter + 1;
end

% Phase 5- Print the solution
fprintf('Optimal Solution x=[%f, %f]\n', x0(1), x0(2))
fprintf('Optimal Value f(x) = %f\n', fobj(x0))

```
