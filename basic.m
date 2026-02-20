%ques 2 for algebraic method

% algebraic method-
% 1. basic solution set
% 2. basic feasible solution set
% 3. degenerate & non-degenerate basic feasible solutions
% 4. optimal solution
% 5. optimal value

% Max Z = x1 + 2x2
% subject to:
% -x1 + x2 <= 1
%  x1 + x2 <= 1
% x1, x2 >= 0

% standard form:
% Max Z = x1 + 2x2 + 0x3 + 0x4
% subject to:
% -x1 + x2 + x3 = 1
%  x1 + x2 + x4 = 2

% X is column vector [x1; x2; x3; x4]

clc
clear all
clear figure

%  Phase 1 : Input Parameters ----------------

c = [1 2 0 0];              % cost coefficients
A = [-1 1 1 0;              % constraint matrix
      1 1 0 1];
b = [1; 1];                 % RHS column vector

z = @(X) c*X;               % objective function
m = size(A,1);              % number of constraints
n = size(A,2);              % number of variables in standard form

% Phase 2 : Basic & Basic Feasible Solutions ----------------

basicsol = [];              % stores all basic solutions
bfsol = [];                 % stores all basic feasible solutions
degenerate = [];            % stores degenerate BFS
nondegenerate = [];         % stores non-degenerate BFS

ncm = nchoosek(n,m);        % total number of combinations
pair = nchoosek(1:n,m);     % all possible selections of basic variables

for i = 1:ncm
    
    y = zeros(n,1);         
    % starting solution vector with all variables zero
    
    basicvar_index = pair(i,:);  
    % selecting one combination of basic variables
    
    B = A(:, basicvar_index);
    if abs(det(B)) < 1e-6
          continue;
    end
    X = B\b;
    
    % X = A(:,basicvar_index)\b;    
    % solving AX = b for chosen basic variables
    
    y(basicvar_index) = X;        
    % forming full basic solution vector
    
    basicsol = [basicsol y];      
    % storing basic solution
    
    if all(X >= 0)                 
        % checking basic feasible solution condition
        
        bfsol = [bfsol y];
        
        if any(X == 0)
            % degenerate BFS (as per modified condition)
            degenerate = [degenerate y];
        else
            % degenerate BFS
            nondegenerate = [nondegenerate y];

        end
        
    end
end

disp('Basic Solution Set:');
disp(basicsol);

disp('Basic Feasible Solution Set:');
disp(bfsol);

disp('Degenerate Basic Feasible Solutions:');
disp(degenerate);

disp('Non-Degenerate Basic Feasible Solutions:');
disp(nondegenerate);

%  Phase 3 : Optimal Solution & Optimal Value ----------------

cost = z(bfsol);             % computing Z for each BFS
[optimum_val, index] = max(cost); % MAX problem
optimum_sol = bfsol(:,index);     % optimal solution

disp('Optimal Solution:');
disp(optimum_sol);

disp('Optimal Value:');
disp(optimum_val);
