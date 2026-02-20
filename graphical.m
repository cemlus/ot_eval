clc;
clear all;
close all;

% minimize Z = 3x1 - 5x2
% x1 + x2 <= 6
% 2*x1 - x2 <= 9
% x1 >= 0, x2 >= 0

% phase 1 => input parameters
C = [3, -5];
A = [1 2; 1 -1];
B = [6 ; 9]  % ensure this is always a column matrix

z = @(x1, x2) 3*x1 - 5*x2;
c1 = @(x1,x2) x1 + x2 - 6;
c2 = @(x1,x2) 2*x1 - x2 - 9;

[m,n] = size(A);


% phase 2 => plotting constraints
x1 = 0:max(B(B>0)./A(B>0,1))

for i= 1:m
	x2 = (B(i) - A(i,1)*x1)/A(i,2);
	plot(x1, x2);
	hold on;
	
end

% phase 3 => find intersection and corner points

A = [A; eye(2)];
B = [B; 0; 0];
m = size(A,1)
pt = []
for i=1:m
	for j=i+1:m
		aa = [A(i, :) ; A(j, :)];
		bb = [B(i) ; B(j)];
		d = det(aa)
		
		if(d ~= 0)
			X = aa\bb;
			if(X>=0)
				pt = [pt X];
			end
		end
	end
end


% phase 4 => find feasible points
FP = []
Z = []

for i = 1:size(pt,2)
	PT1 = pt(1, i);
	PT2 = pt(2, i);
	
	if(c1(PT1,PT2)<=0 && c2(PT1,PT2)<=0)
		FP = [FP pt(:,i)]
		plot(PT1,PT2, 'r*', 'MarkerSize', 10)
		cost = z(PT1,PT2)
		Z = [Z cost]	
	end
end

% phase 5 => find optimal solution
[optimal_val, idx] = max(Z)
optimal_sol = FP(:,idx)

printf('the optimal solution is: %.2f\n', optimal_val)
printf('the optimal points are x1: %.2f and x2: %.2f\n', optimal_sol(1), optimal_sol(2))
