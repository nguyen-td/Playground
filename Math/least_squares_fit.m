%% try out anonymous function
m = 2;
n = 6;
fun_linear = @(x) m*x + n;

x = linspace(0, 10, 10);
y = fun_linear(x);
plot(y)

%% linear least squares problem
% Find the shortest distance from the origin [0 0 0] to the plane x_1 +
% 2*x_2 + 4*x_3 = 7. In other words, minimize f(x) = x_1^2 + x_2^2 + x_3^3
% subject to the constraint x_1 + 2*x_2 + 4*x_3 = 7.

pointtoplane = optimproblem;
x = optimvar('x',3);

obj = sum(x.^2); % objective function
pointtoplane.Objective = obj;

v = [1,2,4]; % linear constraint
pointtoplane.Constraints = dot(x,v) == 7;
show(pointtoplane)

[sol,fval,exitflag,output] = solve(pointtoplane);
disp(sol.x)

%% nonlinear least squares problem
% generate date from the following exponential decay model plus noise: 
% y = exp(-1.3t) + eps. Problem: given data (t,y) find the exponential
% decay rate that best fits the data.

rng default
t = linspace(0,3);
y = exp(-1.3*t) + 0.05*randn(size(t));
fun_nonlinear = @(r) exp(-r*t) - y; % y_hat-y, function to minimize

x0 = 4; % arbitrarily choose initial guess 
x = lsqnonlin(fun_nonlinear,x0);

plot(t,y,'ko')
hold on;
plot(t,exp(-x*t),'b-')
legend('Data','Best fit')
xlabel('t')
ylabel('exp(-tx)')