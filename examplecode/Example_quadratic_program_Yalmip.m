%% Example of solution of an optimization problem with quadratic cost and linear constraints in Matlab
%   minimize   (1/2)*x'*H*x + c'*x
%   subject to A * x <= b
%              lb <= x <= ub
%
% Usage:
%   - Edit the example (H, c, A, b, lb, ub) below and run.
%
% Notes:
%   - x = [x1; x2].
%   - Assumes H is symmetric positive semidefinite (convex QP).

% Note: this version uses Yalmip: 
% Requirements: YALMIP installed (https://yalmip.github.io) and a QP solver (MATLAB's quadprog works).

clear; clc; close all;

%% ===== Example (edit these) =====
H  = [2 0.5; 0.5 1.5];    % symmetric PSD for convexity
c  = [-3; -2];

A  = [1 1;
      2 1;
      1 0];
b  = [4; 5; 3];

lb = [0; 0];
ub = [Inf; Inf];


tol = 1e-9;

% ---- Symmetrize H and (soft) convexity check
H = 0.5*(H + H.');
ev = eig(H);
if min(ev) < -1e-10
    warning('H is not positive semidefinite (min eigenvalue = %.3g). QP may be nonconvex.', min(ev));
end

% ---- YALMIP variables and constraints
%x1 = sdpvar(1,1); x2 = sdpvar(1,1);
x  = sdpvar(2,1);%[x1; x2];

Cons = [];
if ~isempty(A)
    Cons = [Cons; A*x <= b];
end
if isfinite(lb(1)), Cons = [Cons; x(1) >= lb(1)]; end
if isfinite(ub(1)), Cons = [Cons; x(1) <= ub(1)]; end
if isfinite(lb(2)), Cons = [Cons; x(2) >= lb(2)]; end
if isfinite(ub(2)), Cons = [Cons; x(2) <= ub(2)]; end

obj = 0.5*x.'*H*x + c.'*x;

% ---- Solve (choose a solver you have: 'quadprog', 'gurobi', 'mosek', etc.)
opts = sdpsettings('solver','quadprog','verbose',0);
diag = optimize(Cons, obj, opts);

% ---- Report
code = diag.problem; % 0=OK, 1=infeasible, 2=unbounded, 3=inaccurate, 4=numerical
if code ~= 0
    fprintf('YALMIP/solver status: problem=%d, info="%s"\n', code, diag.info);
else
    fprintf('Optimal solution found.\n');
end

x_opt = value(x);
fval  = value(obj);
fprintf('x* = [x1; x2] = [%.6g; %.6g]\n', x_opt(1), x_opt(2));
fprintf('Min objective value = (1/2)*x''H x + c''x = %.6g\n', fval);

% ----------------- Build boundary lines (for plotting) -----------------
lines = {}; rhs = [];
if ~isempty(A)
    for i = 1:size(A,1)
        lines{end+1} = A(i,:); 
        rhs(end+1,1) = b(i);  
    end
end
if isfinite(lb(1)), lines{end+1} = [1 0]; rhs(end+1,1) = lb(1); end
if isfinite(ub(1)), lines{end+1} = [1 0]; rhs(end+1,1) = ub(1); end
if isfinite(lb(2)), lines{end+1} = [0 1]; rhs(end+1,1) = lb(2); end
if isfinite(ub(2)), lines{end+1} = [0 1]; rhs(end+1,1) = ub(2); end

% ----------------- Candidate intersections -----------------
P = [];
for i = 1:numel(lines)
    ai = lines{i};  bi = rhs(i);
    for j = i+1:numel(lines)
        aj = lines{j};  bj = rhs(j);
        Aeq = [ai; aj];
        if rank(Aeq) == 2
            p = Aeq \ [bi; bj];
            P = [P, p]; %#ok<AGROW>
        end
    end
end

% Include finite bound corners (if both coordinates finite)
corners = [];
xs = [lb(1), lb(1), ub(1), ub(1)];
ys = [lb(2), ub(2), lb(2), ub(2)];
for k = 1:4
    if isfinite(xs(k)) && isfinite(ys(k))
        corners = [corners, [xs(k); ys(k)]]; %#ok<AGROW>
    end
end
P = [P, corners];

% ----------------- Keep feasible intersections -----------------
if isempty(P)
    Pf = [];
else
    feasA = true(1, size(P,2));
    if ~isempty(A)
        feasA = all(A * P <= b + tol, 1);
    end
    feasLB = ( ~isfinite(lb(1)) | P(1,:) >= lb(1) - tol ) & ...
             ( ~isfinite(lb(2)) | P(2,:) >= lb(2) - tol );
    feasUB = ( ~isfinite(ub(1)) | P(1,:) <= ub(1) + tol ) & ...
             ( ~isfinite(ub(2)) | P(2,:) <= ub(2) + tol );
    Pf = P(:, feasA & feasLB & feasUB);

    if ~isempty(Pf)
        Pf = unique(round(Pf.', 10), 'rows').';
        ctr = mean(Pf, 2);
        ang = atan2(Pf(2,:) - ctr(2), Pf(1,:) - ctr(1));
        [~, ord] = sort(ang);
        Pf = Pf(:, ord);
    end
end

% ----------------- Figure & window -----------------
figure; hold on; grid on; box on;
xlabel('x_1'); ylabel('x_2');

% Determine view limits from feasible vertices and x_opt
pts = [Pf, x_opt];
if isempty(pts)
    xlim_ = [-1 1]; ylim_ = [-1 1];
else
    xmin = min(pts(1,:)); xmax = max(pts(1,:));
    ymin = min(pts(2,:)); ymax = max(pts(2,:));
    span = max(1, max([xmax - xmin, ymax - ymin]));
    pad  = 0.15 * span;
    xlim_ = [xmin - pad, xmax + pad];
    ylim_ = [ymin - pad, ymax + pad];
end
xlim(xlim_); ylim(ylim_);

% ----------------- Feasible polygon (if exists) -----------------
if size(Pf,2) >= 3
    fill(Pf(1,:), Pf(2,:), 0.92*[1 1 1], 'EdgeColor', 0.5*[1 1 1], 'LineWidth', 1.5, ...
         'DisplayName','Feasible region');
else
    text(mean(xlim_), mean(ylim_), 'Feasible region unbounded or degenerate', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.4 0.4 0.4]);
end

% ----------------- Plot boundary lines -----------------
xx = linspace(xlim_(1), xlim_(2), 400);
for i = 1:numel(lines)
    a = lines{i}; bi = rhs(i);
    if abs(a(2)) > 1e-12
        yy = (bi - a(1)*xx)/a(2);
        plot(xx, yy, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
    elseif abs(a(1)) > 1e-12
        xconst = bi / a(1);
        plot([xconst xconst], ylim_, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
    end
end

% ----------------- Objective contours -----------------
nx = 250; ny = 250;
[X1, X2] = meshgrid(linspace(xlim_(1), xlim_(2), nx), linspace(ylim_(1), ylim_(2), ny));
Z = zeros(size(X1));
for i = 1:nx
    Xi = [X1(:,i), X2(:,i)];
    Z(:,i) = 0.5*sum((Xi*H).*Xi, 2) + Xi*c;
end
base = value(obj);
if isfinite(base)
    levels = base + linspace(-1, 5, 6)*max(1, abs(base)*0.2);
else
    levels = linspace(min(Z, [], 'all'), max(Z, [], 'all'), 8);
end
[C,hc] = contour(X1, X2, Z, levels, 'LineWidth', 1.0);
hc.DisplayName = 'Objective contours';

% ----------------- Plot optimal point -----------------
if all(isfinite(x_opt))
    plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName','Optimal x*');
    text(x_opt(1)+0.02*(xlim_(2)-xlim_(1)), x_opt(2), ...
        sprintf('x* = (%.3g, %.3g)', x_opt(1), x_opt(2)), 'FontWeight','bold');
end

legend('Location','northeastoutside');


