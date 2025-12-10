%% Example of solution of a linear program in Matlab
%   minimize   c' * x
%   subject to A * x <= b
%              lb <= x <= ub
%
% How to used it:
%   - Edit the example (c, A, b, lb, ub) below and run.
%   - Works for 2 decision variables only (x = [x1; x2]).


clear; 
clc; 
close all;

% minimize (c' * [x1; x2]) 
c  = [-3; 2];

A  = [1 1;      
      2 1;      
      1 0];     
b  = [4; 5; 3];

% Example of Ax <= b: 
% A  = [1 1;      % x1 + x2 <= 4
%       2 1;      % 2*x1 + x2 <= 5
%       1 0];     % x1 <= 3
% b  = [4; 5; 3];

lb = [0; 0];              % lower bounds (use -Inf for no lower bound)
ub = [Inf; Inf];          % upper bounds (use  Inf for no upper bound)

 
tol = 1e-9;

% -- Solve LP
opts = optimoptions('linprog','Display','none');
[x_opt, fval, exitflag] = linprog(c, A, b, [], [], lb, ub, opts);

% -- Report
if exitflag ~= 1
    fprintf('linprog exitflag = %d (solution may be infeasible/unbounded/nonoptimal)\n', exitflag);
else
    fprintf('Optimal solution found.\n');
end
fprintf('x* = [x1; x2] = [%.6g; %.6g]\n', x_opt(1), x_opt(2));
fprintf('Min objective value c''x* = %.6g\n', fval);


%% Function to plot the region and the optimal point

% ----------------- Build boundary lines -----------------
% Each inequality a^T x <= b becomes boundary a^T x = b for intersections.
lines = {}; rhs = [];
if ~isempty(A)
    for i = 1:size(A,1)
        lines{end+1} = A(i,:);
        rhs(end+1,1) = b(i);  
    end
end
% Add finite bound lines: x1 = lb1 / ub1, x2 = lb2 / ub2
% Represent x1 = const as [1 0] x = const; x2 = const as [0 1] x = const.
if isfinite(lb(1)), lines{end+1} = [1 0]; rhs(end+1,1) = lb(1); end
if isfinite(ub(1)), lines{end+1} = [1 0]; rhs(end+1,1) = ub(1); end
if isfinite(lb(2)), lines{end+1} = [0 1]; rhs(end+1,1) = lb(2); end
if isfinite(ub(2)), lines{end+1} = [0 1]; rhs(end+1,1) = ub(2); end

% ----------------- Candidate intersections -----------------
P = [];  % columns are points
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

% Also include potential bound corners if both lb/ub finite
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

    % Deduplicate and order
    if ~isempty(Pf)
        Pf = unique(round(Pf.', 10), 'rows').';
        ctr = mean(Pf, 2);
        ang = atan2(Pf(2,:) - ctr(2), Pf(1,:) - ctr(1));
        [~, ord] = sort(ang);
        Pf = Pf(:, ord);
    end
end

% ----------------- Figure & plotting window -----------------
figure; hold on; grid on; box on;
xlabel('x_1'); ylabel('x_2');

% Determine view limits from feasible points and x_opt
pts = [Pf, x_opt];
if isempty(pts)
    xlim_ = [-1 1]; ylim_ = [-1 1];
else
    xmin = min(pts(1,:)); xmax = max(pts(1,:));
    ymin = min(pts(2,:)); ymax = max(pts(2,:));
    pad = 0.1 * max(1, max([xmax - xmin, ymax - ymin]));
    xlim_ = [xmin - pad, xmax + pad];
    ylim_ = [ymin - pad, ymax + pad];
end

% ----------------- Draw feasible polygon (if it exists) -----------------
if size(Pf,2) >= 3
    fill(Pf(1,:), Pf(2,:), 0.92*[1 1 1], 'EdgeColor', 0.5*[1 1 1], 'LineWidth', 1.5, ...
         'DisplayName','Feasible region');
else
    % No polygon (unbounded or degenerate) — just note it
    text(mean(xlim_), mean(ylim_), 'Feasible region unbounded or degenerate', ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.4 0.4 0.4]);
end

% ----------------- Plot constraint/boundary lines (within window) -----------------
% Draw each boundary a^T x = b_i (and finite bound lines) where possible.
xx = linspace(xlim_(1), xlim_(2), 400);
for i = 1:numel(lines)
    a = lines{i}; bi = rhs(i);
    % Try to express as x2 = (bi - a1*x1)/a2 if |a2|>0, else vertical line x1 = bi/a1
    if abs(a(2)) > 1e-12
        yy = (bi - a(1)*xx)/a(2);
        plot(xx, yy, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
    elseif abs(a(1)) > 1e-12
        xconst = bi / a(1);
        plot([xconst xconst], ylim_, 'k-', 'LineWidth', 0.8, 'HandleVisibility','off');
    end
end

% ----------------- Plot optimal point -----------------
if isfinite(x_opt(1)) && isfinite(x_opt(2))
    plot(x_opt(1), x_opt(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName','Optimal x*');
    text(x_opt(1)+0.02*(xlim_(2)-xlim_(1)), x_opt(2), ...
        sprintf('x* = (%.3g, %.3g)', x_opt(1), x_opt(2)), 'FontWeight','bold');
end

% ----------------- Objective level set through optimum -----------------
% Plot {x : c' x = c' x_opt} as a reference.
const = c.' * x_opt;
if abs(c(2)) > 1e-12
    yy_obj = (const - c(1)*xx) / c(2);
    plot(xx, yy_obj, 'r-', 'LineWidth', 1.2, 'DisplayName','Objective level (c''x = const)');
elseif abs(c(1)) > 1e-12
    xvert = const / c(1);
    plot([xvert xvert], ylim_, 'r-', 'LineWidth', 1.2, 'DisplayName','Objective level (c''x = const)');
end

% ----------------- Finalize axes -----------------
xlim(xlim_); ylim(ylim_);
axis tight;
legend('Location','northeastoutside');

% ----------------- Print which inequalities are tight at optimum -----------------
if ~isempty(A)
    res = b - A*x_opt;
    tight = find(abs(res) < 1e-7);
    if ~isempty(tight)
        fprintf('Active (tight) A*x<=b rows at optimum: %s\n', mat2str(tight));
    end
end
% Bound activity
if isfinite(lb(1)) && abs(x_opt(1) - lb(1)) < 1e-7, fprintf('x1 at lower bound.\n'); end
if isfinite(ub(1)) && abs(x_opt(1) - ub(1)) < 1e-7, fprintf('x1 at upper bound.\n'); end
if isfinite(lb(2)) && abs(x_opt(2) - lb(2)) < 1e-7, fprintf('x2 at lower bound.\n'); end
if isfinite(ub(2)) && abs(x_opt(2) - ub(2)) < 1e-7, fprintf('x2 at upper bound.\n'); end

