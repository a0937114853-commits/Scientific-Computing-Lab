









clear;clc;format long g;



f = @(x,y) 1.25*exp(x+y/2);

function p=poisson(a,b,m,f,)
	h = (b-a)/(m+1);
	x = linspace(a,b,m+2);
	y = linspace(a,b,m+2);


	[X,Y] = meshgrid(x,y);	
	X = X';                     	
	Y = Y';                     	

	Iint = 2:m+1;              	
	Jint = 2:m+1;              	
	Xint = X(Iint,Jint);	
	Yint = Y(Iint,Jint);


	rhs = f(Xint,Yint);


	utrue = exp(X+Y/2);
	usoln = utrue;



	rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
	rhs(:,m) = rhs(:,m) - usoln(Iint,m+2)/h^2;
	rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
	rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;



	F = reshape(rhs,m*m,1);

	I = speye(m);
	e = ones(m,1);
	T = spdiags([e -4*e e],[-1 0 1],m,m);
	S = spdiags([e e],[-1 1],m,m);
	A = (kron(I,T) + kron(S,I)) / h^2;

	uvec = A\F;  
	usoln(Iint,Jint) = reshape(uvec,m,m);
	err = max(max(abs(usoln-utrue)));   
	fprintf('Error relative to true solution of PDE = %10.3e \n',err)











% plot results:

clf
hold on


plot(X,Y,'g');  plot(X',Y','g')

contour(X,Y,usoln,30,'k')

axis([a b a b])
daspect([1 1 1])
title('Contour plot of computed solution')
hold off







































clear;clc;format long g;

%%
x1=0;x2=1;y1=0,y2=2;		
m_x=20;m_y=40;
h_x=(x2-x1)/(m_x+1);
h_y=(y2-y1)/(m_x+1);

%%
x=linspace(x1,x2,m_x+2);
y=linspace(y1,y2,m_y+2);
[X,Y]=meshgrid(x,y);
X=X';Y=Y';


%%
Iint = 2:m_x+1;              
Jint = 2:m_y+1;              
Xint = X(Iint,Jint);       
Yint = Y(Iint,Jint);



				
%%
f = @(x,y) 1.25*exp(x+y/2);
rhs = f(Xint,Yint);       
u_true = exp(X+Y/2);
u_soln = u_true;              
                            



%%
rhs(:,1) = rhs(:,1) - u_soln(Iint,1)/(h_y)^2;
rhs(:,m_y) = rhs(:,m_y) - u_soln(Iint,m_y+2)/(h_y)^2;
%%
rhs(1,:) = rhs(1,:) - u_soln(1,Jint)/(h_x)^2;
rhs(m_x,:) = rhs(m_x,:) - u_soln(m_x+2,Jint)/(h_x)^2;



%%
F = reshape(rhs,m_x*m_y,1);

I_x = speye(m_x);
I_y = speye(m_y);

e_x = ones(m_x,1);
T_x = spdiags([e_x -2*e_x e_x],[-1 0 1],m_x,m_x);
A_x =kron(I_y,T_x)/(h_x)^2;	


e_y = ones(m_y,1);
S_y = spdiags([e_y -2*e_y e_y],[-1 0 1],m_y,m_y);
A_y =kron(S_y,I_x)/(h_y)^2;

A=A_x+A_y;



%%
u_vec = A \ F;
u_soln(Iint, Jint) = reshape(u_vec, m_x, m_y);



surf(X,Y,u_soln)























subject to Dirichlet boundary conditions $u(x, y) = g(x, y)$ on $\partial \Omega$. The source term and the exact solution for this test problem are:
\begin{align*}
f(x, y) &= 1.25 \exp(x + y/2) \\
u_{true}(x, y) &= \exp(x + y/2)
\end{align*}



















































clear;clc;format long g;

f=@(x)x;
m=4;




A=spdiags(ones(m,1)*[1,-2,1],[-1,0,1],m,m);
D=diag(A)


for i=2:m
	for j=1:i-1
	L(i,j)=-A(i,j)
	end
end

for j=2:m
	for i=1:j-1
	U(i,j)=-A(i,j)
	end
end
























clear; clc; format long g;

% 定義測試的網格數量序列 (讓 m 每次翻倍)
m_values = [20, 40, 80, 160];
errors = zeros(size(m_values));
h_values = zeros(size(m_values));

fprintf('%5s | %10s | %10s | %10s\n', 'm', 'h', 'Max Error', 'Order');
fprintf('------------------------------------------------------------\n');

for k = 1:length(m_values)
    m = m_values(k);
    a = 0; b = 1;
    h = (b-a)/(m+1);
    h_values(k) = h;
    
    x = linspace(a,b,m+2);
    y = linspace(a,b,m+2);
    [X,Y] = meshgrid(x,y);
    X = X'; Y = Y';

    Iint = 2:m+1; Jint = 2:m+1;
    Xint = X(Iint,Jint); Yint = Y(Iint,Jint);

    f = @(x,y) 1.25*exp(x+y/2);
    rhs = f(Xint,Yint);
    utrue = exp(X+Y/2);
    usoln = utrue; % 設定邊界條件

    % 修正 RHS 邊界項
    rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
    rhs(:,m) = rhs(:,m) - usoln(Iint,m+2)/h^2;
    rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
    rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;

    % 轉換為向量並構建矩陣 A
    F = reshape(rhs, m*m, 1);
    I = speye(m);
    e = ones(m,1);
    T = spdiags([e -4*e e], [-1 0 1], m, m);
    S = spdiags([e e], [-1 1], m, m);
    A = (kron(I,T) + kron(S,I)) / h^2;

    % 求解
    uvec = A\F;
    usoln(Iint,Jint) = reshape(uvec, m, m);

    % 計算最大誤差 (L-infinity norm)
    errors(k) = max(max(abs(usoln - utrue)));

    % 計算收斂階數 (Order) = log(E1/E2) / log(h1/h2)
    if k > 1
        order = log(errors(k-1)/errors(k)) / log(h_values(k-1)/h_values(k));
        fprintf('%5d | %10.3e | %10.3e | %10.2f\n', m, h, errors(k), order);
    else
        fprintf('%5d | %10.3e | %10.3e | %10s\n', m, h, errors(k), '---');
    end
end

% 繪製 Log-Log 圖驗證
figure;
loglog(h_values, errors, '-o', 'LineWidth', 2);
grid on; hold on;
% 畫出一條斜率為 2 的參考線 (y = Cx^2)
loglog(h_values, (errors(1)/h_values(1)^2) * h_values.^2, '--r'); 
xlabel('Grid spacing h'); ylabel('Max Error');
legend('Computed Error', 'Theoretical 2nd Order');
title('Grid Refinement Study (Log-Log Plot)');




























clear; clc; format long g;


x1=0; x2=1; 
y1=0; y2=2;		
m_x=20; m_y=40;


h_x = (x2-x1)/(m_x+1);
h_y = (y2-y1)/(m_y+1);


x = linspace(x1,x2,m_x+2);
y = linspace(y1,y2,m_y+2);
[X,Y] = meshgrid(x,y);
X=X'; Y=Y';


Iint = 2:m_x+1;
Jint = 2:m_y+1;
Xint = X(Iint,Jint);
Yint = Y(Iint,Jint);


f = @(x,y) 1.25*exp(x+y/2);
rhs = f(Xint,Yint);
u_true = exp(X+Y/2);
u_soln = u_true;

%%
rhs(:,1) = rhs(:,1) - u_soln(Iint,1)/(h_y)^2;
rhs(:,m_y) = rhs(:,m_y) - u_soln(Iint,m_y+2)/(h_y)^2;
%%
rhs(1,:) = rhs(1,:) - u_soln(1,Jint)/(h_x)^2;
rhs(m_x,:) = rhs(m_x,:) - u_soln(m_x+2,Jint)/(h_x)^2;

%% 建構大型稀疏矩陣 A
F = reshape(rhs, m_x*m_y, 1);

I_x = speye(m_x);
I_y = speye(m_y);

e_x = ones(m_x,1);
T_x = spdiags([e_x -2*e_x e_x], [-1 0 1], m_x, m_x);
A_x = kron(I_y, T_x) / (h_x^2);	

e_y = ones(m_y,1);
S_y = spdiags([e_y -2*e_y e_y], [-1 0 1], m_y, m_y);
A_y = kron(S_y, I_x) / (h_y^2);

A = A_x + A_y;







%% 求解
u_vec = A \ F;
u_soln(Iint, Jint) = reshape(u_vec, m_x, m_y);

%% 計算誤差與繪圖
err = max(max(abs(u_soln - u_true)));
fprintf('Max error with hx=%f, hy=%f: %e\n', h_x, h_y, err);

surf(X,Y,u_soln)
xlabel('x'); ylabel('y'); zlabel('u');
title('Numerical Solution on Rectangular Domain');





































%(c) final
clear; clc; format long g;

% Domain and Mesh Definition
x1=0; x2=1; y1=0; y2=2;		
m_x=20; m_y=40;
h_x=(x2-x1)/(m_x+1);
h_y=(y2-y1)/(m_y+1);

x=linspace(x1,x2,m_x+2);
y=linspace(y1,y2,m_y+2);
[X,Y]=meshgrid(x,y);
X=X'; Y=Y';

% Interior indices
Iint = 2:m_x+1; Jint = 2:m_y+1;
Xint = X(Iint,Jint); Yint = Y(Iint,Jint);

% RHS and Exact Solution
f = @(x,y) 1.25*exp(x+y/2);
rhs = f(Xint,Yint);
u_true = exp(X+Y/2);
u_soln = u_true;

% Boundary condition adjustments
rhs(:,1) = rhs(:,1) - u_soln(Iint,1)/(h_y)^2;
rhs(:,m_y) = rhs(:,m_y) - u_soln(Iint,m_y+2)/(h_y)^2;
rhs(1,:) = rhs(1,:) - u_soln(1,Jint)/(h_x)^2;
rhs(m_x,:) = rhs(m_x,:) - u_soln(m_x+2,Jint)/(h_x)^2;

% Matrix Construction via Kronecker Products
F = reshape(rhs,m_x*m_y,1);
I_x = speye(m_x); I_y = speye(m_y);
e_x = ones(m_x,1);
T_x = spdiags([e_x -2*e_x e_x],[-1 0 1],m_x,m_x);
e_y = ones(m_y,1);
S_y = spdiags([e_y -2*e_y e_y],[-1 0 1],m_y,m_y);

A = kron(I_y,T_x)/(h_x)^2 + kron(S_y,I_x)/(h_y)^2;

% Solver
u_vec = A \ F;
u_soln(Iint, Jint) = reshape(u_vec, m_x, m_y);

% Error
error_grid = abs(u_soln - u_true);
max_err = max(error_grid(:));
fprintf('Max Error = %e\n', max_err);

% plot(Numerical Solution)
figure(1); clf;
surf(X, Y, u_soln);
colorbar;
xlabel('x'); ylabel('y'); zlabel('u');
title(['Numerical Solution']);
view(-45, 30);
grid on;
drawnow;

% plot(Error Surface)
figure(2); clf;
surf(X, Y, error_grid);
colorbar;
axis tight;
xlabel('x'); ylabel('y'); zlabel('|Error|');
title(['Absolute Error Distribution']);
grid on;
view(-45, 30);
drawnow;





























clear; clc; format long g;

% --- 1. Domain Definition (Same as Part c) ---
x1 = 0; x2 = 1; 
y1 = 0; y2 = 2; 

% --- 2. The "Strange Phenomenon" Handling ---
% We want dx = dy = h. Let's pick a target h.
h_target = 0.0505; 

% To avoid non-integer m_x or m_y, we round the counts:
m_x = round((x2 - x1) / h_target) - 1;
m_y = round((y2 - y1) / h_target) - 1;

% Re-calculate actual hx and hy. 
% In this specific domain [0,1]x[0,2], they will both be exactly 0.05.
hx = (x2 - x1) / (m_x + 1);
hy = (y2 - y1) / (m_y + 1);

fprintf('hx = %f, hy = %f\n', hx, hy);
fprintf('mx = %d, my = %d\n', m_x, m_y);

% --- 3. Mesh Generation ---
x = linspace(x1, x2, m_x+2);
y = linspace(y1, y2, m_y+2);
[X, Y] = meshgrid(x, y);
X = X'; Y = Y'; % Transpose to match matrix indexing

% Interior indices
Iint = 2:m_x+1; Jint = 2:m_y+1;
Xint = X(Iint, Jint); Yint = Y(Iint, Jint);

% --- 4. RHS and Boundary Conditions ---
f = @(x,y) 1.25*exp(x + y/2);
rhs = f(Xint, Yint);
u_true = exp(X + Y/2);
u_soln = u_true; % Dirichlet BCs

% Boundary adjustments for RHS
rhs(:, 1)   = rhs(:, 1)   - u_soln(Iint, 1) / hy^2;
rhs(:, m_y) = rhs(:, m_y) - u_soln(Iint, m_y+2) / hy^2;
rhs(1, :)   = rhs(1, :)   - u_soln(1, Jint) / hx^2;
rhs(m_x, :) = rhs(m_x, :) - u_soln(m_x+2, Jint) / hx^2;

% --- 5. Matrix Construction (b) ---
% Even though hx = hy, we use the general form for robustness
F = reshape(rhs, m_x*m_y, 1);
Ix = speye(m_x); Iy = speye(m_y);
ex = ones(m_x, 1); Tx = spdiags([ex -2*ex ex], [-1 0 1], m_x, m_x);
ey = ones(m_y, 1); Sy = spdiags([ey -2*ey ey], [-1 0 1], m_y, m_y);

A = kron(Iy, Tx)/hx^2 + kron(Sy, Ix)/hy^2;

% --- 6. Solve and Error Calculation ---
u_vec = A \ F;
u_soln(Iint, Jint) = reshape(u_vec, m_x, m_y);

max_err = max(max(abs(u_soln - u_true)));
fprintf('Max Error: %e\n', max_err);

% --- 7. Visualization ---
figure(3); clf;
surf(X, Y, u_soln);
colorbar; axis tight;
xlabel('x'); ylabel('y'); zlabel('u');
title(['Numerical Solution']);

figure(4); clf;
surf(X, Y, abs(u_soln - u_true));
colorbar; axis tight;
xlabel('x'); ylabel('y'); zlabel('|Error|');
title(['Absolute Error Distribution']);























clear; clc; format long g;

% 定義測試的 m 值 (網格點數)
m_values = [20, 40, 80, 160];
errors = zeros(size(m_values));
h_values = zeros(size(m_values));

fprintf(' m    |      h      |    Max Error (E_inf)  |  Order \n');
fprintf('-------------------------------------------------------\n');

for k = 1:length(m_values)
    m = m_values(k);
    a = 0; b = 1;
    h = (b-a)/(m+1);
    h_values(k) = h;
    
    x = linspace(a,b,m+2);
    y = linspace(a,b,m+2);
    [X,Y] = meshgrid(x,y);
    X = X'; Y = Y';
    
    Iint = 2:m+1; Jint = 2:m+1;
    Xint = X(Iint,Jint); Yint = Y(Iint,Jint);
    
    f = @(x,y) 1.25*exp(x+y/2);
    rhs = f(Xint,Yint);
    utrue = exp(X+Y/2);
    usoln = utrue; % 設定邊界條件
    
    % 調整 RHS 以包含邊界項
    rhs(:,1) = rhs(:,1) - usoln(Iint,1)/h^2;
    rhs(:,m) = rhs(:,m) - usoln(Iint,m+2)/h^2;
    rhs(1,:) = rhs(1,:) - usoln(1,Jint)/h^2;
    rhs(m,:) = rhs(m,:) - usoln(m+2,Jint)/h^2;
    
    % 建立 A 矩陣 (老師原本的 kron 邏輯)
    F = reshape(rhs, m*m, 1);
    I = speye(m);
    e = ones(m,1);
    T = spdiags([e -4*e e], [-1 0 1], m, m);
    S = spdiags([e e], [-1 1], m, m);
    A = (kron(I,T) + kron(S,I)) / h^2;
    
    % 解線性系統
    uvec = A \ F;
    usoln(Iint,Jint) = reshape(uvec, m, m);
    
    % 計算最大誤差
    errors(k) = max(max(abs(usoln - utrue)));
    
    % 計算觀察到的收斂階數 (Order)
    if k == 1
        order_str = '---';
    else
        order = log(errors(k)/errors(k-1)) / log(h_values(k)/h_values(k-1));
        order_str = sprintf('%.4f', order);
    end
    
    fprintf('%4d  |  %.6f  |      %.6e      |  %s \n', m, h, errors(k), order_str);
end


% --- 繪製 Log-Log Plot ---
figure(5); clf;
loglog(h_values, errors, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;

% 建立一條斜率為 2 的參考線 (Reference Line)
% 讓參考線通過第一個誤差點，方便觀察斜率是否平行
ref_line = errors(1) * (h_values / h_values(1)).^2;
loglog(h_values, ref_line, 'r--', 'LineWidth', 1.5);

grid on;
xlabel('Grid spacing (h)');
ylabel('Max Error (E_{\infty})');
title('Log-Log Plot of Error vs. Grid Spacing');
legend('Numerical Error', 'Theoretical O(h^2) Slope', 'Location', 'northwest');

% 顯示斜率數值
slope = (log(errors(end)) - log(errors(1))) / (log(h_values(end)) - log(h_values(1)));
text(h_values(2), errors(2), sprintf('  Observed Slope = %.4f', slope), 'FontSize', 12);





























































