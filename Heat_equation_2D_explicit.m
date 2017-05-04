%% 2D Diffusion Equation
% Explicit Method
% The script solves:
%      du/dt = d^2u/dx^2 + d^2u/dy^2
% Using the Explicit method

% Domain
%   x = [a_x,b_x]
%   y = [a_y,b_y]

% BC
%   u(x,y,0) = 0
%   u(a_x,y,t) = f_b(y),     u(b_x,y,t) = g_b(y)
%   du/dy(b_y) = 0
%   u(x,a_y,t) = f_b(a_y) + [(x-a_y)/(b_x - a_x)][g_b(a_y - f_b(a_y)]
%   a_x = a_y = 0,           b_x = b_y = 2pi
%   f_b(y) = (b_y - y)^2*cos(pi*y/b_y)
%   g_b(y) = y(b_y - y)^2

debugg = load('Debuggexp.txt');
D = debugg(1);
steps_x = debugg(2);
steps_y = debugg(3);
T = debugg(4);
eta = debugg(5);
a_x = debugg(6); b_x = debugg(7);
a_y = debugg(8); b_y =debugg(9);

% % Given information
% D = input('Input value of D ([] otherwise):');
% steps_x = input('Input number of x steps (N_x):');
% steps_y = input('Input number of y steps (N_y):');
% T = input('Input time length T:');
% eta = input('Input time-step selection eta(between 0 and 1 or []):');
% 
% %   Domain
% x_int = input('Input the x interval(a_x,b_x):');
% y_int = input('Input the y interval(a_y,b_y):');
% 
% %   Check
% if isempty(D), D = 1; end
% if isempty(eta), eta = 0.9; end
% if isempty(x_int) | isempty(y_int)
%     error('Error: One or more domain intervals are missing.')
% end
% 
% a_x = x_int(1);
% b_x = x_int(2);
% a_y = y_int(1);
% b_y = y_int(2);

% Spartial discritization
dx = (b_x-a_x)/steps_x;
x = a_x:dx:b_x;
nx = length(x);

dy = (b_y-a_y)/steps_y;
y = a_y:dy:b_y;
ny = length(y);

% Time steps.
dt = (.5*(dx*dy)^2*eta)/(D*(dx^2 + dy^2));
t = 0:dt:T;
n_t = length(t);

% BC
f_b = (b_y - y).^2.*cos((pi/b_y).*y);
g_b = y.*(b_y - y).^2;

u_ax = f_b;
u_bx = g_b;
  f_b_ay = (b_y - a_y)^2*cos(pi*a_y/b_y);
  g_b_ay = a_y*(b_y - a_y)^2;
u_ay = f_b_ay + ((x - a_x)/(b_x -a_x))*(g_b_ay - f_b_ay);

% Neumann BC
% Create a ghost node
% du/dx = v
% u_ny+2 = 2dy*v + u_ny
% Plug this in when k = ny + 1
v = 0;

% IC
u_t0 = 0;

% Parameters
lamb_x = dt*D/dx^2;
lamb_y = dt*D/dy^2;

% Initialize u the approximate solution.
%   Order of the form y by x
uex = [u_ay;
    u_ax(2:end)' zeros(ny -1,nx -2) u_bx(2:end)'];

u_new = zeros(ny,nx);   % The new solution

% Create the grid
[X,Y] = meshgrid(x,y);

tym = 0;

for i = 1:n_t % traversing time
    
    for k = 2:(nx-1)   % traversing x direction
        for j = 2:(ny-1)   % traversing the y direction
            u_new(j,k) = lamb_x*(uex(j-1,k) + uex(j+1,k)) + (1 - 2*lamb_x - 2*lamb_y)*uex(j,k) + lamb_y*(uex(j,k-1) + uex(j,k+1));
        end
    end
    
% Enforce the BC
    u_new(1,:) = u_ay;
    u_new(2:ny,1) = u_ax(2:ny)';
    u_new(2:ny,nx) = u_bx(2:ny)';
    
    for k = 2:nx-1
        u_new(ny,k) = lamb_x*(2*uex(ny-1,k) + 2*v*dy) + (1 - 2*lamb_x - 2*lamb_y)*uex(ny,k) + lamb_y*(uex(ny,k-1) + uex(ny,k+1));
        %u(ny,k) + (dt/dx^2)*(2*u(ny-1,k) -2*u(ny,k) + 2*dy*v) + (dt/dy^2)*(u(ny,k-1) -2*u(ny,k) + u(ny,k+1));
    end
    uex = u_new;
    tym = tym +1;
    
    if mod(i,500)==0
        figure (1), clf
        pcolor(X/1e3,Y/1e3,uex); shading interp,  colorbar
        hold on
        contour(X/1e3,Y/1e3,uex);
        xlabel('x')
        ylabel('y')
        zlabel('approximation')
        title('Explicit')
        drawnow
        
        figure (2)
        surf(X,Y,uex)
        xlabel('x')
        ylabel('y')
        zlabel('approximation')
        title('Explicit')
        drawnow
    end
    
    tym = tym + 1;
end

lamb = D*dt*(1/dx^2 + 1/dy^2);  % This is the stability condition: lamb < .5