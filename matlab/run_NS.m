%%% -------------------------------------------------- %%%
%%% 2D Navier-Stokes pseudo-spectral solver            %%%
%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

clear, close all
format longE
addpath 'tools/'

% ------------------------------------------------------------ %

% declaration of global variables:
global dy nu Dx Dy Kx Ky K2 Lx Ly Nx Ny Nxy X Y

% ------------------------------------------------------------ %

%%% Physical parameters:
nu = 1e-4;	% 1/Re or viscosity

% ------------------------------------------------------------ %

%%% Domain definition:
Lx = pi;    % Domain half-length in x-direction
Ly = pi;    % Domain half-length in y-direction

%%% Numerical parameters:
Nx = 1024;  % number of Fourier modes in discrete solution x-dir
Ny = 1024;	% number of Fourier modes in discrete solution y-dir
Nxy = Nx*Ny;

dx = 2*Lx/Nx;   		% distance between two physical points
x = (1-Nx/2:Nx/2)'*dx;  % physical space discretization

dy = 2*Ly/Ny;   		% distance between two physical points
y = (1-Ny/2:Ny/2)'*dy;  % physical space discretization

[X,Y] = meshgrid(x,y);	% 2D composed grid

% ------------------------------------------------------------ %

% vectors of wavenumbers in the transformed space:
kx = [0:Nx/2 1-Nx/2:-1]'*pi/Lx;
ky = [0:Ny/2 1-Ny/2:-1]'*pi/Ly;

% antialising treatment
jx = (Nx/4+2:Nx/4*3);  % the frequencies we sacrify
kx(jx) = 0;

jy = (Ny/4+2:Ny/4*3);  % the frequencies we sacrify
ky(jy) = 0;

% ------------------------------------------------------------ %

%%% Some operators arising in NS equations:
[Kx, Ky] = meshgrid(kx,ky);
K2 = Kx.^2 + Ky.^2;     % to compute the Laplace operator

K2inv = zeros(size(K2));
K2inv(K2 ~= 0) = 1./K2(K2 ~= 0);

Dx = 1i*Kx.*K2inv;		% u velocity component reconstruction from the vorticity
Dy = 1i*Ky.*K2inv;		% v velocity component reconstruction from the vorticity

% ------------------------------------------------------------ %

fftw('planner', 'hybrid');

%%% Set random number generator (for the initial condition)
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% ------------------------------------------------------------ %

%%% Time-stepping parameters:
t = 0.0;           	% the discrete time variable
Tf = 5.0;          	% final simulation time
ds = 0.1;			% write time of the results

ops = odeset('RelTol', 1e-10, 'AbsTol', 1e-10, 'OutputFcn', @odetpbar);

% ------------------------------------------------------------ %

Om = InitCond();

FigHandle = figure(1);
set(FigHandle, 'Position', [100, 100, 600, 550]);

frame = 0;
Plot(Om, t, frame);	% we plot the initial condition

% ------------------------------------------------------------ %

Om_vec = reshape(Om, Nxy, 1);
while (t < Tf) % main loop in time
	[~, v] = ode113(@RHS, [0:ds:ds], Om_vec, ops);

	Om_vec = v(end,:);
	Om = reshape(Om_vec, Nx, Ny);

	t = t + ds; frame = frame + 1;
	Plot(Om, t, frame); % and we plot the solution
end % while (t)