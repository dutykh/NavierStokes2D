%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function Om = InitCond ()

	global Nx Ny

	Om_hat = zeros(Nx, Ny);

	% We take a random distribution of the vorticity:
	Om_hat(1,5) = randn(1) + 1i*randn(1);
	Om_hat(2,2) = randn(1) + 1i*randn(1);
	Om_hat(4,1) = randn(1) + 1i*randn(1);

	Om = real(ifft2(Om_hat));
	Om = Om/max(max(Om));	% normalize to O(1)
end % InitCond ()