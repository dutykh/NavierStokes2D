%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function rhs = RHS (t, Om_vec)

	global nu Dx Dy Kx Ky K2 Nx Ny Nxy

	Om = reshape(Om_vec, Nx, Ny);
	Om_hat = fft2(Om);

	Omx = real(ifft2(1i*Kx.*Om_hat));
	Omy = real(ifft2(1i*Ky.*Om_hat));

	u = real(ifft2(Dy.*Om_hat));
	v = real(ifft2(-Dx.*Om_hat));

	rhs = reshape(real(ifft2(-nu*K2.*Om_hat)) - u.*Omx - v.*Omy, Nxy, 1);
end % RHS ()