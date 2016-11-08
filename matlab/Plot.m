%%% -------------------------------------------------- %%%
%%% Author: Denys Dutykh, CNRS -- LAMA, Univ of Savoie %%%
%%% E-mail: Denys.Dutykh@univ-savoie.fr                %%%
%%% Web:    http://www.denys-dutykh.com/               %%%
%%% Blog:   http://dutykh.github.io/                   %%%
%%% GitHub: https://github.com/dutykh/                 %%%
%%% -------------------------------------------------- %%%

function Plot (Om, t, frame)

	global dy Lx Ly X Y
	
	surf(X, Y, Om), grid off
	shading interp;
	colormap(jet); cc = colorbar;
    xlim([-Lx Lx]); ylim([-Ly+dy Ly]); caxis([-1 1]);
    xlabel('$x$', 'interpreter', 'latex', 'fontsize', 12);
    ylabel('$y$', 'interpreter', 'latex', 'fontsize', 12, 'Rotation', 1);
    xlabel(cc, '$\omega(x,y,t)$', 'interpreter', 'latex', 'fontsize', 12, 'Rotation', 90);
    view([0 90]);

    title (['Vorticity distribution at t = ',num2str(t,'%4.2f')], 'interpreter', 'latex', 'fontsize', 12);

    set(gcf, 'Color', 'w');
    drawnow
end % Plot ()