function [te0_field, te1_field] = get_modes(wv1, wv2, x, z)

c = 3e8;
um = 1e-6;
lambda = 1.55*um;
k0 = 2*pi/lambda;
width = 4*um;

%%% TE0 %%%
%%% CHANGE PARAM HERE
wavevector = wv1;
n_eff = 1.52114;
n_clad = 1.445918;
beta = sqrt(k0^2*n_eff^2-wavevector^2);
gamma = sqrt(k0^2*(n_eff^2 - n_clad^2) - wavevector^2);
x_core = x(x>-width/2 & x<width/2);
x_1 = x(x<=-width/2);
x_2 = x(x>=width/2);
e_core = cos(wavevector*x_core)/cos(wavevector*width/2);
e_1 = exp(-gamma*(abs(x_1) - width/2));
e_2 = exp(-gamma*(abs(x_2) - width/2));

e = [e_1 e_core e_2];
e = 1/max(e) * e * exp(1i*beta*z);
te0_field = e;


%%% TE1 %%%
%%% CHANGE PARAM HERE
wavevector = wv2;
n_eff = 1.52114;
n_clad = 1.445918;
beta = sqrt(k0^2*n_eff^2-wavevector^2);
gamma = sqrt(k0^2*(n_eff^2 - n_clad^2) - wavevector^2);
x_core = x(x>-width/2 & x<width/2);
x_1 = x(x<=-width/2);
x_2 = x(x>=width/2);
e_core = sin(wavevector*x_core)/sin(wavevector*width/2);
e_1 = exp(-gamma*(abs(x_1) - width/2));
e_2 = -exp(-gamma*(abs(x_2) - width/2));

e = [e_1 e_core e_2];
e = 1/max(e) * e * exp(1i*beta*z);
te1_field = e;

end

