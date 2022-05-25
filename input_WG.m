function [output_wave] = input_WG(core_index, width)
PlotEachStep = false; %set true to see the field u^m each step; 
    %much faster when false! useful perhaps for debugging. 
SecondWG = false; 
    %set true to include a second waveguide, e.g. for a co-directional
    %coupler. For your projects you'll make similar modifications to the
    %index profile to describe your structure, including along the z
    %direction (not implemented here). 
um=1e-6;
i = sqrt(-1);

% Parameters of the simulation

% Free space wavelength
lambda = 1.55*um;
k0 = 2*pi/lambda;

% Indices of refraction
nCladding = 1.445918;
nCore = core_index;
nBar = (nCladding + nCore) / 2;

% Length of waveguide (z-direction)
inputWGLength = 500*um;%400*um;

% Width of the domain (x-direction)
widthDomain = 12*um; % 
% Length of simulation domain (z-direction)
lengthDomain = inputWGLength;

% Width of waveguide (x-direction)
inputWGWidth = 4*um;

% gap between waveguides: 
coupGap = 0.5*um; 

% Width of initial Gaussian
sig = width*um;

% Discretization in z-direction
deltaz = .25*um;
Nzpts = round(lengthDomain/deltaz); %excluding z=0; 

% Discretization in x-direction
deltax = 0.01 * um;
% Number of points in x-direction resulting from discretization
N = round(widthDomain / deltax);

% Create the x-axis and z-axis: 0 is at center of domain
x = linspace(-widthDomain/2, widthDomain/2, N); 
z = linspace(0, deltaz*Nzpts, Nzpts+1); 

% Index profile for WG 
n_Input_WG = nCladding*ones(1, N); 
coreinds = find((x<=inputWGWidth/2).*(x>=-inputWGWidth/2)); 
    %find returns nonzero element indices; the boolean expressions 
    %identifies values of x for which we are in the waveguide core. 
n_Input_WG(coreinds) = nCore; 
if(SecondWG)
    core2inds = find((x<=(3*inputWGWidth/2 + coupGap)).* ...
    (x>=(inputWGWidth/2 + coupGap))); 
    n_Input_WG(core2inds) = nCore; 

    %core3inds = find((x>=-(3*inputWGWidth/2 + coupGap)).* ...
    %(x<=-(inputWGWidth/2 + coupGap))); 
    %n_Input_WG(core3inds) = nCore; 
    % this last bit is if we want a third waveguide in, just for
    % illustration. You'll want to remove these bits and describe the index
    % profile particular to your own device. 
end 
% Initial Gaussian that will be launched
u = exp(-(x/sig).^2);

% allField records all steps of u (at increasing values of m): 
allField = zeros(Nzpts+1, N); 
allField(1,:) = u;

%
% Parabolic Absorbing Boundary Domains
% Region on either side of domain that will absorb field
widthAbsEdge = 2*um; % delta
% Max value of absorption coefficient
kappa_max = -0.03;
% Define the kappa vector: it's 0 in the central region,
% and has value kappa in outer edges
kappa = zeros(1, N); 
for j = 1 : N
    if x(j) < -widthDomain/2+widthAbsEdge
        kappa(j) = kappa_max*((widthDomain/2 - ...
            widthAbsEdge + x(j))/widthAbsEdge)^2;
    elseif x(j) > widthDomain/2-widthAbsEdge
        kappa(j) = kappa_max*((x(j) - widthDomain/2 + ...
            widthAbsEdge)/widthAbsEdge)^2;
    else % Central region: no absorption
        kappa(j) = 0;
    end
end

% Adapted from Pedrola, "Beam Propagation Method for Design of Optical
% Waveguide Devices" 
alpha = 0.5; % Scheme Parameter; 0(1) means purely Forward (Backward)

% Crank-Nicolson intermediate parameters (Pedrola, p. 36): 
a = -alpha / (deltax)^2;
c = a;

% MAIN CN LOOP
r = zeros(1, N); 
for m = 1 : Nzpts 
    b = (2*alpha)/(deltax^2) - alpha*(n_Input_WG.^2 - kappa.^2 + ...
        2*i*n_Input_WG.*kappa - nBar^2)*k0^2 + 2*i*k0*nBar/deltaz;
    % rj values
    r(1) = (1-alpha)/(deltax^2) * (0+u(2)) + ...
        ((1-alpha)*((n_Input_WG(1))^2-kappa(1)^2 + 2*i*n_Input_WG(1)*kappa(1) - ...
        nBar^2)*k0^2 - 2*(1-alpha)/(deltax^2) + 2*i*k0*nBar/deltaz) * u(1);
    r(N) = (1-alpha)/(deltax^2) * (u(N-1)+0) + ...
        ((1-alpha)*(n_Input_WG(N)^2-kappa(N)^2 + 2*i*n_Input_WG(N)*kappa(N) - ...
        nBar^2)*k0^2 - 2*(1-alpha)/(deltax^2) + 2*i*k0*nBar/deltaz) * u(N);
    for j = 2 : N-1
        r(j) = (1-alpha)/(deltax^2) * (u(j-1)+u(j+1)) + ...
            ((1-alpha)*(n_Input_WG(j)^2-kappa(j)^2 + 2*i*n_Input_WG(j)*kappa(j)- ...
            nBar^2)*k0^2 - 2*(1-alpha)/(deltax^2) + 2*i*k0*nBar/deltaz) * u(j);
    end
    % Thomas algorithm to solve tridiagonal Crank Nicolson scheme
    beta = b(1);
    u(1) = r(1) / beta;
    for j = 2 : N
        gamma(j) = c / beta;
        beta = b(j) - a * gamma(j);
        u(j) = (r(j) - a * u(j-1)) / beta;
    end
    for j = 1 : N-1
        ktemp = N - j;
        u(ktemp) = u(ktemp) - gamma(ktemp+1) * u(ktemp+1);
    end
    titlestring = strcat('Intensity in Input WG after ', '{ }', ...
        num2str(m*deltaz/um), ' microns');
    
    % Keep track of all fields every recordFieldStep microns
    allField(m+1, :) = u;
    
    if(PlotEachStep)
        figure(2)
        plot(x, abs(u).^2, 'linewidth', 3);
        hold on
        title(titlestring)
        axis([min(x) max(x) -0.1 2])
        hold off
        shading interp
        drawnow
    end
end % end of the loop: for m = 1 : round(lengthDomain/deltaz)

output_wave = u;

% figure;
% surf(x, z, abs(allField).^2)
% view(0, 90); shading interp; axis tight; colorbar 
% xlabel('z (m)', 'fontSize', 14);
% ylabel('x (m)', 'fontSize', 14);
end 