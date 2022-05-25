%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Toggles
%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotEachStep = false;
PlotEachStepN = false;
SecondWG = false; 
PlotNProfile = true;
PlotKappa = false;
consistent = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%
um=1e-6;
% Free space wavelength
lambda = 1.55*um;
k0 = 2*pi/lambda;

% Indices of refraction 
nCladding = 1.445918; 
nCore = 1.52114;
nBar = (nCladding + nCore) /2;

% Length of waveguide (z-direction)
inputWGLength = 3000*um;

% Width of the domain (x-direction)
widthDomain = 12*um; % 
% Length of simulation domain (z-direction)
lengthDomain = 2500*um;

% Width of waveguide (x-direction)
inputWGWidth = 4*um;

% Width of initial Gaussian
sig = 2.3;

%%%%%%%%%%%%%%%%
% Discretization
%%%%%%%%%%%%%%%%%

% Discretization in z-direction
deltaz = 0.25*um;
Nzpts = round(lengthDomain/deltaz);

% Discretization in x-direction
deltax = 0.01*um;
% Number of points in x-direction resulting from discretization
N = round(widthDomain / deltax);

% Create the x-axis and z-axis: 0 is at center of domain
x = linspace(-widthDomain/2, widthDomain/2, N); 
z = linspace(0, deltaz*Nzpts, Nzpts+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refractive Index Profile
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Index profile for WG 
n_Input_WG = nCladding*ones(Nzpts + 1, N);
counter = 0;
L = 500*um;
A = 63.8*um;
width_perturbation = 0.325*um;
offset = 500; % Offset for where perturbation starts
last_perturbation_index = 1;
for j = 1 : Nzpts + 1
    if (j < offset)
        coreinds = find((x<=inputWGWidth/2).*(x>=-inputWGWidth/2)); 
        n_Input_WG(j,coreinds) = nCore; 
        continue
    end
    if((counter<(A/deltaz)*0.5)&&(j<L/deltaz + offset))
        coreinds = find((x<=inputWGWidth/2).*(x>=-inputWGWidth/2));  
        n_Input_WG(j,coreinds) = nCore; 
        perturbationinds = coreinds(1:round(width_perturbation/deltax));
        n_Input_WG(j,perturbationinds) = nCladding;
        last_perturbation_index = j;
    else
        coreinds = find((x<=inputWGWidth/2).*(x>=-inputWGWidth/2)); 
        n_Input_WG(j,coreinds) = nCore; 
    end
    if(counter>=(A/deltaz))
         counter = 0;
    end
    counter = counter+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Refractive Index Profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PlotNProfile
    [X, Z] = meshgrid(x, z);
    figure;
    surf(X, Z, n_Input_WG);
    view(0,90); shading interp; axis tight; colorbar; 
    title("Refractive Index Profile");
    xlabel("x (m)");
    ylabel("z (m)");
end

% Initial Gaussian that will be launched
% u = exp(-(x/sig).^2);
u = input_WG(nCore, sig);

% allField records all steps of u
allField = zeros(Nzpts+1, N);  
allField(1,:) = u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parabolic Absorbing Boundary Domains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Region on either side of domain that will absorb field
widthAbsEdge = 2*um; 
% Define the kappa vector: it's 0 in the central region,
% and has value kappa in outer edges
kappa_max = -0.03;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Kappa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(PlotKappa)
    figure;
    plot(x, kappa, '*')
    title('Representation of absorbing region')
end

% Adapted from Pedrola, "Beam Propagation Method for Design of Optical
% Waveguide Devices" 
alpha = 0.5; % Scheme Parameter; 0(1) means purely Forward (Backward)

% Crank-Nicolson intermediate parameters (Pedrola, p. 36): 
a = -alpha / (deltax)^2;
c = a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crank–Nicolson Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
te0_power = zeros(1, Nzpts);
te1_power = zeros(1, Nzpts);
r = zeros(1, N);

for m = 1 : Nzpts
    n_m_plus_1 = n_Input_WG(m + 1, :);
    n_m = n_Input_WG(m, :);
    b = (2*alpha)/(deltax^2) - alpha*(n_m_plus_1.^2 - kappa.^2 + ...
        2*1i*n_m_plus_1.*kappa - nBar^2)*k0^2 + 2*1i*k0*nBar/deltaz;
    % rj values
    r(1) = (1-alpha)/(deltax^2) * (0+u(2)) + ...
        ((1-alpha)*((n_m(1))^2-kappa(1)^2 + 2*1i*n_m(1)*kappa(1) - ...
        nBar^2)*k0^2 - 2*(1-alpha)/(deltax^2) + 2*1i*k0*nBar/deltaz) * u(1);
    r(N) = (1-alpha)/(deltax^2) * (u(N-1)+0) + ...
        ((1-alpha)*(n_m(N)^2-kappa(N)^2 + 2*1i*n_m(N)*kappa(N) - ...
        nBar^2)*k0^2 - 2*(1-alpha)/(deltax^2) + 2*1i*k0*nBar/deltaz) * u(N);
    for j = 2 : N-1
        r(j) = (1-alpha)/(deltax^2) * (u(j-1)+u(j+1)) + ...
            ((1-alpha)*(n_m(j)^2-kappa(j)^2 + 2*1i*n_m(j)*kappa(j)- ...
            nBar^2)*k0^2 - 2*(1-alpha)/(deltax^2) + 2*1i*k0*nBar/deltaz) * u(j);
    end
    % Thomas algorithm to solve tridiagonal Crank-Nicolson scheme
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
    
    allField(m+1, :) = u;
    
    % Overlap integral calculation
    % TE0
    [te0_field, te1_field] = get_modes(629108, 1.77991e6, x, 500*um + z(m));
    self_overlap_u = abs(trapz(x, u.*conj(u)));
    self_overlap_0 = abs(trapz(x, te0_field.*conj(te0_field)));
    self_overlap_1 = abs(trapz(x, te1_field.*conj(te1_field)));
    c_u = 1/sqrt(self_overlap_u);
    c_0 = 1/sqrt(self_overlap_0);
    c_1 = 1/sqrt(self_overlap_1);
    
%     te0_coupled_power = abs(trapz(x, u.*conj(te0_field)))^2 / trapz(x, abs(conj(te0_field)).^2);
    te0_coupled_power =  abs(trapz(x, c_u*u.*conj(c_0*te0_field)));
    te0_power(m) = te0_coupled_power;
    
%     te1_coupled_power = abs(trapz(x, u.*conj(te1_field)))^2 / trapz(x, abs(conj(te1_field)).^2);
    te1_coupled_power =  abs(trapz(x, c_u*u.*conj(c_1*te1_field)));
    te1_power(m) = te1_coupled_power;
    
    if(PlotEachStep)
        figure(2);
        plot(x, abs(c_u*u).^2, 'linewidth', 3);
        hold on
        plot(x, abs(c_0*te0_field).^2, 'linewidth', 3);
        hold on
        plot(x, abs(c_1*te1_field).^2, 'linewidth', 3);
        legend("U", "TE0", "TE1");
        title(titlestring)
        hold off
        shading interp
        drawnow
    end
    
    if(PlotEachStepN)
        index_string = strcat('Index Slice after ', '{ }', ...
            num2str(m*deltaz/um), ' microns');
        figure(3);
        plot(x, n_m, 'linewidth', 3);   
        hold on
        plot(x, n_Input_WG(end, :), 'linewidth', 3); 
        title(index_string);
        hold off
        xlabel("x(m)");
        ylabel("n");
        legend("n_m", "n_{last}");
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot E-Field Profile in X-Z Plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, ~, ~] = mkdir(strcat('index_', num2str(nCore)));
fname = strcat('/Users/robinlin/Desktop/ECE4370_Final_Project/', 'index_', num2str(nCore));
filename = strcat('index_', num2str(nCore), '_width_', num2str(sig), '.png');
fullfile_name = fullfile(fname, filename);

if(consistent)
    figure;
    surf(z,x, abs(allField').^2)
    hold on
    xline(z(offset),'-',{'Grating Start'}, 'lineWidth', 3, 'Color', 'w');
    hold on
    xline(z(last_perturbation_index),'-',{'Grating End'}, 'lineWidth', 3, 'Color', 'w');
    view(0, 90); shading interp; axis tight; colorbar 
    xlabel('x (m)', 'fontSize', 20);
    ylabel('y (m)', 'fontSize', 20);
    saveas(gcf, fullfile_name, 'png');
else
    figure;
    surf(x, z, abs(allField).^2)
    view(0, 90); shading interp; axis tight; colorbar 
    xlabel('z (m)', 'fontSize', 20);
    ylabel('x (m)', 'fontSize', 20);
	saveas(gcf, fullfile_name, 'png');
end

figure;
plot(z(1:end-1), te0_power);
hold on
plot(z(1:end-1), te1_power);
hold on
xline(z(offset),'-',{'     Grating Start'}, 'lineWidth', 3, 'Color', 'b', 'fontSize', 15);
hold on
xline(z(last_perturbation_index),'-',{'     Grating End'}, 'lineWidth', 3, 'Color', 'b', 'fontSize', 15);
xlabel("Propagation Direction, z (m)", 'fontSize', 15);
ylabel("Mode Power", 'fontSize', 15);
title("Mode Power in Waveguide", 'fontSize', 15);
legend("TE0", "TE1", 'fontSize', 15);