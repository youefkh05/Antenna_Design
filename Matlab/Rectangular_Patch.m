% MATLAB Program to find the basic parameters of a Probe-Fed Rectangular Microstrip Patch Antenna (by Ramann Mantha).

% Input Parameters:
epsilon_r = 3.55;
tan_delta = 0.0027;
h = 2e-3;
L =  2.795e-3;
W = 4.969e-3;
sigma = 5.8*(10^7);
a = 0.635*(10^(-3));
x0 = 1.85*(10^(-2));

% Fundamental Constants:
c = 2.99792458*(10^8);
mu0 = 4*pi*(10^(-7));
epsilon0 = 1/((c^2)*mu0);
eta0 = sqrt(mu0/epsilon0);

% Formulas to Calculate Resonant Frequency & Fringing:
epsilonr_effective = ((epsilon_r+1)/2)+(((epsilon_r-1)/2)*(1/(sqrt(1+12*(h/W)))));
delta_L = 0.412*h*(((W/h)+0.262)/((W/h)+0.813))*((epsilonr_effective+0.3)/(epsilonr_effective-0.258));
delta_W = h*(log(4)/pi);
W_effective = W+(2*delta_W);
L_effective = L+(2*delta_L);
lambda0 = 2*(L_effective)*(sqrt(epsilon_r));
f0_GHz = c/(lambda0*(10^9));
omega = 2*pi*(f0_GHz)*(10^9);
k0 = omega*(sqrt(mu0*epsilon0));
k1 = k0*(sqrt(epsilon_r));

% Formulas to Calculate Probe Reactance:
Gamma = 0.577216;
Xp = (eta0/(2*pi))*(k0*h)*(-Gamma+(log(2/((k0*a)*(sqrt(epsilon_r))))));

% CAD Formulas for Patch Properties:
c1 = 1-(1/epsilon_r)+((2/5)/((epsilon_r)^2));
c2 = -0.0914153;
a2 = -0.16605;
a4 = 0.00761;
intermediate1 = 1+((a2/10)*((k0*W_effective)^2))+((((a2)^2)+(2*a4))*(3/560)*((k0*W_effective)^4));
intermediate2 = (c2*(1/5)*((k0*L_effective)^2))+(a2*c2*(1/70)*((k0*W_effective)^2)*((k0*L_effective)^2));
p = intermediate1+intermediate2;
erhed = 1/(1+((3/4)*pi*(k0*h)*(1/c1)*((1-(1/epsilon_r))^3)));
Rs = sqrt((omega*mu0)/(2*sigma));
intermediate3 = 1+(erhed*(tan_delta+((Rs/(pi*eta0))*(1/(h/lambda0))))*((3/16)*(epsilon_r/(p*c1))*(L_effective/W_effective)*(1/(h/lambda0))));
er = erhed/intermediate3;
er_percentage = er*100;
BW = (1/sqrt(2))*(tan_delta+((Rs/(pi*eta0))*(1/(h/lambda0)))+((16/3)*((p*c1)/epsilon_r)*(h/lambda0)*(W_effective/L_effective)*(1/erhed)));
BW_percentage = BW*100;
Q = 1/(sqrt(2)*BW);
intermediate4 = (tan_delta+((Rs/(pi*eta0))*(1/(h/lambda0)))+((16/3)*((p*c1)/epsilon_r)*(h/lambda0)*(W_effective/L_effective)*(1/erhed)));
R_edge = ((4/pi)*(eta0)*(L_effective/W_effective)*(h/lambda0))/intermediate4;
R = R_edge*((cos((pi*(x0+delta_L))/L_effective))^2);
D = 3*(1/(p*c1))*((tan(k1*h)/(k1*h))^2)*(1/(1+((1/epsilon_r)*((tan(k1*h))^2))));
D_dB = 10*log10(D);

% Formulas for Circuit Elements:
R_circuit = R;
L_circuit = R/(omega*Q);
C_circuit = 1/(L_circuit*(omega^2));
Lp_circuit = Xp/omega;

% To Display the Ouput Parameters:
fprintf('\nf0_GHz = %f\n', f0_GHz);
fprintf('BW_percentage = %f\n', BW_percentage);
fprintf('Q = %f\n', Q);
fprintf('er_percentage = %f\n', er_percentage);
fprintf('R_edge = %f\n', R_edge);
fprintf('R = %f\n', R);
fprintf('Xp = %f\n', Xp);
fprintf('D_dB = %f\n', D_dB);

fprintf('\nR_circuit = %f\n', R_circuit);
fprintf('L_circuit = %d\n', L_circuit);
fprintf('C_circuit = %d\n', C_circuit);
fprintf('Lp_circuit = %d\n', Lp_circuit);

% To Plot the Input Impedance:
f_GHz = (1.52:0.002:1.62);
intermediate5 = (1/R_circuit)+(1./(1j*2*pi*f_GHz*(10^9)*L_circuit))+(1j*2*pi*f_GHz*(10^9)*C_circuit);
Zin = (1j*2*pi*f_GHz*(10^9)*Lp_circuit)+(1./intermediate5);
Rin = real(Zin);
Xin = imag(Zin);
plot(f_GHz,Rin,'r',f_GHz,Xin,'b');
title('Plot of Input Impedance vs Frequency');
xlabel('Frequency in GHz');
ylabel('Input Resistance & Input Reactance');
legend('Rin','Xin');
grid on;