clc
clear
close all

% Parameters
f_center = 2.4e9;                                               %center frequency
Rp = input('passband ripple(in dB)= ');                         % 0.04321 for n = 5
Rs = input('Enter minimum attenuation in stopband(in dB)= ');   % 30 for n= 5
BW = 144e6;                                                     % bandwidth
R0 = 50;                                                        % source impedance in ohm
Wsamp = 10e9;
Wp = [f_center-(BW/2) f_center+(BW/2)]/(Wsamp/2);
Ws = [f_center-BW f_center+BW]/(Wsamp/2);
FBW = (BW)/f_center;

p = input('Enter number of of points= ');
step_size = 2*BW/p;


% Calculate the ripple factor
epsilon = sqrt(10^(Rp/10) - 1);

%Design the filter
[n,Wp] = cheb1ord(Wp, Ws, Rp, Rs);
[b,a] = cheby1(n,Rp,Wp);

string = "The order of the filter is ";
disp(string + int2str(n))

freqz(b,a,512,Wsamp);

% Calculate the beta and gamma parameter
beta = log(coth(Rp/17.37));
gamma = sinh(beta/(2*n));

% Calculate the g-elements
disp('g(0) = 1.0')

for i = 1:n
    if i ==1
        g(i) = (2/gamma)*sin(pi/(2*n));
    else
        g(i) = (1/g(i-1))*(4*sin((2*i-1)*pi/(2*n)))*sin((2*i-3)*pi/(2*n))...
            /(gamma.^2 + (sin((i-1)*pi/n)).^2);
    end
    disp(['g(',int2str(i),') = ', num2str(g(i))]);
end

if mod(n,2) == 0
    Rl = (coth(beta/4))^2;
else
    Rl = 1.0;
end
disp(['g(',int2str(n+1),')= ',num2str(Rl)])


% Calculating coupling coefficients
m = zeros(n);
for i = 1:n-1
    m(i,i+1) = FBW/sqrt(g(i)*g(i+1));
    m(i+1,i) = m(i,i+1);
end
m1 = m/FBW;
disp('coupling coefficient matrices(m)')
disp(m)

% Calculate Quality factor
Q = zeros(n);
Q(1,1) = 1*g(1)/FBW;
Q(n,n) = g(n)*1/FBW;

% q factor
qe1 = Q(1,1)*FBW;
qen = Q(n,n)*FBW;
q = zeros(n);
q(1,1) = 1/qe1;
q(n,n) = 1/qen;
disp('Quality factor(Q)')
disp(Q)

% s parameter calculation
U = eye(n);

for k = 1:p
    f(k) = (f_center-(BW)) + k*step_size;
    p_val = (1i/FBW)*((f(k)/f_center)-(f_center/f(k)));
    A = q +(p_val).* U - 1j * m1;
    B = inv(A);
    s21(k) = 2*(1/sqrt(qe1*qen))*B(n,1);
    s21_db(k) = 20*log(abs(s21(k)));

    s11(k) = ( 1-(2*(1/qe1))*B(1,1));
    s11_db(k) = 20*log(abs(s11(k)));
end

% Plot the results
figure;
%subplot(211);
plot(f, s11_db, 'r');
title('S11 Parameter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%subplot(212);
figure;
plot(f, s21_db, 'b');
title('S21 Parameter');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');



