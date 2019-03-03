%% Function-graph application section
f = @(t) 2.0*cos(2*pi*5*t)+ 0.8*sin(2*pi*12*t) + 0.3*cos(2*pi*47*t);
domain=[0,1];
N = 128; % number of samples

% sample f N times over the domain
delta = (domain(2)/128);
Xs = domain(1):delta:(domain(2)-delta);
Ys = f(Xs);
% plot our samples over the graph
hold on
plot(0:.001:1, f(0:.001:1));
plot(Xs,Ys,'kx');
hold off

%% Decomposition into basic waveforms using the definitions (N^2)
close all
c = zeros(length(Xs));
for i=1:N
    for j=1:N
        c(i) = c(i) + Ys(j)*exp(-2*pi* 1i * (i-64)*j/128);
    end
    c(i) = c(i) / N; 
end

plot((-N/2+1):(N/2),abs(c));
title('Frequency Response (Energy at each Frequency)');

%% Decomposition into basic waveforms using simpler one-sided definition (b/c symmetry)
close all
H = zeros(length(Xs),1);
h=Ys;
Z = @(k) exp(-(1i)*2*k*pi/N);
freqs = Z(0:127);

for v=0:(N-1)
    sum = 0;
    for k=0:(N-1)
        sum = sum + h(k+1)*Z(v*k);
    end
    H(v+1) = sum;
end

plot(Xs,abs(H));
title('Frequency Response (Energy at each Frequency)');

input('Hit enter to continue');

% now we plot the graph with our harmonics and coefficients
Xs=0:(1/N):1;
Ys=zeros(length(Xs),1);
for k=0:(N-1)
    sum = 0;
    for v=0:(N-1)
        sum = sum + H(v+1) * Z(-k*v);
    end
    Ys(k+1) = sum/N;
end

plot(Xs,abs(Ys));
title('reconstruction with more samples');
%% Audio-fun Section!!
load('audio_data.mat');
soundsc(v1,fs);
%populates v1 (sound sample) and fs(sample-rate)
