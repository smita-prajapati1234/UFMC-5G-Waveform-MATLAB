%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MASSIVE MIMO BEAMFORMING
% COPYRIGHT RAYMAPS (C) 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
rng(111)
% SETTING THE PARAMETERS FOR THE SIMULATION
f=1e9;        %Carrier frequency
c=3e8;        %Speed of light
l=c/f;        %Wavelength
d=l/2;        %Rx array spacing
N=20;         %Receive array size
M=16;         %Transmit array size (users)
theta=2*pi*(rand(1,M));         %Angular separation of users
EbNo=0;                        %Energy per bit to noise PSD
sigma=1/sqrt(2*EbNo);           %Standard deviation of noise

n=1:N;                          %Rx array number
n=transpose(n);                 %Row vector to column vector

% RECEIVE SIGNAL MODEL (LINEAR)
s=2*(round(rand(M,1))-0.5);           %BPSK signal of length M
H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);  %Channel matrix of size NxM

H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);
wn=sigma*(randn(N,1)+i*randn(N,1));   %AWGN noise of length N
x=H*s+wn;                             %Receive vector of length N

% LINEAR ARRAY PROCESSING - METHODS
% 1-MATCHED FILTER
 y=H'*(x);%/sqrt(norm((H*H')));

% 2-PINV without tol
%   y=pinv(H)*x;

% 3-PINV with tol
%  y=pinv(H,0.1)*x;

% 4-Minimum Mean Square Error (MMMSE)
%  y=(H'*H+(2*sigma^2)*eye([M,M]))^(-1)*(H')*x;

% DEMODULATION AND BER CALCULATION
s_est=sign(real(y));          %Demodulation
ber=sum(s ~= s_est)/length(s)  %BER calculation

%Note: Please select the array processing technique
%you want to implement (1-MF, 2-LS1, 3-LS2, 4-MMSE)
return


H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);

H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);