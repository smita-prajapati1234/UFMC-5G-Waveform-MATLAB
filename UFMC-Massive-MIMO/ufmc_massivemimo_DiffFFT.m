% clc 
% clear all
% close all
s = rng(211);       % Set RNG state for repeatability
numFFT = [512];        % number of FFT points
subbandSize = [30] ;    % must be > 1 
numSubbands = [10];    % numSubbands*subbandSize <= numFFT

clr={'b','r','k','g'}
p=[67]% 67 131 131]%[35 35 35]
% Dolph-Chebyshev window design parameters
filterLen = 83;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB
prototypeFilter = chebwin(filterLen, slobeAtten);
%for sloop= size(filterLen)%
%prototypeFilter = kaiser(filterLen, 2.5);
%prototypeFilter = chebwin(filterLen(sloop), slobeAtten);

bitsPerSubCarrier = 8;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');
for sloop=1:1
    
    subbandOffset = (numFFT(sloop)/2)-((subbandSize(sloop)* numSubbands(sloop))/2 )    %subbandSize*numSubbands/2 for band center
% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize(sloop), numSubbands(sloop));
txSig = complex(zeros(numFFT(sloop)+filterLen-1, 1));
snrdB = 1:1:25;              % SNR in dB


hFig = figure;
axis([-0.5 0.5 -100 20]);
hold on; 
grid on

xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['UFMC, ' num2str(numSubbands(sloop)) ' Subbands, '  ...
    num2str(subbandSize(sloop)) ' Subcarriers each'])

%  Loop over each subband
for bandIdx = 1:numSubbands(sloop)

    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize(sloop), 1);
    symbolsIn = qamMapper(bitsIn);
    inpData(:,bandIdx) = bitsIn; % log bits for comparison
    
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize(sloop); 
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT(sloop)-offset-subbandSize(sloop), 1)];
    ifftOut = ifft(ifftshift(symbolsInOFDM));
    
    % Filter for each subband is shifted in frequency
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT(sloop)* ...
                 ((bandIdx-1/2)*subbandSize(sloop)+0.5+subbandOffset+numFFT(sloop)/2) );    
    filterOut = conv(bandFilter,ifftOut);
    
    % Plot power spectral density (PSD) per subband
    [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                          numFFT(sloop)*2, 1, 'centered'); 
   % subplot(2,2,sloop)
    plot(f,10*log10(psd)); 
    
    % Sum the filtered subband responses to form the aggregate transmit
    % signal
    txSig = txSig + filterOut;     
end
% set(hFig, 'Position', figposition([20 50 25 30]));
% hold off;

% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txSig);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);


% SETTING THE PARAMETERS FOR THE SIMULATION
fc=1e9;        %Carrier frequency
c=3e8;        %Speed of light
l=c/fc;        %Wavelength
d=l/2;        %Rx array spacing
N=20;         %Receive array size
M=16;         %Transmit array size (users)
theta=2*pi*(rand(1,M));         %Angular separation of users

%
% n=1:N(loop);                          %Rx array number
% n=transpose(n);                 %Row vector to column vector
sigma=0.1;

% RECEIVE SIGNAL MODEL (LINEAR)
%H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);  %Channel matrix of size NxM
%H=H/norm(H'*H);
H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);
%w111=reshape([txSig;zeros(6,1)],M,560);
w111=reshape([txSig;zeros(14,1)],M,70);
x=H*w111;                             %Receive vector of length N
% Add WGN
for loop=1:length(snrdB)

rxSig = awgn(x, snrdB(loop), 'measured');

% LINEAR ARRAY PROCESSING - METHODS
% 1-MATCHED FILTER
 %y=pinv(H)*rxSig;
 
y=pinv(H'*H)*H'*rxSig;

% 2-PINV without tol
%  y=pinv(H)*rxSig;

% 3-PINV with tol
%  y=pinv(H,0.1)*rxSig;

% 4-Minimum Mean Square Error (MMMSE)
%y=(H'*H+(2*sigma^2)*eye([M,M]))^(-1)*(H')*rxSig;

% H=exp(-i*(n-1)*2*pi*d*cos(theta)/l);
% H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);

y1=y(:);
 y1(end-13:end)=[];

yRxPadded = [y1; zeros(2*numFFT(sloop)-numel(txSig),1)];

% Perform FFT and downsample by 2
RxSymbols2x = fftshift(fft(yRxPadded));
RxSymbols = RxSymbols2x(1:2:end);

% Select data subcarriers
dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands(sloop)*subbandSize(sloop)));

% Plot received symbols constellation
% constDiagRx = comm.ConstellationDiagram('ShowReferenceConstellation', ...
%     false, 'Position', figposition([20 15 25 30]), ...
%     'Title', 'UFMC Pre-Equalization Symbols', ...
%     'Name', 'UFMC Reception', ...
%     'XLimits', [-150 150], 'YLimits', [-150 150]);
% constDiagRx(dataRxSymbols);

% Use zero-forcing equalizer after OFDM demodulation
rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT(sloop)); ...
       zeros(numFFT(sloop)-filterLen,1)];
prototypeFilterFreq = fftshift(fft(rxf));
prototypeFilterInv = 1./prototypeFilterFreq(numFFT(sloop)/2-subbandSize(sloop)/2+(1:subbandSize(sloop)));

% Equalize per subband - undo the filter distortion
dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize(sloop),numSubbands(sloop));
EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
EqualizedRxSymbols = EqualizedRxSymbolsMat(:);

% Plot equalized symbols constellation
% constDiagEq = comm.ConstellationDiagram('ShowReferenceConstellation', ...
%     false, 'Position', figposition([46 15 25 30]), ...
%     'Title', 'UFMC Equalized Symbols', ...
%     'Name', 'UFMC Equalization');
% constDiagEq(EqualizedRxSymbols);

% Demapping and BER computation
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
BER = comm.ErrorRate;

% Perform hard decision and measure errors
rxBits = qamDemod(EqualizedRxSymbols);
ber = BER(inpData(:), rxBits)
ber1(loop)=ber(1);
disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdB(loop)) ' dB']);
end

% Restore RNG state
rng(s);
figure (100)
hold on
semilogy(snrdB,smooth(ber1),'r')
set(gca,'Yscale','log')
grid on
xlabel('SNR')
ylabel('BER')

end
