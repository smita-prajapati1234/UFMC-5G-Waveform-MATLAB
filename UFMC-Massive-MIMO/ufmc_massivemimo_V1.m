 clc 
  clear all
%  close all
 s = rng(211);       % Set RNG state for repeatability
numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB
prototypeFilter = chebwin(filterLen, slobeAtten);
%prototypeFilter = kaiser(filterLen, 2.5);

%prototypeFilter = hamming (filterLen);

bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
map=2^bitsPerSubCarrier;
% QAM Symbol mapper
% qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
%     2^bitsPerSubCarrier, 'BitInput', true, ...
%     'NormalizationMethod', 'Average power');

% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
txSig = complex(zeros(numFFT+filterLen-1, 1));
           % SNR in dB
snrdB = 1:2:32;

for bandIdx = 1:numSubbands

    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);
    symbolsIn = qammod(bitsIn, map,'InputType', 'bit', 'UnitAveragePower', true);
    inpData(:,bandIdx) = bitsIn; % log bits for comparison
    
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize; 
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-subbandSize, 1)];
    ifftOut = ifft(ifftshift(symbolsInOFDM));
    
    % Filter for each subband is shifted in frequency
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                 ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );    
    filterOut = conv(bandFilter,ifftOut);
    
    % Plot power spectral density (PSD) per subband
%     [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
%                           numFFT*2, 1, 'centered'); 
%     plot(f,10*log10(psd)); 
%     
    % Sum the filtered subband responses to form the aggregate transmit
    % signal
    txSig = txSig + filterOut;     
end


% SETTING THE PARAMETERS FOR THE SIMULATION
fc=1e9;        %Carrier frequency
c=3e8;        %Speed of light
l=c/fc;        %Wavelength
d=l/2;        %Rx array spacing
N=64;         %Receive array size
M=16;         %Transmit array size (users)
theta=2*pi*(rand(1,M));         %Angular separation of users

%
% n=1:N(loop);                          %Rx array number
% n=transpose(n);                 %Row vector to column vector
sigma=0.3;

% RECEIVE SIGNAL MODEL (LINEAR)

H=(1/sqrt(2))*randn(N,M)+i*(1/sqrt(2))*randn(N,M);

w111=reshape([txSig;zeros(6,1)],M,35);
x=H*w111;                             %Receive vector of length N
% Add WGN
for loop=1:length(snrdB)

rxSig = awgn(x, snrdB(loop), 'measured');

% LINEAR ARRAY PROCESSING - METHODS
% 1-MATCHED FILTER
% y=pinv(H)*rxSig;
 %y=H'*rxSig;
y=pinv(H'*H)*H'*rxSig;

% 2-PINV without tol
  %y=pinv(H)*rxSig;

% 3-PINV with tol
%y=pinv(H,0.1)*rxSig;

% 4-Minimum Mean Square Error (MMMSE)
%y=(H'*H+(2*sigma^2)*eye([M,M]))^(-1)*(H')*rxSig;

y1=y(:);
 y1(end-5:end)=[];

yRxPadded = [y1; zeros(2*numFFT-numel(txSig),1)];

% Perform FFT and downsample by 2
RxSymbols2x = fftshift(fft(yRxPadded));
RxSymbols = RxSymbols2x(1:2:end);

% Select data subcarriers
dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));

% Plot received symbols constellation
% constDiagRx = comm.ConstellationDiagram('ShowReferenceConstellation', ...
%     false, 'Position', figposition([20 15 25 30]), ...
%     'Title', 'UFMC Pre-Equalization Symbols', ...
%     'Name', 'UFMC Reception', ...
%     'XLimits', [-150 150], 'YLimits', [-150 150]);
% constDiagRx(dataRxSymbols);

% Use zero-forcing equalizer after OFDM demodulation
rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
       zeros(numFFT-filterLen,1)];
prototypeFilterFreq = fftshift(fft(rxf));
prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));

% Equalize per subband - undo the filter distortion
dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
EqualizedRxSymbols = EqualizedRxSymbolsMat(:);

% Plot equalized symbols constellation
% constDiagEq = comm.ConstellationDiagram('ShowReferenceConstellation', ...
%     false, 'Position', figposition([46 15 25 30]), ...
%     'Title', 'UFMC Equalized Symbols', ...
%     'Name', 'UFMC Equalization');
% constDiagEq(EqualizedRxSymbols);

% Demapping and BER computation

BER = comm.ErrorRate;

% Perform hard decision and measure errors
rxBits = qamdemod(EqualizedRxSymbols, map, ...
    'OutputType', 'bit', ...         % Equivalent to 'BitOutput', true
    'UnitAveragePower', true);;
ber = BER(inpData(:), rxBits)
ber1(loop)=ber(1);
disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdB(loop)) ' dB']);

snr_linear = 10^(snrdB(loop)/10);
capacity(loop) = real(log2(det(eye(N) + snr_linear/M * (H * H'))));
% capacity(loop) = log2(det(eye(N) + snrdB(loop)/M * (H * H')));

end

% Restore RNG state
% rng(s);
figure (100)
hold on
semilogy(snrdB,smooth(ber1),'g')
set(gca,'Yscale','log')
grid on
xlabel('SNR')
ylabel('BER')

%  Plot capacity versus SNR
figure (101);
hold on
plot(snrdB, capacity, 'y-o', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Capacity (bits/symbol)');
