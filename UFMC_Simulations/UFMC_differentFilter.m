 clc 
 clear all
% close all
s = rng(211);       % Set RNG state for repeatability
numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1 
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center
%clr={'b','g','k','r','m','c','--gs','o','y'}
% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB
% prototypeFilter = chebwin(filterLen, slobeAtten);
%for sloop= size(filterLen)%
% prototypeFilter = kaiser(filterLen, 2.5);
%prototypeFilter = chebwin(filterLen(sloop), slobeAtten);
prototypeFilter = hamming (filterLen);
bitsPerSubCarrier = 8;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
txSig = complex(zeros(numFFT+filterLen-1, 1));
snrdB = 1:2:22;              % SNR in dB


hFig = figure;
axis([-0.5 0.5 -100 20]);
hold on; 
grid on

xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['UFMC, ' num2str(numSubbands) ' Subbands, '  ...
    num2str(subbandSize) ' Subcarriers each'])

%  Loop over each subband
for bandIdx = 1:numSubbands

    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);
    symbolsIn = qamMapper(bitsIn);
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
    [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                          numFFT*2, 1, 'centered'); 
    plot(f,10*log10(psd)); 
    
    % Sum the filtered subband responses to form the aggregate transmit
    % signal
    txSig = txSig + filterOut;     
end
set(hFig, 'Position', figposition([20 50 25 30]));
hold off;

% Compute peak-to-average-power ratio (PAPR)
PAPR = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprUFMC] = PAPR(txSig);
disp(['Peak-to-Average-Power-Ratio (PAPR) for UFMC = ' num2str(paprUFMC) ' dB']);


for loop=1:length(snrdB)
% Add WGN
rxSig = awgn(txSig, snrdB(loop), 'measured');

%%

yRxPadded = [rxSig; zeros(2*numFFT-numel(txSig),1)];

% Perform FFT and downsample by 2
RxSymbols2x = fftshift(fft(yRxPadded));
RxSymbols = RxSymbols2x(1:2:end);

% Select data subcarriers
dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));

% Plot received symbols constellation
constDiagRx = comm.ConstellationDiagram('ShowReferenceConstellation', ...
    false, 'Position', figposition([20 15 25 30]), ...
    'Title', 'UFMC Pre-Equalization Symbols', ...
    'Name', 'UFMC Reception', ...
    'XLimits', [-150 150], 'YLimits', [-150 150]);
constDiagRx(dataRxSymbols);

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
constDiagEq = comm.ConstellationDiagram('ShowReferenceConstellation', ...
    false, 'Position', figposition([46 15 25 30]), ...
    'Title', 'UFMC Equalized Symbols', ...
    'Name', 'UFMC Equalization');
constDiagEq(EqualizedRxSymbols);

% Demapping and BER computation
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
BER = comm.ErrorRate;

% Perform hard decision and measure errors
rxBits = qamDemod(EqualizedRxSymbols);
ber = BER(inpData(:), rxBits)
ber1(loop)=ber(1);
disp(['UFMC Reception, BER = ' num2str(ber1(loop)) ' at SNR = ' ...
    num2str(snrdB(loop)) ' dB']);
end
% Restore RNG state
rng(s);
figure (100)
hold on
semilogy(snrdB,smooth(ber1),'g')
set(gca,'Yscale','log')
grid on
xlabel('SNR')
ylabel('BER')
title('BEr Vs SNR for UFMC at 256 QAM for different filters')
