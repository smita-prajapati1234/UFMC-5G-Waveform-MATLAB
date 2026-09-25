% clc;
% clear all;
% close all;


s = rng(211);       % Set RNG state for repeatability
numFFT = 512;        % number of FFT points
subbandSize = 40;    % must be > 1
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 56; % numFFT/2-subbandSize*numSubbands/2 for band center

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB

bitsPerSubCarrier = 8;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM

snrdB = 15;              % SNR in Db
numRuns = 10000;
PAPRarray=[];
                                                  % Design window with specified attenuation
prototypeFilter = chebwin(filterLen, slobeAtten);

  clr={'b','k','r'}                                                                     % QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');
%for sloop=1:3
% Transmit-end processing
%  Initialize arrays
for n = 1:numRuns
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
txSig = complex(zeros(numFFT+filterLen-1, 1));

% hFig = figure;
% axis([-0.5 0.5 -100 20]);
% hold on;
% grid on

% xlabel('Normalized frequency');
% ylabel('PSD (dBW/Hz)')
% title(['UFMC, ' num2str(numSubbands) ' Subbands, '  ...
%     num2str(subbandSize) ' Subcarriers each'])

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
    %plot(f,10*log10(psd));

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

PAPRarray(n)=paprUFMC;
end

[N,X] = histcounts(PAPRarray, 200,'Normalization','count');
figure (100)
hold on
 semilogy(X(1:end-1),1-cumsum(N)/max(cumsum(N)),'r')
  grid on
 xlabel('PAPR')
 ylabel('CCDF')
 title('UFMC PAPR')
