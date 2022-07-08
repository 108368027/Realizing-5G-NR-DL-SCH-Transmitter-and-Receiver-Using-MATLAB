clear all;

rv = 0;
NFrames = 1;
SNRdB = 500;                % SNR in dB
perfectEstimation = true; % Perfect synchronization and channel estimation
rng("default");            % Set default random number generator for repeatability

Carrier = nrCarrierConfig;         % Carrier resource grid configuration
Carrier.NSizeGrid = 132;            % Bandwidth in number of resource blocks (51 RBs at 30 kHz SCS for 20 MHz BW)
Carrier.SubcarrierSpacing = 60;    % 15, 30, 60, 120 (kHz)
Carrier.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
%Carrier.NCellID = 1;               % Cell identity
carrier = Carrier;

PDSCH = nrPDSCHConfig;
PDSCH.Modulation = "64QAM";
PDSCH.NumLayers = 2;
PDSCH.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation
% Scrambling identifiers
PDSCH.NID = Carrier.NCellID;
PDSCH.RNTI = 1;


% DM-RS and antenna port configuration (TS 38.211 Section 7.4.1.1)
PDSCH.DMRS.DMRSPortSet = 0:PDSCH.NumLayers-1; % DM-RS ports to use for the layers
PDSCH.DMRS.DMRSTypeAPosition = 2;      % Mapping type A only. First DM-RS symbol position (2,3)
PDSCH.DMRS.DMRSLength = 1;             % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
PDSCH.DMRS.DMRSAdditionalPosition = 2; % Additional DM-RS symbol positions (max range 0...3)
PDSCH.DMRS.DMRSConfigurationType = 2;  % DM-RS configuration type (1,2)
PDSCH.DMRS.NumCDMGroupsWithoutData = 3;% Number of CDM groups without data
PDSCH.DMRS.NIDNSCID = 1;               % Scrambling identity (0...65535)
PDSCH.DMRS.NSCID = 0;                  % Scrambling initialization (0,1)

% PT-RS configuration (TS 38.211 Section 7.4.1.2)
simParameters.PDSCH.EnablePTRS = 0;                  % Enable or disable PT-RS (1 or 0)

pdsch = PDSCH;

% Coding rate
if PDSCH.NumCodewords == 1
    codeRate = 517/1024;
else
    codeRate = [517 517]./1024;
end

% Create DL-SCH encoder object
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = false;
encodeDLSCH.TargetCodeRate = codeRate;

% Create DLSCH decoder object
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = false;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;

%
nTxAnts = 8;
nRxAnts = 8;

% Check that the number of layers is valid for the number of antennas
if pdsch.NumLayers > min(nTxAnts,nRxAnts)
    error("The number of layers ("+string(pdsch.NumLayers)+") must be smaller than min(nTxAnts,nRxAnts) ("+string(min(nTxAnts,nRxAnts))+")")
end

channel = nrCDLChannel;
channel.DelayProfile = "CDL-C";
[channel.TransmitAntennaArray.Size,channel.ReceiveAntennaArray.Size] = ...
        hArrayGeometry(nTxAnts,nRxAnts);

ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;

% Total number of slots in the simulation period
NSlots = NFrames * carrier.SlotsPerFrame;


%Set up a loop to simulate the transmission and reception of slots. Create a comm.ConstellationDiagram to display the constellation of the equalized signal.
constPlot = comm.ConstellationDiagram;                                          % Constellation diagram object
constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation); % Reference constellation values
constPlot.EnableMeasurements = 1;                                               % Enable EVM measurements

% Initial timing offset
offset = 0;

estChannelGrid = getInitialChannelEstimate(channel,carrier);
newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

for nSlot = 0:NSlots-1
    % New slot
    carrier.NSlot = nSlot;

    % Generate PDSCH indices info, which is needed to calculate the transport
    % block size
    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

    % Calculate transport block sizes
    Xoh_PDSCH = 0;
    %trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);
    trBlkSizes = 50000;
    for cwIdx = 1:pdsch.NumCodewords
        % Create and store a new transport block for transmission
            trBlk = transpose(repelem([1],trBlkSizes(cwIdx)));
            %trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1);
    end

    %Encode the transport blocks. The transport block is already stored in one of the internal soft buffers of the DL-SCH encoder object.
    codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G,rv);

    %Generate PDSCH symbols from the coded transport blocks.
    pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlock);

    %Get the precoding weights.
    precodingWeights = newPrecodingWeight;

    %Precode the PDSCH symbols.
    pdschSymbolsPrecoded = pdschSymbols*precodingWeights;

    %Generate DM-RS symbols and indices.
    dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
    
    %Generate an empty resource grid. This grid represents a slot.
    pdschGrid = nrResourceGrid(carrier,nTxAnts);

    [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
    pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;

    % PDSCH DM-RS precoding and mapping
    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
        dmrsSymbolsPrecoded = dmrsSymbols(:,p)*precodingWeights(p,:);
        pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbolsPrecoded;
    end
    ofdmInfo = nrOFDMInfo(carrier);
    [txWaveform1,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);
%%%%Propagation Channel%%%%
    %Pad the input signal with enough zeros to ensure that the generated signal is flushed out of the channel filter.
    chInfo = info(channel);
    maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
    %maxChDelay = 9;
    txWaveform = [txWaveform1; zeros(maxChDelay,size(txWaveform1,2))];
    %Send the signal through the channel and add noise.
    [rxWaveform1,pathGains] = channel(txWaveform);
    pathGains_size = size(pathGains);
    %pathGains: Channel path gains of the fading process
    noise = generateAWGN(SNRdB,nRxAnts,waveformInfo.Nfft,size(rxWaveform1));
    rxWaveform2 = rxWaveform1 + noise;
    %Perform perfect or practical timing estimation and synchronization.
%%%%Timing Synchronization%%%%
    if perfectEstimation
        % Get path filters for perfect timing estimation
        % Get path filter impulse response for link-level MIMO fading
        % channel (獲取 MIMO 衰落信道的路徑濾波器脈衝響應)
        pathFilters = getPathFilters(channel);
        %pathFilters : Path filter impulse response

        [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters); 
        %offset: Timing offset in samples
        %mag   : Channel impulse response magnitude
        
    else
        [t,mag] = nrTimingEstimate(carrier,rxWaveform2,dmrsIndices,dmrsSymbols);
        offset = hSkipWeakTimingOffset(offset,t,mag);
    end
    rxWaveform = rxWaveform2(1+offset:end,:);

    %OFDM-demodulate the synchronized signal.
    rxGrid = nrOFDMDemodulate(carrier,rxWaveform);

    %Perform perfect or practical channel estimation.
    if perfectEstimation
        % Perform perfect channel estimation between transmit and receive
        % antennas.
        estChGridAnts = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset);
        estChGridAnts1 = estChGridAnts(:,:,1,1);
        estChGridAnts_size = size(estChGridAnts);
        % Get perfect noise estimate (from noise realization)
        noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end ,:));
        noiseGrid_size = size(noiseGrid);
        noiseEst = var(noiseGrid(:));

        % Get precoding matrix for next slot
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);

        % Apply precoding to estChGridAnts. The resulting estimate is for
        % the channel estimate between layers and receive antennas.
        [estChGridLayers,estChannelGrid_size] = precodeChannelEstimate(estChGridAnts,precodingWeights.');
        estChGridLayers_size = size(estChGridLayers);
    else
        % Perform practical channel estimation between layers and receive
        % antennas.
        [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);

        % Remove precoding from estChannelGrid before precoding
        % matrix calculation
        estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));

        % Get precoding matrix for next slot
        newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
    end

    %%繪製第一層和第一個接收天線之間的信道估計。
    mesh(abs(estChGridLayers(:,:,1,1)));
    title('Channel Estimate');
    xlabel('OFDM Symbol');
    ylabel("Subcarrier");
    zlabel("Magnitude");


    %The equalizer uses the channel estimate to compensate for the distortion introduced by the channel.

    %Extract the PDSCH symbols from the received grid and associated channel estimates. 
    %The csi output has channel state information (CSI) for each of the equalized PDSCH symbols. 
    %The CSI is a measure of the channel conditions for each PDSCH symbol. 
    %Use the CSI to weight the decoded soft bits after PDSCH decoding,
    %effectively increasing the importance of symbols experiencing better channel conditions.
    [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
    total = real(pdschRx).^2+imag(pdschRx).^2;
    pdschRx_mean = mean(sum(total,2));
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
    csi1 = csi;
    
    %Plot the constellation of the equalized symbols. The plot includes the constellation diagrams for all layers.
    constPlot.ChannelNames = "Layer "+(pdsch.NumLayers:-1:1);
    constPlot.ShowLegend = true;
    % Constellation for the first layer has a higher SNR than that for the
    % last layer. Flip the layers so that the constellations do not mask
    % each other.
    constPlot(fliplr(pdschEq));


    %PDSCH Decoding
    [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);
    dlschLLRs1 = dlschLLRs;
    % Scale LLRs by CSI
    csi = nrLayerDemap(csi);                                    % CSI layer demapping
    csi2 = csi;
    for cwIdx = 1:pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % Bits per symbol
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % Expand by each bit per symbol
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % Scale
    end

    %DL-SCH Decoding
    decodeDLSCH.TransportBlockLength = trBlkSizes;
    [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,rv);

    disp("Slot "+(nSlot));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noise = generateAWGN(SNRdB,nRxAnts,Nfft,sizeRxWaveform)
% Generate AWGN for a given value of SNR in dB (SNRDB), which is the
% receiver SNR per RE and antenna, assuming the channel does
% not affect the power of the signal. NRXANTS is the number of receive
% antennas. NFFT is the FFT size used in OFDM demodulation. SIZERXWAVEFORM
% is the size of the receive waveform used to calculate the size of the
% noise matrix.

    % Normalize noise power by the IFFT size used in OFDM modulation, as
    % the OFDM modulator applies this normalization to the transmitted
    % waveform. Also normalize by the number of receive antennas, as the
    % channel model applies this normalization to the received waveform by
    % default. The SNR is defined per RE for each receive antenna (TS
    % 38.101-4).
    SNR = 10^(SNRdB/10); % Calculate linear noise gain
    N0 = 1/sqrt(2.0*nRxAnts*double(Nfft)*SNR);
    noise = N0*complex(randn(sizeRxWaveform),randn(sizeRxWaveform));
end
    
function wtx = getPrecodingMatrix(PRBSet,NLayers,hestGrid)
% Calculate precoding matrix given an allocation and a channel estimate
    
    % Allocated subcarrier indices
    allocSc = (1:12)' + 12*PRBSet(:).';
    allocSc = allocSc(:);
    
    % Average channel estimate
    [~,~,R,P] = size(hestGrid);
    estAllocGrid = hestGrid(allocSc,:,:,:);
    Hest = permute(mean(reshape(estAllocGrid,[],R,P)),[2 3 1]);
    
    % SVD decomposition
    [~,~,V] = svd(Hest);
    
    wtx = V(:,1:NLayers).';
    wtx = wtx/sqrt(NLayers); % Normalize by NLayers
end

function estChannelGrid = getInitialChannelEstimate(channel,carrier)
% Obtain an initial channel estimate for calculating the precoding matrix.
% This function assumes a perfect channel estimate

    % Clone of the channel
    chClone = channel.clone();
    chClone.release();

    % No filtering needed to get channel path gains
    chClone.ChannelFiltering = false;    
    
    % Get channel path gains
    [pathGains,sampleTimes] = chClone();
    
    % Perfect timing synchronization
    pathFilters = getPathFilters(chClone);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
end

function refPoints = getConstellationRefPoints(mod)
% Calculate the reference constellation points for a given modulation
% scheme.
    switch mod
        case "QPSK"
            nPts = 4;
        case "16QAM"
            nPts = 16;
        case "64QAM"
            nPts = 64;
        case "256QAM"
            nPts = 256;            
    end
    binaryValues = int2bit(0:nPts-1,log2(nPts));
    refPoints = nrSymbolModulate(binaryValues(:),mod);
end

function [estChannelGrid,estChannelGrid_size] = precodeChannelEstimate(estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate.

    % Linearize 4-D matrix and reshape after multiplication
    K = size(estChannelGrid,1);
    L = size(estChannelGrid,2);
    R = size(estChannelGrid,3);
    estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
    estChannelGrid_size = size(estChannelGrid);
    estChannelGrid = estChannelGrid*W;
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end
