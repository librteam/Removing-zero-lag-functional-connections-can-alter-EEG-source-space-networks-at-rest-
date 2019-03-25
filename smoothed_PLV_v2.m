function plv=smoothed_PLV_v2(Signal,srate,window)

%% Inputs/Outputs

numChannels = size(Signal, 1);
numSamples= size(Signal,2);
winSamples= window * srate; 
N=floor(numSamples/winSamples);
plv = zeros(N,numChannels, numChannels);

%% Extract the instantaneous phase
for channelCount = 1:numChannels
    Signal(channelCount, :) = angle(hilbert((Signal(channelCount, :))));
end

%% Compute PLV

No=winSamples/2+1;

for count=1:N
    for channelCount = 1:numChannels-1
        channelData = squeeze(Signal(channelCount, No-winSamples/2:No+winSamples/2-1));
        for compareChannelCount = channelCount+1:numChannels
            compareChannelData = squeeze(Signal(compareChannelCount, No-winSamples/2:No+winSamples/2-1));
                diff=channelData(:, :) - compareChannelData(:, :);
                diff=diff';
               plv(count,channelCount,compareChannelCount) =abs(sum(exp(i*diff)))/length(diff);  
        end
    end
    
    No=No+winSamples;
end