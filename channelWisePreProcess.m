
function out = channelWisePreProcess(in)
% As input has 4 channels (modalities), remove the mean and divide by the
% standard deviation of each modality independently.
chn_Mean = mean(in,[1 2 ]);
chn_Std = std(in,0,[1 2 ]);
out = (in - chn_Mean)./chn_Std;

end