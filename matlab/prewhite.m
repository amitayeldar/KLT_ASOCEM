function [noiseMcPreWhite] = prewhite(noiseMc,apprxNoisePsd,apprxCleanPsd,mcSz,patchSz,R)
% prewhite the micrograph using the noise RPSD
%  
%  Input:
%   noiseMc            the micrograph we want to prewhiten.
%   apprxNoisePsd      RPSD of the noise
%   apprxCleanPsd      RPSD of the particle
%   patchSz            The patch size from whom we estimate the noise RPSD.
%   R                  The radii samples of the RPSD
%  Output:
%   noiseMcPreWhite    the prewhitened micrograph           
% 
%
% Amitay Eldar November 2019.

r = floor((mcSz(2)-1)/2);
c = floor((mcSz(1)-1)/2);
T = 1; % Nyquist sampling rate
bandLimit = pi/T;
col = (-c:1:c)*(bandLimit/c);
row = (-r:1:r)*(bandLimit/r);
[Row,Col] = meshgrid(row,col);
radMat = sqrt(Col.^2+Row.^2);
[radSamp,~,idx] = unique(radMat);
radSampTmp(radSamp<=max(R*bandLimit)) = radSamp(radSamp<=max(R*bandLimit));
noisePsdNodes = abs(triginterp(radSampTmp,R*bandLimit,apprxNoisePsd));
noisePsdNodesExtend = zeros(size(radSamp,1),1);
noisePsdNodesExtend(1:size(noisePsdNodes,2))= noisePsdNodes;
noisePsdNodesExtend(size(noisePsdNodes,2):end)= noisePsdNodes(end);
noisePsdMat = reshape(noisePsdNodesExtend(idx),length(col),length(row)); % This mat will use for whitening.
noiseMcPreWhite = cryo_prewhiten(noiseMc , noisePsdMat); % we dont want zeros for whitening apprxCleanPsd.

end