function LBP= efficientLBP(inImg, varargin) % filtR, isRotInv, isChanWiseRot
%% efficientLBP
% The function implements LBP (Local Binary Pattern analysis).
%
%% Syntax
%  LBP= efficientLBP(inImg);
%
%% Description
% The LBP tests the relation between pixel and it's neighbors, encoding this relation into
%   a binary word. This allows detection of patterns/features.
% The function is inpired by materials published by Matti Pietik?inen in
%   http://www.cse.oulu.fi/CMV/Research/LBP . This implementation hovewer is not totally
%   allighned with the mthods proposed by Professor Pietik?inen (see Issues & Comments).
%
%% Input arguments (defaults exist):
% inImg- input image, a 2D matrix (3D color images will be converted to 2D intensity
%     value images)
% filtR- a 2D matrix representing a round/radial filter. It can be generated using
%   generateRadialFilterLBP function.
% isRotInv- a logical flag. When enabled generated rotation invariant LBP accuired via
%     fining an angle at whihc the LBP og a given pixelis minimal. Icreases run time, and
%     results in a relatively sparce hsitogram (as many combinations disappear).
% isChanWiseRot- a logical flag, when enabled (default value) allowes channel wise
%     rotation. When disabled/false rotation carried out based on roation of first color
%     channel. Supported only when "isEfficent" is enabled. When  "isEfficent" is
%     disabled "isChanWiseRot" is true.
%
%% Output arguments
%   LBP-    LBP image UINT8/UINT16/UINT32/UINT64/DOUBLE of same dimentions
%     [Height x Width] as inImg.
%
%% Issues & Comments
% - Currenlty, all neigbours are treated alike. Basically, we can use wighted/shaped
%     filter.
% - The rotation invariant LBP histogram includes less then bins then regular LBP BY
%     DEFINITION the zero trailing binary words are excluded for example, so it can be
%     reduced to a mush more component representation. Actually for 8 niegbours it's 37
%     bins, instead of 256. An efficnet way to calculate those bins value is needed.
%
%% Example
% img=imread('peppers.png');
% filtR=generateRadialFilterLBP(8, 1);
% tic;
%  % note this filter dimentions aren't legete...
% effLBP= efficientLBP(img, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false);
% effTime=toc;
%
% % verify pixel wise implementation returns same results
% tic;
% % same parameters as before
% pwLBP=pixelwiseLBP(img, 'filtR', filtR, 'isRotInv', true, 'isChanWiseRot', false); 
% inEffTime=toc;
% fprintf('\nRun time ratio %.2f. Same result eqaulity chesk: %o.\n', inEffTime/effTime,...
%    isequal(effLBP, pwLBP));
%
% figure;
% subplot(1, 3, 1)
% imshow(img);
% title('Original image');
%
% subplot(1, 3, 2)
% imshow( effLBP );
% title('Efficeint LBP image');
%
% subplot(1, 3, 3)
% imshow( pwLBP );
% title('Pixel-wise LBP image');
%
%% See also
% pixelwiseLBP  % a straigh forward iplmenetation of LBP, should achive same results
% generateRadialFilterLBP   % custom function generating circulat filters
%
%% Revision history
% First version: Nikolay S. 2012-05-01.
% Last update:   Nikolay S. 2014-01-09.
%
% 
%% Deafult params
isRotInv=false;
isChanWiseRot=false;
filtR=generateRadialFilterLBP(8, 1);

%% Get user inputs overriding default values
funcParamsNames={'filtR', 'isRotInv', 'isChanWiseRot'};
assignUserInputs(funcParamsNames, varargin{:});

if ischar(inImg) && exist(inImg, 'file')==2 % In case of file name input- read graphical file
    inImg=imread(inImg);
end

nClrChans=size(inImg, 3);

inImgType=class(inImg);
calcClass='single';

isCalcClassInput=strcmpi(inImgType, calcClass);
if ~isCalcClassInput
    inImg=cast(inImg, calcClass);
end
imgSize=size(inImg);

nNeigh=size(filtR, 3);

if nNeigh<=8
    outClass='uint8';
elseif nNeigh>8 && nNeigh<=16
    outClass='uint16';
elseif nNeigh>16 && nNeigh<=32
    outClass='uint32';
elseif nNeigh>32 && nNeigh<=64
    outClass='uint64';
else
    outClass=calcClass;
end

if isRotInv
    nRotLBP=nNeigh;
    nPixelsSingleChan=imgSize(1)*imgSize(2);
    iSingleChan=reshape( 1:nPixelsSingleChan, imgSize(1), imgSize(2) );
else
    nRotLBP=1;
end

nEps=-3;
weigthVec=reshape(2.^( (1:nNeigh) -1), 1, 1, nNeigh);
weigthMat=repmat( weigthVec, imgSize([1, 2]) );
binaryWord=zeros(imgSize(1), imgSize(2), nNeigh, calcClass);
LBP=zeros(imgSize, outClass);
possibleLBP=zeros(imgSize(1), imgSize(2), nRotLBP);
for iChan=1:nClrChans  
    % Initiate neighbours relation filter and LBP's matrix
    for iFiltElem=1:nNeigh
        % Rotate filter- to compare center to next neigbour
        filtNeight=filtR(:, :, iFiltElem);
        
        % calculate relevant LBP elements via filtering
        binaryWord(:, :, iFiltElem)=cast( ...
            roundnS(filter2( filtNeight, inImg(:, :, iChan), 'same' ), nEps) >= 0,...
            calcClass );
        % Without rounding sometimes inaqulity happens in some pixels
        % compared to pixelwiseLBP
    end % for iFiltElem=1:nNeigh

    for iRot=1:nRotLBP
        % find all relevant LBP candidates
        possibleLBP(:, :, iRot)=sum(binaryWord.*weigthMat, 3);
        if iRot < nRotLBP
            binaryWord=circshift(binaryWord, [0, 0, 1]); % shift binaryWord elements
        end
    end
    
    if isRotInv
        if iChan==1 || isChanWiseRot
            % Find minimal LBP, and the rotation applied to first color channel
            [minColroInvLBP, iMin]=min(possibleLBP, [], 3);
            
            % calculte 3D matrix index
            iCircShiftMinLBP=iSingleChan+(iMin-1)*nPixelsSingleChan;
        else
            % the above rotation of the first channel, holds to rest of the channels
            minColroInvLBP=possibleLBP(iCircShiftMinLBP);
        end % if iChan==1 || isChanWiseRot
    else
        minColroInvLBP=possibleLBP;
    end % if isRotInv
    
    if strcmpi(outClass, calcClass)
        LBP(:, :, iChan)=minColroInvLBP;
    else
        LBP(:, :, iChan)=cast(minColroInvLBP, outClass);
    end
end % for iChan=1:nClrChans