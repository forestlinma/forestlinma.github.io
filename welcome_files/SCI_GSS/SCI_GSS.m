function [score, SimMatrix] = SCI_GSS(Im1, Im2)
% ========================================================================
% GSS Index for screen content image, Version 1.0
% Copyright(c) 2016 Zhangkai Ni, Lin Ma, Huanqiang Zeng, Canhui Cai, and 
% Kai-Kuang Ma
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for calculating the
% gradient similarity score (GSS) index between two images.
%
% For more information, please refer to the following paper:
%
% Zhangkai Ni, Lin Ma, Huanqiang Zeng, Canhui Cai, and Kai-Kuang Ma,"Gradient 
% Direction for Screen Content Image Quality Assessment", IEEE Signal Processing 
% Letters
% 
%----------------------------------------------------------------------
%
% Input : (1) Im1: the first image being compared (grayscale image, double type, 0~255)
%         (2) Im2: the second image being compared (grayscale image, double type, 0~255)
%
% Output: (1) GSS: is the similarty score calculated using GSS algorithm. 
%	          GSS only considers the luminance component of images. 
%
%         (2) SimMatrix: is the local quality map of the distorted image
%        
%-----------------------------------------------------------------------
%
% Usage:
% Given two test images Im1 and Im2. For gray-scale images, their dynamic 
% range should be 0-255. For colorful images, the dynamic range of each 
% color channel should be 0-255.
%
%        [score, SIMap] = GSS(Im1, Im2);
%-----------------------------------------------------------------------
%% parameters 
dirNum = 12;     % the number of directions 
ks = 13;         % the length of convolution line 

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Downsample the image
%%%%%%%%%%%%%%%%%%%%%%%%%
Down_step = 2;
aveKernel = fspecial('average',2);
aveIm1 = conv2(Im1, aveKernel,'same');
aveIm2 = conv2(Im2, aveKernel,'same');
Im1 = aveIm1(1:Down_step:end,1:Down_step:end);
Im2 = aveIm2(1:Down_step:end,1:Down_step:end);

%%  image gradient magnitude
[H, W] = size(Im1);
Im1X = [abs(Im1(:,1:(end-1)) - Im1(:,2:end)),zeros(H,1)];
Im1Y = [abs(Im1(1:(end-1),:) - Im1(2:end,:));zeros(1,W)];
Im1Edge = Im1X + Im1Y;

[H, W] = size(Im2);
Im2X = [abs(Im2(:,1:(end-1)) - Im2(:,2:end)),zeros(H,1)];
Im2Y = [abs(Im2(1:(end-1),:) - Im2(2:end,:));zeros(1,W)];
Im2Edge = Im2X + Im2Y;

%% convolution kernel with horizontal direction
kerRef = zeros(ks*2+1);
kerRef(ks+1,:) = 1;

%% image gradient direction
response1 = zeros(H,W,dirNum);
for ii = 0 : (dirNum-1)
    ker = imrotate(kerRef, ii*180/dirNum, 'bilinear', 'crop');
    response1(:,:,ii+1) = conv2(Im1Edge, ker, 'same');
end
[~ , index] = sort(response1, 3);
Dirs1 = (index(:,:,end)-1)*180/dirNum;

response2 = zeros(H,W,dirNum);
for ii = 0 : (dirNum-1)
    ker = imrotate(kerRef, ii*180/dirNum, 'bilinear', 'crop');
    response2(:,:,ii+1) = conv2(Im2Edge, ker, 'same');
end
[~ , index] = sort(response2, 3);
Dirs2 = (index(:,:,end)-1)*180/dirNum;

T1 = 205;
T2 = 160;         
MagnitudeSimMatrix = (2*Im1Edge .* Im2Edge + T2) ./(Im1Edge.^2 + Im2Edge.^2 + T2);
DirectionSimMatrix = (2*Dirs1 .* Dirs2 + T1) ./(Dirs1.^2 + Dirs2.^2 + T1);
SimMatrix = DirectionSimMatrix .* MagnitudeSimMatrix;
score = std2(SimMatrix);
