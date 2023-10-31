function [regMovie, x, y, dx, dy, target] = rapidRegNonPar(movie, target, varargin)
% rapidReg Movie frame registration to a target frame
%
%   [REG, X, Y, DX, DY, TARGET] = rapidReg(MOVIE, TARGET,...) register 
%   MOVIE (an (X,Y,T) array) to TARGET (either the target image frame (X,Y)
%   or 'auto', to find one automatically). Returns the registered movie, REG,
%   the x and y coordinates of the registered movies pixels relative to the
%   target (useful for when clipping is used), the translations required to
%   register each frame, DX and DY, and TAGRET, the target frame. 

[h, w, nFrames] = size(movie);

%% Setup
%convert movie data to an appropriate floating point type if necessary
dataType = class(movie);
switch dataType
  case {'int8' 'int16' 'uint8' 'uint16' 'int32' 'uint32'}
    %convert all integer types up to 32-bits to single
    movie = single(movie);
    origTypeFun = str2func(dataType);
  case {'int64'  'uint64'}
    %convert 64-bit integer types to double
    movie = double(movie);
    origTypeFun = str2func(dataType);
  case {'single' 'double'}
    %no conversion
    origTypeFun = @identity;
  otherwise
    error('''%s'' is not a recognised data type', dataType);
end

%create a Gaussian filter for filtering frames
hGauss = fspecial('gaussian', [5 5], 3);
%look for flag on whether to use parallel codepath
% if any(cell2mat(strfind(varargin, 'nopar')) == 1)
%   parallel = false;
% else
%   parallel = true;
% end
parallel = false;
if length(varargin)>0
    for iArg = 1:2:length(varargin)
        argNames(ceil(iArg/2)) = varargin(iArg);
        argVales(ceil(iArg/2)) = varargin(iArg+1);
    end
else
    argNames = [];
    argVales = [];
end

if ismember('filter',lower(argNames))
    [~,idx] = ismember('filter',argNames);
    pair = cell2mat(argVales(idx));
    hGauss = pair;
end

if strcmpi(target, 'auto')
  %% Compute the best target frame
  %first compute a smoothed mean of each frame
  meanF = smooth(mean(reshape(movie, h*w, nFrames)));
  %now look in the middle third of the image frames for the minimum
  fromFrame = round(nFrames*1/3);
  toFrame = round(nFrames*2/3);
  [~, idx] = min(meanF(fromFrame:toFrame));
  minFrame = fromFrame + idx;
  %Gaussian filter the target image
  target = imfilter(movie(:,:,minFrame), hGauss, 'same', 'replicate');
  % AR added display output of which frame number is being used as target
  disp(['target frame is ',num2str(minFrame),'..']);
end

ftTarget = fft2(target);

%% Fourier transform the movie frames, unfiltered and filtered
ftMovie = fft2(movie);
movie = imfilter(movie, hGauss, 'same', 'replicate');

%% Compute required displacement and register each frame
dx = zeros(1, nFrames);
dy = zeros(1, nFrames);
nr = ifftshift((-fix(h/2):ceil(h/2) - 1));
nc = ifftshift((-fix(w/2):ceil(w/2) - 1));
[nc, nr] = meshgrid(nc, nr);
regMovie = zeros(h, w, nFrames, class(movie));


  %% Register sequentially
  for t = 1:nFrames
    %find the best registration translation
    currentFrame = fft2(movie(:,:,t));
    output = dftregistration(fft2(target), currentFrame, 20);
    dx(t) = output(4);
    dy(t) = output(3);
    %translate the original (i.e. unfiltered) frame
    ftRegFrame = ftMovie(:,:,t).*exp(sqrt(-1)*2*pi*(-dy(t)*nr/h - dx(t)*nc/w));
    regMovie(:,:,t) = abs(ifft2(ftRegFrame));
  end

% %% If requested, clip the frames to the maximum fully valid region
% if any(cell2mat(strfind(varargin, 'clip')) == 1)
%   disp('clipping..');
%   dxMax = max(0, ceil(max(dx)));
%   dxMin = min(0, floor(min(dx)));
%   dyMax = max(0, ceil(max(dy)));
%   dyMin = min(0, floor(min(dy)));
%   x = (1 + dxMax):(h + dxMin);
%   y = (1 + dyMax):(h + dyMin);
%   regMovie = regMovie(y,x,:);
% else
  x = 1:w;
  y = 1:h;
% end

% regMovie = [];

function [output Greg] = dftregistration(buf1ft,buf2ft,usfac)
% function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
% Efficient subpixel image registration by crosscorrelation. This code
% gives the same precision as the FFT upsampled cross correlation in a
% small fraction of the computation time and with reduced memory 
% requirements. It obtains an initial estimate of the crosscorrelation peak
% by an FFT and then refines the shift estimation by upsampling the DFT
% only in a small neighborhood of that estimate by means of a 
% matrix-multiply DFT. With this procedure all the image points are used to
% compute the upsampled crosscorrelation.
% Manuel Guizar - Dec 13, 2007

% Portions of this code were taken from code written by Ann M. Kowalczyk 
% and James R. Fienup. 
% J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued 
% object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458 
% (1990).

% Citation for this algorithm:
% Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup, 
% "Efficient subpixel image registration algorithms," Opt. Lett. 33, 
% 156-158 (2008).

% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)

% Outputs
% output =  [error,diffphase,net_row_shift,net_col_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0,
    CCmax = sum(sum(buf1ft.*conj(buf2ft))); 
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero); 
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    output=[error,diffphase];
        
% Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
% peak
elseif usfac == 1,
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc); 
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n); 
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax)); 
    md2 = fix(m/2); 
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end

    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
% Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    CC=zeros(mlarge,nlarge);
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));
  
    % Compute crosscorrelation and locate the peak 
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak 
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2 
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;

    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2,
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac; 
        col_shift = round(col_shift*usfac)/usfac;     
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid 
        [max1,loc1] = max(CC);   
        [max2,loc2] = max(max1); 
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);  
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;    

    % If upsampling = 2, no additional pixel shift refinement
    else    
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end  

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
return

function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1) 
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the 
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return
