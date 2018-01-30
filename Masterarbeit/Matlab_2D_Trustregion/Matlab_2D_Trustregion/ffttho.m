%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the Fast-Fourier-Transformation or the inverse FFT Y of the 
% input vector y.
%                                                                         
% Optional input argument/ value pairs:                                   
% 'nfft','[value]': Number of FFT acquisition points. Default: Next power 
%                   of two to the length of the input vector.             
% 'window','[value]': Window function to be applied to the input vector.  
%                     Possible values: 
%                     'barthannwin','bartlett',          
%                     'blackman','blackmanharris','bohmanwin','chebwin',  
%                     'flattopwin','gausswin','hamming','hann','kaiser',  
%                     'nuttallwin','parzenwin','rectwin','taylorwin',     
%                     'triang','tukeywin'                                 
% 'winopt','[value]': Optional input argument of window function. See the 
%                     documentation of the desired window for details
% 
%                                                                         
% written by: T. Hoffmann               
%
% See also: fft, ifft, fftshift, window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = ffttho(y,ikey,varargin) 

%Parse input arguments
p = inputParser;
defaultNFFT = 0;
defaultWindow = 'rectwin';
expectedWindow = {'barthannwin','bartlett','blackman','blackmanharris','bohmanwin','chebwin','flattopwin',...
    'gausswin','hamming','hann','kaiser','nuttallwin','parzenwin','rectwin','taylorwin','triang','tukeywin'};
addParameter(p,'nfft',defaultNFFT,@isnumeric);
addParameter(p,'window',defaultWindow, @(x) any(validatestring(x,expectedWindow)));
addParameter(p,'winopt',[]);
parse(p,varargin{:});
NFFT = p.Results.nfft;
winhandle = str2func(p.Results.window);
winopt = p.Results.winopt;

%Define FFT length
if NFFT == 0
    NFFT = 2^nextpow2(length(y));
end

%Define window
if isempty(winopt)
    win = window(winhandle,length(y));
else
    win = window(winhandle,length(y),winopt);
end

if ikey == 1
    %Use window on inputdata
    try
        y = y.*win;
    catch
       y = y.*win';
    end
    
    %Calculate FFT and shift zero frequency to the center of the output vector
    Y = fftshift(fft(y,NFFT))/length(y);
elseif ikey == -1
    Y = ifft(ifftshift(y),NFFT)*length(y);
else
    error('ikey is invalid');
end

