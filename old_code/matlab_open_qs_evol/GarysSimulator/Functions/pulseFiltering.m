function outPulse = pulseFiltering(inPulse,Ts,filterType,Ffilt,varargin)
%inPulse = pulse (filter on first dim)
%Ts = sampling time
%filterType = 'low pass', 'high pass', 'band pass'
%Ffilt = filter frequency (middle for band)
%varargin = band width in band pass

nx = size(inPulse,1);
nf = 2*nx-1; %zero padding
filter = zeros(nf,1);
minFreq = 1/(nx*Ts);
maxFreq = 1/(2*Ts);

function num = freq2num(freq)
    num = round(nx*(1+freq/maxFreq));
end

%Make perfect filter
switch(filterType)
    case 'low pass'
        if(Ffilt >= maxFreq)
            filter = ones(nf,1);
        elseif(Ffilt < maxFreq && Ffilt > minFreq)
            fm = freq2num(-Ffilt)+1;
            fp = freq2num(Ffilt)-1;
            filter(fm:fp) = 1;
        end
    case 'high pass'
        if(Ffilt <= minFreq)
            filter = ones(nf,1);
        elseif(Ffilt > minFreq && Ffilt < maxFreq)
            fm = freq2num(-Ffilt)+1;
            fp = freq2num(Ffilt)-1;
            filter([1:fm+1 fp:nf]) = 1;
        end
    case 'band pass'
        Df = round(varargin{1}/2);
        if(Df~=0) %need ALL conditions
            fmm = freq2num(-Ffilt-Df)+1;
            fmp = freq2num(-Ffilt+Df)+1;
            fpm = freq2num(Ffilt-Df)-1;
            fpp = freq2num(Ffilt+Df)-1;
            filter([fmm:fmp fpm:fpp]) = 1;
        end
end
if(~any(filter))
    error('Pulse filtering on pulse => 0');
end

%Make pulse fft
outPulse = zeros(size(inPulse));
for dim=1:size(inPulse,2)
    yfft = fftshift(fft(inPulse(:,dim),nf)); %fft + shift with 0 @ nx
    yfftF = filter.*yfft;
    yF = real(ifft(ifftshift(yfftF),nf));
    outPulse(:,dim) = yF(1:nx);
end

end