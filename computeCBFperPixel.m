function [modemap,ampmode,Peaks,Locations,specAmps,maskSpectrum]=computeCBFperPixel(ax,filt,L,Fs,cutStoN,promCut);

clear modemap ampmode Peaks Locations specAmps maskSpectrum
maskSpectrum=ones(size(filt,1),size(filt,2));
X=squeeze(filt(1,1,:));
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
%
spectouse=P1(3:end-1);
spectouse=spectouse-mean(spectouse(end-10:end-2));
specAmps=nan(size(filt,1),size(filt,2),length(spectouse));
cla(ax)
ylim(ax,[0,1])
xlim(ax,[0,1])
ph = patch(ax,[0 0 0 0],[0 0 1 1],[0.67578 1 0.18359]); %greenyellow
th = text(ax,1,1,'0%','VerticalAlignment','bottom','HorizontalAlignment','right');
for i=1:size(filt,1)
    for jj=1:size(filt,2)
        warning('off')
        X=squeeze(filt(i,jj,:));
        Y = fft(X);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        
        %
        spectouse=P1(3:end-1);
        spectouse=spectouse-mean(spectouse(end-10:end-2));
        if sum(diff(spectouse))~=0
            freqtouse=f(3:end-1);
            minHeight=median(spectouse)+2*std(spectouse);
            [peak,loc]=findpeaks(spectouse,freqtouse,'MinPeakHeight',minHeight,'MinPeakDistance',3);
            if ~isempty(peak)
                indfilt=find(spectouse<min(peak));
                minHeight=median(spectouse(indfilt))+2*std(spectouse(indfilt));
                [peak,loc,w,p]=findpeaks(spectouse,freqtouse,'MinPeakHeight',minHeight,'MinPeakDistance',3);
                ind1=find((peak/std(spectouse(indfilt)))>cutStoN);
                ind2=find(p>promCut);
                indp=intersect(ind1,ind2);
                % minHeight=median(P1(2:end))+2.5*std(P1(2:end));
                % [peak,loc]=findpeaks(P1(2:end),f(2:end),'MinPeakHeight',minHeight,'MinPeakDistance',1.5);
                peak=peak(indp);
                loc=loc(indp);
                Peaks{i,jj}=peak;
                Locations{i,jj}=loc;
                specAmps(i,jj,:)=spectouse;
                
                % keep only highest 3 modes
                [peaks,ind]=sort(peak);
                if length(loc)==1
                    modemap(i,jj,1)=loc;
                    ampmode(i,jj,1)=peak;
                elseif length(loc)==2                 
                for k=1:2
                    curfreq=loc(ind(length(loc)+1-k));
                    modemap(i,jj,k)=curfreq;
                    curamp=peak(ind(length(loc)+1-k));
                    ampmode(i,jj,k)=curamp;
                end
                elseif length(loc)==3                 
                for k=1:3
                    curfreq=loc(ind(length(loc)+1-k));
                    modemap(i,jj,k)=curfreq;
                    curamp=peak(ind(length(loc)+1-k));
                    ampmode(i,jj,k)=curamp;
                end
                end
                
                if isempty(Locations{i,jj})
                    maskSpectrum(i,jj)=0;
                end
                
            end
            
        end %check if non-zero spectrum
        
    end %parfor loop
    ph.XData = [0 i/size(filt,1)  i/size(filt,1) 0];
    th.String = sprintf('%.0f%%',round(i/(size(filt,1))*100));
    drawnow %update graphics
end% for loop