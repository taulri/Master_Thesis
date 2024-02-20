function [Dictionary,loc,frq,wd] = Create_Gabor_Dictionary(N,wnds,Freq,Overlap,phase,fs,alpha)
%[Dictionary] = Create_Gabor_Dictionary(Sample_Size,wnd,Freq,Overlap,phase,fs)
% input: N: size of the Atom
%        wnd: Windowing: 
%        Freq: Frequency range
%        Overlap: number of shifts
%        phase: phae of the DCT signal
%        fs: Sampling Frequency

   if nargin==6
        alpha = 0.95; % cosine fraction of the window, lower amplitude towards the ends --> the nearer to 0, the more rectangular 
   end

   if rem(log2(N),1)~=0
        error('Signal Length is not Dyadic. Please Choose a Dyadic Signal Length');
   end

    
   n = 1;
   for k = 1:length(wnds)
        wnd = wnds(k);
        s = length(wnd);
        for i=1:s
           Nw = N*wnd;
           Nw = floor(Nw);
           w=tukeywin(Nw,alpha);
           D = dctmtxphase(Nw,phase); % short atom snippets, length of the window 
           D=D'.*w; % width in time --> multiplying with window 
           Nf1=floor(Freq(1)/(fs/2)*Nw)+1;
           Nf2=floor(Freq(2)/(fs/2)*Nw)+1;  
               for f=Nf1:Nf2 % all the frequencies 
                  for k=1:Overlap+1-1 % shift in time 
                     h=zeros(1,N);
                     h((k)*N/Overlap)=1; %shift in time here, by 2 points 
                     wm=conv(h,D(:,f));
                     Dictionary(:,n) = wm([ceil(length(D)/2):(ceil(length(D)/2)+N-1)]);
                     if wnd<=1/4
                         nz = find(wm>0.0001);
                         loc(n) = (nz(end)-nz(1))/2+nz(1)-ceil(length(D)/2);
                         frq(n) = fs*(f-1)/(2*Nw);
                         wd(n) = wnd;
                     else
                         loc(n) = 0;
                         frq(n) = fs*f/(2*Nw);
                         wd(n) = wnd;
                     end
                     n=n+1;
                  end
                end
        end
    end


