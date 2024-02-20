function [Dict,High,Medium,Low] = Create_Dictionary(sz,fs)
    % sz = size atom
   
    % Dictionary Elements
    Dict.A = [];
    Dict.B = [];
    Dict.C = [];
    Dict.D = [];
    Dict.E = [];
    Dict.F = [];
    Dict.G = [];
    Dict.H = [];
    Dict.I = [];
    Dict.J = [];
    Dict.K = [];
    Dict.L = [];

    % Dictionary Elements
    Dict.A = []; Dict.Aloc = []; Dict.Afrq = []; Dict.Awnd = []; 
    Dict.B = []; Dict.Bloc = []; Dict.Bfrq = []; Dict.Bwnd = []; 
    Dict.C = []; Dict.Cloc = []; Dict.Cfrq = []; Dict.Cwnd = []; 
    Dict.D = []; Dict.Dloc = []; Dict.Dfrq = []; Dict.Dwnd = []; 
    Dict.E = []; Dict.Eloc = []; Dict.Efrq = []; Dict.Ewnd = []; 
    Dict.F = []; Dict.Floc = []; Dict.Ffrq = []; Dict.Fwnd = []; 
    Dict.G = []; Dict.Gloc = []; Dict.Gfrq = []; Dict.Gwnd = []; 
    Dict.H = []; Dict.Hloc = []; Dict.Hfrq = []; Dict.Hwnd = []; 
    Dict.I = []; Dict.Iloc = []; Dict.Ifrq = []; Dict.Iwnd = []; 
    Dict.J = []; Dict.Jloc = []; Dict.Jfrq = []; Dict.Jwnd = []; 
    Dict.K = []; Dict.Kloc = []; Dict.Kfrq = []; Dict.Kwnd = []; 
    Dict.L = []; Dict.Lloc = []; Dict.Lfrq = []; Dict.Lwnd = []; 


    for ph = [0, pi/2]

        [TempA,TemplocA,TempfrqA,TempwdA] = Create_Gabor_Dictionary(sz,[1/4,1/8],[60 160],256,ph,fs,0.2);
        Dict.A=[TempA,Dict.A];
        Dict.Aloc = [TemplocA,Dict.Aloc];
        Dict.Awnd = [TempwdA,Dict.Awnd];
        Dict.Afrq = [TempfrqA,Dict.Afrq];

        [TempB,TemplocB,TempfrqB,TempwdB] = Create_Gabor_Dictionary(sz,[1/4,1/8],[60 160],256,ph,fs,0.8);
        Dict.B=[TempB,Dict.B];
        Dict.Bloc = [TemplocB,Dict.Bloc];
        Dict.Bwnd = [TempwdB,Dict.Bwnd];
        Dict.Bfrq = [TempfrqB,Dict.Bfrq];



        [TempC,TemplocC,TempfrqC,TempwdC] = Create_Gabor_Dictionary(sz,[1/8,1/16],[140 240],256,ph,fs,0.2);
        Dict.C=[TempC,Dict.C];
        Dict.Cloc = [TemplocC,Dict.Cloc];
        Dict.Cwnd = [TempwdC,Dict.Cwnd];
        Dict.Cfrq = [TempfrqC,Dict.Cfrq];

        [TempD,TemplocD,TempfrqD,TempwdD] = Create_Gabor_Dictionary(sz,[1/8,1/16],[140 240],256,ph,fs,0.8);
        Dict.D=[TempD,Dict.D];
        Dict.Dloc = [TemplocD,Dict.Dloc];
        Dict.Dwnd = [TempwdD,Dict.Dwnd];
        Dict.Dfrq = [TempfrqD,Dict.Dfrq];



        [TempE,TemplocE,TempfrqE,TempwdE] = Create_Gabor_Dictionary(sz,[1/8,1/16],[220 320],256,ph,fs,0.2);
        Dict.E=[TempE,Dict.E];
        Dict.Eloc = [TemplocE,Dict.Eloc];
        Dict.Ewnd = [TempwdE,Dict.Ewnd];
        Dict.Efrq = [TempfrqE,Dict.Efrq];

        [TempF,TemplocF,TempfrqF,TempwdF] = Create_Gabor_Dictionary(sz,[1/8,1/16],[220 320],256,ph,fs,0.8);
        Dict.F=[TempF,Dict.F];
        Dict.Floc = [TemplocF,Dict.Floc];
        Dict.Fwnd = [TempwdF,Dict.Fwnd];
        Dict.Ffrq = [TempfrqF,Dict.Ffrq];



        [TempG,TemplocG,TempfrqG,TempwdG] = Create_Gabor_Dictionary(sz,1/16,[300 400],256,ph,fs,0.2);
        Dict.G=[TempG,Dict.G];
        Dict.Gloc = [TemplocG,Dict.Gloc];
        Dict.Gwnd = [TempwdG,Dict.Gwnd];
        Dict.Gfrq = [TempfrqG,Dict.Gfrq];

        [TempH,TemplocH,TempfrqH,TempwdH] = Create_Gabor_Dictionary(sz,1/16,[300 400],256,ph,fs,0.8);
        Dict.H=[TempH,Dict.H];
        Dict.Hloc = [TemplocH,Dict.Hloc];
        Dict.Hwnd = [TempwdH,Dict.Hwnd];
        Dict.Hfrq = [TempfrqH,Dict.Hfrq];



        [TempI,TemplocI,TempfrqI,TempwdI] = Create_Gabor_Dictionary(sz,[1/16,1/32],[380 480],256,ph,fs,0.2);
        Dict.I=[TempI,Dict.I];
        Dict.Iloc = [TemplocI,Dict.Iloc];
        Dict.Iwnd = [TempwdI,Dict.Iwnd];
        Dict.Ifrq = [TempfrqI,Dict.Ifrq];

        [TempJ,TemplocJ,TempfrqJ,TempwdJ] = Create_Gabor_Dictionary(sz,[1/16,1/32],[380 480],256,ph,fs,0.8);
        Dict.J=[TempJ,Dict.J];
        Dict.Jloc = [TemplocJ,Dict.Jloc];
        Dict.Jwnd = [TempwdJ,Dict.Jwnd];
        Dict.Jfrq = [TempfrqJ,Dict.Jfrq];



        [TempK,TemplocK,TempfrqK,TempwdK] = Create_Gabor_Dictionary(sz,1/32,[460 520],256,ph,fs,0.2);
        Dict.K=[TempK,Dict.K];
        Dict.Kloc = [TemplocK,Dict.Kloc];
        Dict.Kwnd = [TempwdK,Dict.Kwnd];
        Dict.Kfrq = [TempfrqK,Dict.Kfrq];

        [TempL,TemplocL,TempfrqL,TempwdL] = Create_Gabor_Dictionary(sz,1/32,[460 520],256,ph,fs,0.8);
        Dict.L=[TempL,Dict.L];
        Dict.Lloc = [TemplocL,Dict.Lloc];
        Dict.Lwnd = [TempwdL,Dict.Lwnd];
        Dict.Lfrq = [TempfrqL,Dict.Lfrq];


    end

    Low.Atom = [Dict.A,Dict.B,Dict.C,Dict.D]; 
    Lowsize = size(Low,1);
    Low.Atom = normalize(Low.Atom,2);
    Low.frq = [Dict.Afrq,Dict.Bfrq,Dict.Cfrq,Dict.Dfrq];
    
    Medium.Atom  = [Dict.E,Dict.F, Dict.G, Dict.H]; 
    Mediumsize = size(Medium,2);
    Medium.Atom = normalize(Medium.Atom,2);
    Medium.frq = [Dict.Efrq,Dict.Ffrq, Dict.Gfrq, Dict.Hfrq]; 
    
    High.Atom = [Dict.I,Dict.J, Dict.K, Dict.L];
    Highsize = size(High,2);
    High.Atom = normalize(High.Atom,2);
    High.frq = [Dict.Ifrq,Dict.Jfrq, Dict.Kfrq, Dict.Lfrq];

end

