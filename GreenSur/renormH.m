function [Hs,Hb]=renormH(E,H,T,Ni)
	%iteration process
    TT  = T/(E-H);
    TTP = T'/(E-H);
	Hs  = H+TT*T';
	Hb  = H+TT*T'+TTP*T;
	Tt  = TT*T;
	Ttp = TTP*T';
	for i=1:Ni-1
        TTt  = Tt/(E-Hb);
        TTtp = Ttp/(E-Hb);
		Hst  = Hs+TTt*Ttp;
		Hbt  = Hb+TTt*Ttp+TTtp*Tt;
		Ttt  = TTt*Tt;
		Ttpt = TTtp*Ttp;
		Hs   = Hst;
		Hb   = Hbt;
		Tt   = Ttt;
		Ttp  = Ttpt;
	end
end