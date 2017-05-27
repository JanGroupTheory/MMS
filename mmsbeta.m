function beta = mmsbeta(B,Te,ne) %mms units (nT, eV, cm^-3)
beta = B.^(-2).*ne.*Te*0.404385806370078;