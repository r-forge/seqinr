computePI = function(seq){
	
	compoAA=table(factor(seq,levels=toupper(letters)))
	nTermR = which(toupper(letters)==seq[1])
	cTermR = which(toupper(letters)==seq[length(seq)])

	computeCharge <- function(pH ,compoAA ,pK = SEQINR.UTIL$pK , nTermResidue, cTermResidue){
		cter = 10^(-pK[cTermResidue,1]) / (10^(-pK[cTermResidue,1]) + 10^(-pH))	
		nter = 10^(-pH) / (10^(-pK[nTermResidue,2]) + 10^(-pH))
		carg = as.vector(compoAA['R'] * 10^(-pH) / (10^(-pK['R',3]) + 10^(-pH)))
		chis = as.vector(compoAA['H'] * 10^(-pH) / (10^(-pK['H',3]) + 10^(-pH)))
		clys = as.vector(compoAA['K'] * 10^(-pH) / (10^(-pK['K',3]) + 10^(-pH)))
		casp = as.vector(compoAA['D'] * 10^(-pK['D',3]) /(10^(-pK['D',3]) + 10^(-pH)))
		cglu = as.vector(compoAA['E'] * 10^(-pK['E',3]) / (10^(-pK['E',3]) + 10^(-pH)))
		ccys = as.vector(compoAA['C'] * 10^(-pK['C',3]) / (10^(-pK['C',3]) + 10^(-pH)))
		ctyr = as.vector(compoAA['Y'] * 10^(-pK['Y',3]) / (10^(-pK['Y',3]) + 10^(-pH)))
		charge = carg + clys + chis + nter - (casp + cglu + ctyr + ccys + cter)
	}

	critere <- function(pH,compoAA=compoAA,Pk=Pk,nTermResidue=nTermR,cTermResidue=cTermR){ 
		computeCharge(pH, compoAA, Pk, nTermResidue, cTermResidue)^2
		 }

	return(nlm(critere, 7,compoAA=compoAA,Pk=Pk,nTermResidue=nTermR,cTermResidue=cTermR)$estimate)
}


