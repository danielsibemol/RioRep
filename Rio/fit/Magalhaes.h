#include<TComplex.h>

double rhoBC (double &sab){
	return sqrt(1.0 - 2.0*(mpi*mpi+mK*mK)/sab + (1.0/(sab*sab))*pow((mK*mK -mpi*mpi),2));
}

double ImL (double &sab){ // Parte Imagin·ria do loop
	double limiar;
	limiar = pow((mK +mpi),2);
	if (sab>limiar)
		return M_PI*rhoBC(sab);
	else
		return 0;
}

double I(double &sab){
	return (-1.0/(16.0*M_PI*M_PI))*ImL(sab);
}

double gamma2 (double &sab, double &M0, double &cd, double &cm){ // Kernel do espalhamento el·stico
	return ((M0*M0 - sab)/(Fkpi*Fkpi))*( sab - sab*(3.0/8.0)*pow(rhoBC(sab),2.0)  - (mpi*mpi + mK*mK)) 
		+ (3.0/(2.0*pow(Fkpi,4)))*pow((cd*sab - (cd - cm)*(mpi*mpi + mK*mK)),2);
}

double PartialWidth (double &sab, double &M0, double &cd, double &cm){ // Largura
	return gamma2(sab,M0,cd,cm)*I(sab);
}

double ReL (double &sab){ // Parte Real do loop
	double prelimiar,limiar,Sigma,L0, mK2,mpi2, M;
	limiar = pow((mK +mpi),2);
	prelimiar = pow( (mK - mpi),2 );
	mK2=mK*mK;
	mpi2= mpi*mpi;
	Sigma= mK2 +mpi2;
	M=1.0;//escala para manter adimencional  a princÌpio 1 GeV.
	L0=1.0 +((mK2+mpi2)/(mpi2-mK2))*log(mpi/mK) -  ((mpi2-mK2)/sab)*log(mpi/mK); // NOvo Valor subtraido em L(0).
	
	if (sab< prelimiar)
		return  L0 +(sqrt(lambda(sab,mK2,mpi2))/(2.0*sab))*log( (mK2 + mpi2 - sab  + sqrt(lambda(sab,mK2,mpi2)))/(mK2 + mpi2 - sab  - sqrt(lambda(sab,mK2,mpi2))));
	if (sab>=prelimiar && sab<Sigma)
		return  L0 - (sqrt(-lambda(sab,mK2,mpi2))/sab)*atan( sqrt(-lambda(sab,mK2,mpi2))/(mK2 +mpi2 -sab)) ;
	if (sab>=Sigma && sab<limiar)
		return L0 -(sqrt(-lambda(sab,mK2,mpi2))/sab)*(atan( sqrt(-lambda(sab,mK2,mpi2))/(mK2 +mpi2 -sab))
												 + M_PI) ;
	if (sab>=limiar)
		return L0 - (sqrt(lambda(sab,mK2,mpi2))/(2.0*sab))*log( (sab - mK2 - mpi2 + sqrt(lambda(sab,mK2,mpi2)))/ (sab - mK2 - mpi2 - sqrt(lambda(sab,mK2,mpi2)) ));
	else{
		cout << "ERROR: ReL in Magalhaes.h - conditionals do not apply" <<  endl;
		exit(-1);
		return 0;
	}
}

double Rbarra(double &sab){ // = 0 MATRIZ K
	return  (-1.0/(16.0*M_PI*M_PI))*ReL(sab);
}

double M2 (double &sab, double &M0, double &cd, double &cm){ //Running mass
	return M0*M0 + gamma2(sab,M0,cd,cm)*(Rbarra(sab)+ C);
}

TComplex Magalhaes(double &sab, double &M0, double &cd, double &cm){
	TComplex c_M2(M2(sab,M0,cd,cm),0);
	TComplex c_PartialWidth(PartialWidth(sab,M0,cd,cm),0);
	TComplex c_sab(sab,0);
	TComplex one;
	TComplex Amp = one.One()/(c_sab - c_M2 - one.I()*c_PartialWidth);
	return Amp;
}





