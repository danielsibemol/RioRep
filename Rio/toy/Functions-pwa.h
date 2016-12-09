/****************************************************************************************/
/*               Function to define the limits of the mKK bins                          */
/****************************************************************************************/

void KK_edges(double m1, double m2, double m3, double B_mass, int npt, double *mKKlimits){

	double mKKmin, mKKmax, dmKK;
	int i;

	mKKmin = m1+m2+0.000001;
	mKKmax = B_mass - m3;
	dmKK = (mKKmax-mKKmin)/(npt-1);

	for(i=0;i<npt;i++)  mKKlimits[i] = mKKmin + i*dmKK;

	return;
}
/****************************************************************************************/
/*           Function to find out in which bin lies the given value of mKK              */
/****************************************************************************************/

void find_bin_KKP(double mKK, vector<double> mKKlimits, int npt, int &khi, int &klo){

	khi=-1;
	klo=-1;
	int i=0;

	for(i=1;i<npt;i++){
		if(mKK > mKKlimits[i-1] && mKK < mKKlimits[i]){
			klo=i-1;
			khi=i;
		}
	}

	if(klo==-1||khi==-1) cout << "ops, bin = -1, mKK = " << mKK << endl;

	return;

}

void find_bin( double m12, double m13, vector<double> mKKlimits, int npt, int &khi12, int &klo12, int &khi13, int &klo13 ){


	if( m12 > 1.73 || m12 < 0.279 || m13 > 1.73 || m13 < 0.279 ) cout<<" Houston, temos um problema. "<<endl;

	khi12=-1;
	klo12=-1;
	khi13=-1;
	klo13=-1;

	for( int  i = 0;i < npt - 1; i++ ){
		if( m12 >= mKKlimits[ i ] && m12 < mKKlimits[ i + 1 ]){
			klo12 = i;
			khi12 = i + 1;
		}
		if(m13 >= mKKlimits[ i ] && m13 < mKKlimits[ i + 1 ]){
			klo13 = i;
			khi13 = i + 1;
		} 										


	}

	if(klo12==-1||khi12==-1) cout << "Heyyyy!!!, bin = -1, m12 = " << m12 << endl;
	if(klo13==-1||khi13==-1) cout << "Heyyyy!!!, bin = -1, m13 = " << m13 << endl;   

	return;

}



/****************************************************************************************/
/*        Function to compute the second derivatives needed for interpolation           */
/*          Taken from Numerical Recipes in C - second edition - page 115               */
/*                 (adapted to consider zero-offset arrays)                             */
/****************************************************************************************/

void derivative(vector<double> x, vector<double> y, int n, double yp1, double ypn, vector<double> &y2){
	int i,k;
	const int N=n;
	double p, qn, sig, un, u[N];

	/* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/

	y2.resize(n);

	if(yp1 > 0.99e30) {
		y2[0] = 0.;
		u[0] = 0.;
	}
	else{
		y2[0]=-0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}

	/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors*/

	for(i=1;i<n-1;i++){
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(ypn > 0.99e30) {
		qn = 0.;
		un = 0.;
	}
	else{
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

	/* This is the backsubstitution loop of the tridiagonal algorithm */

	for(k=n-2;k>=0;k--) y2[k]=y2[k]*y2[k+1]+u[k];

}

/*
   double fprime(vector<double> x, vector<double> y, Int pm) {

   Doub s1 = x[0]-x[pm*1], s2 = x[0]-x[pm*2], s3 = x[0]-x[pm*3], s12 = s1-s2, s13 = s1-s3, s23 = s2-s3;
   return -(s1*s2/(s13*s23*s3))*y[pm*3]+(s1*s3/(s12*s2*s23))*y[pm*2]
   -(s2*s3/(s1*s12*s13))*y[pm*1]+(1./s1+1./s2+1./s3)*y[0];
   }
 */

void fprime(vector<double> x, vector<TComplex> y, TComplex &yp1 ) {

	double s1 = x[0]-x[1],
	       s2 = x[0]-x[2], 
	       s3 = x[0]-x[3], 
	       s12 = s1-s2, 
	       s13 = s1-s3, 
	       s23 = s2-s3;

	yp1 = -( s1*s2/( s13*s23*s3 ) )*y[3]+( s1*s3/( s12*s2*s23 ) )*y[2]
		-( s2*s3/( s1*s12*s13 ) )*y[1]+( 1./s1+1./s2+1./s3 )*y[0];
}

void fprime(vector<double> x, vector<TComplex> y, int n, TComplex &yp2 ) {

	double s1 = x[n-1] - x[n-2], 
	       s2 = x[n-1] - x[n-3], 
	       s3 = x[n-1] - x[n-4], 
	       s12 = s1-s2, 
	       s13 = s1-s3, 
	       s23 = s2-s3;

	yp2 = -( s1*s2/( s13*s23*s3 ) )*y[n-4]+( s1*s3/( s12*s2*s23 ) )*y[n-3]
		-( s2*s3/( s1*s12*s13 ) )*y[n-2]+( 1./s1+1./s2+1./s3 )*y[n-1];
}

/****************************************************************************************/
/*        Function to compute the second derivatives needed for interpolation           */
/*          Taken from Numerical Recipes in C - second edition - page 115               */
/*                 (adapted to consider zero-offset arrays)                             */
/****************************************************************************************/


void NotComplex_Derivative(vector<double> x, vector<TComplex> y, int n, vector<TComplex> &y2){
	int i,k;
	const int N=n;
	double p, qn, sig, un, u[N];
	TComplex yp1, ypn;
	vector<double> y2re, y2im;

	yp1 = 2.*(y[1] - y[0]) / (x[1] - x[0]);
	ypn = 2.*(y[n-1] - y[n-2]) / (x[n-1] - x[n-2]);

	y2.resize(n);
	y2re.resize(n);
	y2im.resize(n);	    

	/* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(yp1.Re() > 0.99e30) {
		y2re[0] = 0.;
		u[0] = 0.;
	}
	else{
		y2re[0]=-0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1].Re()-y[0].Re())/(x[1]-x[0])-yp1.Re());
	}


	/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors*/

	for(i=1;i<n-1;i++){
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2re[i-1]+2.0;
		y2re[i]=(sig-1.0)/p;
		u[i]=(y[i+1].Re()-y[i].Re())/(x[i+1]-x[i]) - (y[i].Re()-y[i-1].Re())/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(ypn.Re() > 0.99e30) {
		qn = 0.;
		un = 0.;
	}
	else{
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn.Re()-(y[n-1].Re()-y[n-2].Re())/(x[n-1]-x[n-2]));
	}
	y2re[n-1]=(un-qn*u[n-2])/(qn*y2re[n-2]+1.0);

	/* This is the backsubstitution loop of the tridiagonal algorithm */

	for(k=n-2;k>=0;k--) y2re[k]=y2re[k]*y2re[k+1]+u[k];

	/* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(yp1.Im() > 0.99e30) {
		y2im[0] = 0.;
		u[0] = 0.;
	}
	else{
		y2im[0]=-0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1].Im()-y[0].Im())/(x[1]-x[0])-yp1.Im());
	}

	/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors*/

	for(i=1;i<n-1;i++){
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2im[i-1]+2.0;
		y2im[i]=(sig-1.0)/p;
		u[i]=(y[i+1].Im()-y[i].Im())/(x[i+1]-x[i]) - (y[i].Im()-y[i-1].Im())/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}

	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(ypn.Im() > 0.99e30) {
		qn = 0.;
		un = 0.;
	}
	else{
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn.Im()-(y[n-1].Im()-y[n-2].Im())/(x[n-1]-x[n-2]));
	}
	y2im[n-1]=(un-qn*u[n-2])/(qn*y2im[n-2]+1.0);

	/* This is the backsubstitution loop of the tridiagonal algorithm */

	for(k=n-2;k>=0;k--) y2im[k]=y2im[k]*y2im[k+1]+u[k];

	//	y2[0] = yp1;
	//	y2[n-1] = ypn;

	y2[0] = 0;
	y2[n-1] = 0;
	for(k=n-2;k>0;k--) y2[k](y2re[k],y2im[k]);
}


void Complex_Derivative(vector<double> x, vector<TComplex> y, int n, vector<TComplex> &y2){
	int i,k;
	const int N=n;
	double p, qn, sig, un, u[N];
	TComplex yp1, ypn;
	vector<double> y2re, y2im;

	fprime( x, y, yp1 );
	fprime( x, y, n, ypn );	

	//cout<<"yp1= "<<yp1<<" ypn = "<<ypn<<endl; 

	y2.resize(n);
	y2re.resize(n);
	y2im.resize(n);	    

	/* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(yp1.Re() > 0.99e30) {
		y2re[0] = 0.;
		u[0] = 0.;
	}
	else{

		y2re[0] = -0.5;
		u[0] = (3.0/( x[1] - x[0] ) )*( ( y[1].Re() - y[0].Re() )/( x[1] - x[0] ) - yp1.Re() );

		//u[0] = (3.0 / (x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
	}


	/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors*/

	for( i = 1;i < n - 1; i++ ){

		sig       = ( x[i] - x[ i - 1 ] )/( x[ i + 1 ] - x[ i - 1 ] );
		p         =   sig*y2re[ i - 1 ] + 2.0;
		y2re[ i ] = ( sig - 1.0 )/p;
		u[ i ]    = ( y[ i + 1 ].Re() - y[ i ].Re())/( x[ i + 1 ] - x[ i ] ) - ( y[ i ].Re() - y[ i - 1 ].Re())/( x[ i ] - x[ i - 1 ] );
		u[ i ]    = ( 6.0*u[ i ]/( x[ i + 1 ] - x[ i - 1 ] ) - sig*u[ i - 1 ] )/p;

		//	u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
		//    u[i] = (6.0*u[i] / (x[i+1] - x[i-1]) - sig*u[i-1]) / p;

	}

	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(ypn.Re() > 0.99e30) {
		qn = 0.;
		un = 0.;
	}
	else{
		qn = 0.5;
		un = ( 3.0/( x[ n - 1 ] - x[ n - 2 ] ) )*( ypn.Re() - ( y[ n - 1 ].Re() - y[ n - 2 ].Re() )/( x[ n - 1 ] - x[ n - 2 ] ) );
		//	un = (3.0 / (x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));

	}

	y2re[ n - 1 ] = ( un - qn*u[ n - 2 ] )/( qn*y2re[ n - 2 ] + 1.0 );

	/* This is the backsubstitution loop of the tridiagonal algorithm */

	for( k = n - 2; k >= 0; k-- ) y2re[k]=y2re[k]*y2re[k+1]+u[k];

	/* The lower boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(yp1.Im() > 0.99e30) {

		y2im[0] = 0.;
		u[0] = 0.;

	}
	else{

		y2im[0]=-0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1].Im()-y[0].Im())/(x[1]-x[0])-yp1.Im());

	}

	/* This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors*/

	for(i=1;i<n-1;i++){

		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2im[i-1]+2.0;
		y2im[i]=(sig-1.0)/p;
		u[i]=(y[i+1].Im()-y[i].Im())/(x[i+1]-x[i]) - (y[i].Im()-y[i-1].Im())/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;

	}

	/* The upper boundary condition is set either to be "natural" or else to have specified first derivative*/

	if(ypn.Im() > 0.99e30) {

		qn = 0.;
		un = 0.;

	}
	else{

		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn.Im()-(y[n-1].Im()-y[n-2].Im())/(x[n-1]-x[n-2]));
	}

	y2im[n-1]=(un-qn*u[n-2])/(qn*y2im[n-2]+1.0);

	/* This is the backsubstitution loop of the tridiagonal algorithm */


	for(k=n-2;k>=0;k--) y2im[k]=y2im[k]*y2im[k+1]+u[k];

	/****************************************************************************************/

	//	y2[0] = yp1;
	//y2[0](y2re[0],y2im[0]);

	y2[0](0,0);
	y2[n-1](0,0);
	for(k=n-2;k>=1;k--) y2[k](y2re[k],y2im[k]);


}

/****************************************************************************************/
/*                                                                                      */
/****************************************************************************************/

void binning(){


}
