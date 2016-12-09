
bool useVeto = 1;

void veto()
{

cout<<" vetoing sij..."<<endl;

}

bool veto(double sij)
{
	
if( useVeto == 0 ) return kTRUE;

if( sij>0.255 ) return kTRUE;
if( sij<0.235 ) return kTRUE;
 
 return kFALSE;

}