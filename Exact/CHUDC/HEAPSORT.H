#include <math.h>
#define faux	0
#define vrai	1
#define ignore	2
#define	planch(x)	((int)(floor(x)))
#define INFINI	2147483647
typedef union{
	int 	*pint;
	double 	*pdoub;
	} PT;
typedef union{
	int	vint;
	double	vdoub;
	} VSW;
	
typedef	struct {
		int	ncrt;
		short	*signe;/*vrai->entier,faux->double*/
		short	*tend; /* *tend=0 ->croissant, 1 -> decroissant */
		PT	*crit;	//! Criterion, value to sort
		} Echelle;
int	EXP2[31]={1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,
		65536,131072,262144,524288,1048576,2097152,4194304,8388608,
		16777216,33554432,67108864,134217728,268435456,536870912,
		1073741824};

char	*malloc();

trie(code,n,nes,hie_crit)
short	code;
int	n,*nes;
Echelle hie_crit;
{

/*
	Code est vrai, On garde l'ancien tableaux, mais on met le nouvel 
		ordre dans NES;
	Sinon, On met les tableaux dans l'ordre qu'il faut 
*/

	int	i,l,j,depth,deb;

	if(n<=0) return;
	if(code) for(i=0;i<n;++i) *(nes+i)=i;
	depth=planch(log((double)n)/log(2.0));
	for(l=depth-1;l>=0;--l)
		for(deb=EXP2[l]-1,j=deb;j<=2*deb && j<n/2;++j)
			if(code)comspg(*(nes+j),j,n,nes,hie_crit);
			else	comspg1(j,j,n,hie_crit);
	if(code)enleve(n,nes,hie_crit);
	else enleve1(n,hie_crit);
}

	
comspg(memb,j,n,ele,hie)
int	memb,j,n,*ele;
Echelle hie;
{
	int m;

	while(j<n/2){
		m=selgr(j,n,ele,hie);
		if(test(*(ele+m),memb,hie))
		{
			*(ele+j)=(*(ele+m));
			j=m;
		}
		else	break;
	}
	*(ele+j)=memb;
}


comspg1(j1,j,n,hie)
int	j1,j,n;
Echelle hie;
{
	int	memb,m;
	VSW	*val;

	val=(VSW *)malloc((unsigned)(hie.ncrt)*sizeof(VSW));
	reservoir(val,j1,hie);
	while(j<n/2){
		m=selgr1(j,n,hie);
		if(test1(m,val,hie))
		{
			for(memb=0;memb<hie.ncrt;++memb)
				if(*(hie.signe+memb))
					*((hie.crit+memb)->pint+j)
					=(*((hie.crit+memb)->pint+m));
				else	*((hie.crit+memb)->pdoub+j)
					=(*((hie.crit+memb)->pdoub+m));
			j=m;
		}
		else break;
	}
	reprendre(val,j,hie);
	free((char *)val);
}

enleve(n,ele,hie)
int	n,*ele;
Echelle hie;
{
	int	i,remp;

	for(i=1;i<n;++i)
	{
		remp=(*ele);
		comspg(*(ele+n-i),0,n-i+1,ele,hie);
		*(ele+n-i)=remp;
	}
}

enleve1(n,hie)
int	n;
Echelle	hie;
{
	int	i;
	VSW	*remp;

	remp=(VSW *)malloc((unsigned)(hie.ncrt)*sizeof(VSW));
	for(i=1;i<n;++i)
	{
		reservoir(remp,0,hie);
		comspg1(n-i,0,n-i+1,hie);
		reprendre(remp,n-i,hie);
	}
	free((char *)remp);
}
		
reservoir(val,node,hie)
VSW	*val;
int	node;
Echelle hie;
{
	int	memb;

	for(memb=0;memb<hie.ncrt;++memb)
		if(*(hie.signe+memb))
			(val+memb)->vint=(*((hie.crit+memb)->pint+node));
		else	(val+memb)->vdoub=(*((hie.crit+memb)->pdoub+node));
}

reprendre(val,node,hie)
VSW	*val;
int	node;
Echelle hie;
{
	int	memb;

	for(memb=0;memb<hie.ncrt;++memb)
		if(*(hie.signe+memb))
			*((hie.crit+memb)->pint+node)=(val+memb)->vint;
		else	*((hie.crit+memb)->pdoub+node)=(val+memb)->vdoub;
}	

selgr1(x,n,hie)
int	x,n;
Echelle	hie;
{
	int	x1,x2,y;
	VSW	*val;

	x1=2*x+1;
	x2=2*x+2;
	val=(VSW *)malloc((unsigned)(hie.ncrt)*sizeof(VSW));
	reservoir(val,x2,hie);
	y=(x2>=n || test1(x1,val,hie))?x1:x2;
	free((char *)val);
	return(y);
}

selgr(x,n,ele,hie)
int	x,n,*ele;
Echelle hie;
{
	int	x1,x2;

	x1=2*x+1;
	x2=2*x+2;
	return((x2>=n || test(*(ele+x1),*(ele+x2),hie))?x1:x2);
}

test(j1,j2,hie)
int	j1,j2;
Echelle	hie;
{
	int	i,rep;
	VSW	val1,val2;
	PT	*pt;

	for(rep=ignore,i=0;i<hie.ncrt && rep==ignore;++i)
	{
		pt=hie.crit+i;
		if(*(hie.signe+i)==vrai)
		{
			val1.vint=(*(pt->pint+j1));
			val2.vint=(*(pt->pint+j2));
			if(val1.vint==val2.vint) rep=ignore;
			else	if(!*(hie.tend+i) && val1.vint<val2.vint ||
			 	(*(hie.tend+i)) && val1.vint>val2.vint)
					rep=faux;
				else	rep=vrai;
		}
		else{	
			val1.vdoub=(*(pt->pdoub+j1));
			val2.vdoub=(*(pt->pdoub+j2));
			if(val1.vdoub==val2.vdoub) rep=ignore;
			else	if(!*(hie.tend+i) && val1.vdoub<val2.vdoub ||
			 	(*(hie.tend+i)) && val1.vdoub>val2.vdoub)
					rep=faux;
				else	rep=vrai;
		}
	}
	return(rep);
}
 
test1(j1,j2,hie)
int	j1;
VSW	*j2;
Echelle	hie;
{
	int	i,rep;
	VSW	val1,val2;
	PT	*pt;

	for(rep=ignore,i=0;i<hie.ncrt && rep==ignore;++i)
	{
		val2=(*(j2+i));
		pt=hie.crit+i;
		if(*(hie.signe+i)==vrai)
		{
			val1.vint=(*(pt->pint+j1));
			if(val1.vint==val2.vint) rep=ignore;
			else	if(!*(hie.tend+i) && val1.vint<val2.vint ||
			 	(*(hie.tend+i)) && val1.vint>val2.vint)
					rep=faux;
				else	rep=vrai;
		}
		else{	
			val1.vdoub=(*(pt->pdoub+j1));
			if(val1.vdoub==val2.vdoub) rep=ignore;
			else	if(!*(hie.tend+i) && val1.vdoub<val2.vdoub ||
			 	(*(hie.tend+i)) && val1.vdoub>val2.vdoub)
					rep=faux;
				else	rep=vrai;
		}
	}
	return(rep);
}

ecritab(n,tab)
int	n,*tab;
{
	int	i;

	for(i=1;i<=n;++i)
		printf("%d%c",tab[i-1],(i%10==0 || i==n)?'\n':'\t');
}

