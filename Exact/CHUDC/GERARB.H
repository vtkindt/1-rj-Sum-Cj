ajarb(x,p2,lq,bint,fatt)
int x,*p2,*lq,*bint,*fatt;
{
	int i,y,j;

	for(i=fatt[*lq],bint[i]=x;i>0;)
	{
		j=(i-1)/2;
		if(p2[bint[i]]<p2[bint[j]] || p2[bint[i]]==p2[bint[j]] && bint[i]<bint[j])
		{
			y=bint[i];
			bint[i]=bint[j];
			bint[j]=y;
			i=j;
		}
		else	break;
	}
	if(*lq==0) ++fatt[0];
	else --(*lq);
}

modarb(n,p2,lq,bint,fatt)
int n,*p2,*lq,*bint,*fatt;
{
	int i,k1,k2,l;

	for(i=0;i<n && bint[i]!=-1;)
	{
		l=((k1=2*i+1)>=n || bint[k1]==-1)?2:0;
		l+=((k2=k1+1)>=n || bint[k2]==-1)?1:0;
		switch(l){
		case 0: if(p2[bint[k1]]<p2[bint[k2]] || p2[bint[k1]]==p2[bint[k2]] && bint[k1]<bint[k2]){
				bint[i]=bint[k1];
				i=k1;
			}
			else{
				bint[i]=bint[k2];
				i=k2;
			}
			break;
		case 1:	bint[i]=bint[k1];
			i=k1;
			break;
		case 2:	bint[i]=bint[k2];
			i=k2;
			break;
		case 3:	bint[i]=-1;
			++(*lq);
			fatt[*lq]=i;
			i=k1;
			break;
		}
	}
}
