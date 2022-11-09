#include <R.h>
#include <Rinternals.h>
#include <math.h>


void blasso(double * B, double * S, double *M, double *Om, double * soft, int * pin,
            int * qin, int * nin, double * lam, double * tol, int * maxit, int * totalit, double * objective)
{
    double obj=*objective;
    int n = *nin;
    int p=*pin;
    int q=*qin;
    int kit,r,c, a, b;
    double  bnew, H, AH, tmp, tmp2, thisdiff, this_obj_diff, max_diff;
    double *b_rc, *b_ab, *lam_rc, *soft_rc, * M_rc , *Om_cc, *Om_c0, *S_rr, *S_ar, *soft_ab, *Om_cb, *S_0r;

    kit=0;
    max_diff=*tol+1.0;
    while( (max_diff > *tol) && (kit < *maxit) )
    {
        kit++;
        max_diff = 0.0;
        soft_rc=soft;
        lam_rc=lam;
        M_rc=M;
        Om_cc=Om;
        b_rc=B;
        Om_c0=Om;
        for(c=0; c <q; c++)
        {
            S_rr=S;
            S_0r=S;
            for(r=0; r < p; r++)
            {
                H = *soft_rc;
                AH=fabs(H);
                tmp=AH - *lam_rc;
                bnew=0.0;
                if(tmp  > 0.0 )
                {
                    if(H > 0.0)
                        bnew = tmp;
                    else if( H < 0.0 )
                        bnew = -tmp;
                    else
                        bnew=0.0;
                }
                bnew=bnew/(*S_rr * * Om_cc);
                if(bnew != *b_rc)
                {
                    thisdiff=*b_rc-bnew;
                    this_obj_diff=0.0;
                    soft_ab=soft;

                    Om_cb=Om_c0;
                    b_ab=B;
                    for(b=0; b<q;b++) //column  loop
                    {
                        S_ar=S_0r;
                        for(a=0; a<p; a++) //row loop
                        {
                            if(!( (a==r) && (b==c)))
                            {
                                tmp2=*S_ar * *Om_cb * thisdiff;
                                *soft_ab+= tmp2;
                                this_obj_diff+= *b_ab *tmp2;
                            }
                            soft_ab++;
                            b_ab++;
                            S_ar++;
                        }
                        Om_cb++;
                    }
                    this_obj_diff+= thisdiff*(-*M_rc - *soft_rc + *S_rr * * Om_cc * (*b_rc + bnew) );
                    this_obj_diff+=  2.0* *lam_rc *(fabs(*b_rc) - fabs(bnew));
                    this_obj_diff=this_obj_diff/n;
                    *b_rc=bnew;
                    obj-=this_obj_diff;
                    if ( this_obj_diff > max_diff )
                        max_diff=this_obj_diff;
                }
                M_rc++;
                lam_rc++;
                soft_rc++;
                b_rc++;
                S_0r+=p;
                S_rr+=p+1;
            } // end r loop
            Om_c0+=q;
            Om_cc+=q+1;
        } // end c loop
    }
    totalit[0]=kit;
    objective[0]=obj;
}
