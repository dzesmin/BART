
#include <transit.h>

/* initial version March 2nd, 2014 Jasmina Blecic */


/* \fcnfh
   Computes optical depth for eclipse geometry.
   Returns: Optical depth divided by rad.fct:  \frac{tau}{units_{rad}} */       


/* JB instead of prinf to transitprinf verbose level */
static PREC_RES
totaltau_eclipse( PREC_RES *rad,   /* Equispaced (part) of radius array from current layer to outmost layer */     
                  PREC_RES *ex,    /* Extinction[rad], extinction array for all radii and one wn */  /* This often used as ex+ri rather than &ex[ri] */
                    long nrad)    /* Number of radii from current layer to outmost layer */   
{
  PREC_RES res;  /* Optical depth  result in units of radius */
  PREC_RES x3[3],r3[3]; /*  interpolation variables          */

  /* Distance between two radii. Radius array needs to be equispaced.  */
  const PREC_RES dr=rad[1]-rad[0];

  /* Distance along the path */
  PREC_RES s[nrad];

  /* Returns 0 if at outermost layer. No distance travelled. */
  if(nrad == 1)  
    return 0.;


  /* Providing 3 necessary points for spline integration */
  const PREC_RES tmpex=*ex;
  const PREC_RES tmprad=*rad;
  if(nrad==2) *ex=interp_parab(rad-1,ex-1, rad[0]);
  else *ex=interp_parab(rad,ex, rad[0]);

  if(nrad==2){
    x3[0]=ex[0];
    x3[2]=ex[1];
    x3[1]=(ex[1]+ex[0])/2.0;
    r3[0]=rad[0];
    r3[2]=rad[1];
    r3[1]=(rad[0]+rad[1])/2.0;
    *rad=tmprad;
    *ex=tmpex;
    rad=r3;
    ex=x3;
    nrad++;
  }


  /* Distance along the path*/
  s[0] = 0.;
//  printf("s[0]=%g, dr=%g, nrad= %li", s[0], dr, nrad);

  for(int i=1; i< nrad; i++){
    s[i]= s[i-1] + dr;          
//    printf("s[%i]=%g= s[%i] + %g= %g + %g, dr=%g\n", i, s[i], i-1, dr, s[i-1], dr, dr);  
  }


//  printf("Pre integracije: %li\n", nrad);

  /* Integrate extinction along the path: */                          
  /* Use spline if GSL is available along with at least 3 points: */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl = gsl_interp_alloc(gsl_interp_cspline, nrad);
  gsl_interp_init(spl, s, ex, nrad);
  res = gsl_interp_eval_integ(spl, s, ex, 0, s[nrad-1], acc);    
  gsl_interp_free(spl);
  gsl_interp_accel_free(acc);
#else
#error non equispaced integration is not implemented without GSL
#endif /* _USE_GSL */

//  printf("Posle integracije: %li\n", nrad);

  return res;
}





/*
 * tau.c   - Fills out the 2D array of optical depths for each radius and each wavelengh
 */

#include <transit.h>
#include <extraext.h>

#define CIA_DOFLOAT  2
#define CIA_RADFIRST 1

/* \fcnfh
   Calculates optical depth as a function of radii and wavenumber.

   Return: 0 on success   */


int
tau_eclipse(struct transit *tr){
  /* Different structures */
  static struct optdepth tau;            /* Defines optical depth structure */
  tr->ds.tau = &tau;                               
  prop_samp *rad = &tr->rads;              /* Radius sample*/   
  prop_samp *wn  = &tr->wns;               /* Wavenumber sample */ 
  struct extinction *ex = tr->ds.ex;      /* Defines extinction struct */ 
  PREC_RES **e = ex->e[tr->tauiso];        /* Defines 2D extinction array for one isotope all radii and all wavenumbers */         
  PREC_RES (*fcn)() = tr->sol->tauEclipse; /* Defines function fcn that calls totaltau_eclipse transit path solution to calculate optical depth for one ray for all radii at one wn, see structure_tr.h line 91*/ 



  long wi, ri; /* Indices for wavenumber, and radius */
  int rn;   /* JB variable that return error code from computeextradiaus=>exstradisu */


  PREC_RES *r  = rad->v;         /* Values of radii */
  PREC_RES *tau_wn;              /* Declares a pointer to an array of optical depths for 1 wn at all radii*/         
  PREC_ATM *temp = tr->atm.t,    /* Takes temperatures from the atmospheric file        */
           tfct  = tr->atm.tfct; /* Temperature units  from the atmospheric file  */ 

  /* Number of elements: see line 30 in stucture_tr.h*/
  long int wnn = wn->n;   /* Wavenumbers      */   
  long int rnn = rad->n;  /* Radii            */
  double wfct  = wn->fct; /* Wavenumber units factor to cgs */


/* Jasmina include this without refraction index */
//  transitcheckcalled(tr->pi,"tau",2,
//		     "idxrefrac",TRPI_IDXREFRAC,
//		     "extwn",TRPI_EXTWN
//		     );

/* Jasmina ask Pato roho, what this TRU_TAUBITS does, and that it is used just once in the code */
//  transitacceptflag(tr->fl,tr->ds.th->fl,TRU_TAUBITS);


  /* Sets transit maximum optical depth to calculate: */
  struct transithint *trh=tr->ds.th;
  tau.toomuch = 10;  /* Default value */
  if(tr->ds.th->toomuch>0)   /* JB put space before and after the >< sign */
    tau.toomuch = trh->toomuch; 


  /* Aloccates array to store radii indices where tau reached toomuch: */
  tau.last = (long      *)calloc(wnn,       sizeof(long));


  /* JB Allocates tau.t as 2D array (wn, rad) */
  tau.t    = (PREC_RES **)calloc(wnn,       sizeof(PREC_RES *));
  tau.t[0] = (PREC_RES  *)calloc(wnn*rad->n, sizeof(PREC_RES));


  /* Connects tau.t with correct wn and rad */
  for(wi=1; wi<wnn; wi++)
    tau.t[wi] = tau.t[0] + wi*rad->n;

  /* Set cloud structure:   JB this should be excluded         */
  static struct extcloud cl;
  cl.maxe = tr->ds.th->cl.maxe; /* Maximum opacity  */
  cl.rini = tr->ds.th->cl.rini; /* Top layer radius */
  cl.rfin = tr->ds.th->cl.rfin; /* Radius of maxe   */
  if(tr->ds.th->cl.rfct==0)
    cl.rfct = rad->fct;
  else
    cl.rfct = tr->ds.th->cl.rfct;
  tr->ds.cl = &cl;


  /* Tests if ex is calculated for the particular radius, if true/False: */
  _Bool *comp = ex->computed;


  /* JB prints ex values before outemost layer is computed: */
  for(wi=0; wi<wnn; wi++)
    if(wi<2){
      printf("Exctinction before outermost layes is computed for w[%li]=%9.4g; has value of excintion is:", wi, wn->v[wi]);
      for (int i=rnn-1; i>rnn-3; i--)
        printf("%9.4g;", e[i][wi]);
      printf("\n"); 
  }


  /* Restoring ex from the savefile if given: */
  if(tr->save.ext)
    restfile_extinct(tr->save.ext, e, comp, rnn, wnn);


  /* Computes extinction at the outermost layer: */
  /* JB put Pato explaniation in */
  if(!comp[rnn-1]){
    transitprint(1, verblevel,
                 "Computing extinction in the outtermost layer.\n");
    if((rn=computeextradius(rnn-1, tr->atm.t[rnn-1]*tr->atm.tfct, ex))!=0)
      transiterror(TERR_CRITICAL,
                   "computeextradius() returned error code %i.\n", rn);
  }


  /* JB prints ex values after outemost layer is computed: */
  for(wi=1; wi<wnn; wi++)
    if(wi<2){
      printf("Exctinction after outermost layes is computed for w[%li]=%9.4g; has value of excintion is:", wi, wn->v[wi]);
      for (int i=rnn-1; i>rnn-3; i--)
        printf("%9.4g;", e[i][wi]);
      printf("\n"); 
  }



  /* Gets a copy of the radius units factor */
  double rfct = rad->fct;


  /* Requests at least four radii samples to calculate
     a spline interpolation:                                     */
  if(rnn<4)
    transiterror(TERR_SERIOUS,
                 "tau(): At least four radius points are "
                 "required! (three for spline and one for the analitical "
                 "part)");

  transitprint(1, verblevel,
                  "Calculating optical depth at various radius ...\n");



  /* Note that it works only for one isotope */
  if(ex->periso)
    transitprint(2, verblevel,
                 "Computing only for isotope '%s', others were ignored.\n",
                 tr->ds.iso->isof[tr->tauiso].n);

  PREC_RES er[rnn];        /* Array of extinction per radius           */


  int lastr = rnn-1;       /* Radius index of last computed extinction, starts from the outmost layer */
  int wnextout = (long)(wnn/10.0); /* Tenth of wavenumber sample size,
                                      used for progress printing       */


  /* Extinction from scattering and clouds: */
  double e_s[rnn], e_c[rnn];
  /* Extinction from CIA:                   */
  PREC_CIA **e_cia = tr->ds.cia->e;
  /* Extinction from scattering:            */
  struct extscat *sc = tr->ds.sc;



  /* JB pasted from the thesis, explains the flow of this part
For each wavenumber of the selected scale, the optical depth along a tangential
path with a given impact parameter is computed using Eq. (3.27). The
computation starts with an impact parameter, ρ0, equal to the radius of the
uppermost layer and descends into the atmosphere until the tangential optical
depth rises above a certain user-defined limit (--toomuch). We find that
that value should be at least 5 for convergent results. Each time this procedure
reaches a radius whose extinction has not been computed, it calls the extinctioncomputing
routine, which uses Eq. (3.32) to compute the extinction spectrum
over the whole wavenumber range at the specified radius. */
  

  /* For each wavenumber: */
  for(wi=0; wi<wnn; wi++){

    /* Pointing to the declared array all radii for one wn */
    tau_wn = tau.t[wi];

    /* Print output every 10\% that is ready: */
    if(wi>wnextout){
      transitprint(2, verblevel, "%i%%\n", (int)(100*(float)wi/wnn+0.5));
      wnextout += (long)(wnn/10.0);
    }

    /* Calculate extinction from scattering, clouds, and CIA at each level just for the wn of interest: */
    computeextscat(e_s,  rnn, sc, rad->v, rad->fct, temp, tfct, wn->v[wi]*wfct);
    computeextcloud(e_c, rnn, &cl, rad, temp, tfct, wn->v[wi]*wfct);


    /* Calculate exctintion for all radii if excition is already exists, equation 3.32 */
    for(ri=0; ri<rnn; ri++)
      er[ri] = e[ri][wi] + e_s[ri] + e_c[ri] + e_cia[wi][ri];


//    for(ri=rnn-1; ri>-1; ri--)
//      printf ("Radius (r[%i]=%9.4g\n",ri, r[ri]);

    /* For each radii: */
    for(ri=rnn-1; ri>-1; ri--){

//      printf ("Radius r[%li]=%9.4g\n", ri, r[ri]);

      /* Computes extinction only if radius is smaller then the radius of the last calucalted exctinction. He also converts to the same units ip and r */
      /* JB in my case, it is simplier to say if (ri < lastr){  */
      if(r[ri]*rad->fct < r[lastr]*rfct){    


        printf ("Last Tau (r=%9.4g, wn=%9.4g): %10.4g.\n",
                                     r[ri-1], wn->v[wi], tau_wn[ri-1]);

        if(ri)
          transitprint(3, verblevel, "Last Tau (r=%9.4g, wn=%9.4g): %10.4g.\n",
                                     r[ri-1], wn->v[wi], tau_wn[ri-1]);

        /* If extinction was not computed, compute it using computextradius that calls extradius for one radius, for all wn, but then update extintion for the particular wn. Starts from the layer below outermost*/
        do{
          if(!comp[--lastr]){
            /* Compute a new extinction at given radius printing error if
               something happen: */
            transitprint(2, verblevel, "Radius %i: %.9g cm ... ",
                                       lastr+1, r[lastr]);
            if((rn=computeextradius(lastr, temp[lastr]*tfct, ex))!=0)
              transiterror(TERR_CRITICAL,
                           "computeextradius() return error code %i while "
                           "computing radius #%i: %g\n", rn, r[lastr]*rfct);
            /* Update the value of the extinction at the right place: */
            er[lastr] = e[lastr][wi] + e_s[lastr] +
                        e_c[lastr] + e_cia[wi][lastr];
          }
        }while(r[ri]*rad->fct < r[lastr]*rfct);
      }



      /* Check if tau reached toomuch, and fills out tau_wn[ri] up to tau=toomuch: */   /* JB check with Patricio about the units for bb fct/rfct */
      if( (tau_wn[rnn-ri-1] = rfct * fcn(r+ri, er+ri, rnn-ri)) > tau.toomuch){
        tau.last[wi] = rnn-ri-1;
/* JB ovde se tau puni tako da je tau[0] je outer most layer i ide do tau[?] dokle je tauIndex=toomuch */
//  double rfct=rad->fct;
//  double riw=ip->fct/rfct;
        if (ri<3){
          transitprint(1, verblevel,
                       "WARNING: At wavenumber %g (cm-1), the critical TAU "
                       "value (%g) was exceeded with tau=%g at the radius "
                       "level %li (%g km), this should have "
                       "happened in a deeper layer (check IP sampling or ATM "
                       "file).\n", wn->v[wi], tau.toomuch, tau_wn[ri],
                       ri, r[ri]*rfct/1e5);
        }
        break;
      }

      /* Sets that tau of the outermost layer is zero */
      tau_wn[0] = 0;   

      transitDEBUG(22, verblevel,
                   "Tau(lambda %li=%9.07g, r=%9.4g) : %g  (toomuch: %g)\n",
                   wi, wn->v[wi], r[ri], tau_wn[ri], tau.toomuch);
    }

    if(ri==rnn){
      transitprint(1, verblevel,
                   "WARNING: At wavenumber %g cm-1, the bottom of the "
                   "atmosphere was reached before obtaining the critical TAU "
                   "value of %g.\nMaximum TAU reached: %g.\n",
                   wn->v[wi], tau.toomuch, tau_wn[ri]);
      tau.last[wi] = ri-1;
    }

  }

  transitprint(1, verblevel, " Done.\nOptical depth calculated up to %g.\n",
               tr->ds.tau->toomuch);

  /* Print detailed output if requested: */
  if(tr->ds.det->tau.n)
    detailout(&tr->wns, &tr->rads,  &tr->ds.det->tau, tau.t, 0);
  if(tr->ds.det->ext.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->ext, e, CIA_RADFIRST);
  if(tr->ds.det->cia.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->cia, (double **)e_cia,
              CIA_DOFLOAT);

  if(tr->save.ext)
    savefile_extinct(tr->save.ext, e, comp, rnn, wnn);

//  printf(" Radius Index,   Radius,       Wn,      Tau\n");
  for(wi=0; wi<wnn; wi++){
    int rIndex= tau.last[wi]; 
    double radiusLastWn= r[rIndex];
    double wNum= wn->v[wi];

//    printf("%13i,%9.4g,%9.4g,%9g\n", rIndex, radiusLastWn, wNum, tau_wn[rIndex]);
  }

  /* Print lowest radius before optical depth gets too big: */
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch, tr->ds.tau, &tr->wns, &tr->rads);

  /* Free memory that is no longer needed: */
  freemem_lineinfotrans(tr->ds.li, &tr->pi);
  freemem_localextinction();

  /* Set progress indicator and output tau if requested, otherwise return
     success: */
  tr->pi |= TRPI_TAU;
  if(tr->fl & TRU_OUTTAU)
    printtau(tr);
  return 0;
}



static PREC_RES
eclipse_intens(struct transit *tr, /* Main structure */
            PREC_RES *tau,    /* Optical depth array taken from tau.c==>tau             */ 
            PREC_RES w,              /* Current wavenumber value    */      
            long last,              /* Index where tau = toomuch        */
            double toomuch,         /* Maximum optical depth calculated */
            prop_samp *rad){         /* Radii array */


  /* General variables: */
  PREC_RES res;
  PREC_ATM *temp = tr->atm.t;    /* Temperatures        */
//  PREC_ATM tfct  = tr->atm.tfct; /* Temperature units   */ /* JB see about these units */
//  printf("tfct%g\n", tfct);

  prop_samp *wn = &tr->wns; /* JB takes samling propertioes for wavenumber from tr */
  double wfct  = wn->fct; /* Wavenumber units factor to cgs */
//  printf("wfct%9.4g\n", wfct);

  /* Radius parameter variables: */
  long rnn  = rad->n;
  long i;

  /* Blackbody variables: */
  PREC_RES B[rnn];

  /* Integration parts*/
  PREC_RES tauInteg[rnn],  /* Integrand function */
           tauIV[rnn];    /* Tau integration variable, axis on which integral  */

  /* Integrate for each of the planet's layer starting from the           
     outermost until the closest layer. The order is oposite since tau starts from the top and rad arrya starts from the bottom  */  
  /*JB B_freq = \frac{2*PlankConst*freq^3}{c^2} * \frac{1}{exp^(\frac{PlankConst*freq}{Kb*T} - 1}    */
  /*JB B_wn = \frac{2*PlankConst*wn^3*1e8*c2} * \frac{1}{exp^(\frac{PlankConst*wn*100*c}{Kb*T} - 1}    */
                           
  for(i=0;i<=last; i++){   
    tauIV[i] = tau[i];  /* JB do I need radius factor here since it is in radius units, I guess not because it takes this array from tau.c */ 
    B[i] =  2. * H * w * w * w * wfct * wfct * wfct * LS*LS / ( exp(H * w * wfct * LS / (KB * temp[rnn-1-i])) - 1.);         
    tauInteg[i] = B[i] * exp(-tau[i]);
  }


  /* JB this is added as a 0 at the end when tau reach toomuch, so the spline looks nice */
  /* Add all other layers to be 0. Only two to have a nice ending
    spline and not unnecessary values: */
  for(; i<rnn; i++){
    tauInteg[i] = 0;
    tauIV[i]= tauIV[i-1] + 1; /* JB: Use geometic progression to have enough elements for integral to work, will have same result. */
   }
  /* Increment last to represent the number of elements, check that we
     have enough: */
  /* Adding additional 0 layer, and last represent number of elements -1
                so we need to add one more. */
 last+=2;  
    /* JB: If atmosphere is transparent, and at level tau has not reached tau.toomuch, last is set to max number of layers (rnn, instead of rnn-1, because we added 2 on the previous step), we never should go over it. */
  if(last > rnn)    
    last= rnn;

  /* JB check if we have at lkeast 3 layes to do the spline */
  if(last<3)
    transiterror(TERR_CRITICAL,
                 "Condition failed, less than 3 items (only %i) for radial "
                 "integration.\n", last);

//  printf("Pre integracije u INTENSITIJU: %li\n", rnn);

 
//  printf("Last == %li, rnn == %li\n", last, rnn);
  for(i=0; i<rnn; i++)
//    if(B[i]<0)
    {
      printf("JB current w=%9.4g and current temperature T=%9.4g and current tau[%li]=%9.4g\n and exp(-tau[%li])=%9.4g\n and B[%li]=%9.4g\n", w, temp[rnn-1-i], i, tau[i], i, exp(-tau[i]), i, B[i]);
//      printf(" Jasmina LS=%9.4g, H=%9.4g, KB=%9.4g\n", LS, H, KB);
      printf(" (2. * H * w * w * w * wfct * wfct * wfct * LS*LS)=%9.4g\n", (2. * H * w * w * w * wfct * wfct * wfct * LS*LS)) ;
      printf(" ( exp(H * w * wfct * LS / (KB * temp[rnn-1-i])) - 1.)=%9.4g\n", ( exp(H * w * wfct * LS / (KB * temp[rnn-1-i])) - 1.));
      printf("tauIV[%li]=%9.4g tauInteg[%li]=%9.4g\n", i, tauIV[i], i, tauInteg[i]);
    }
  /* Integrate along tau up to tau = toomuch: */                  
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl       = gsl_interp_alloc(gsl_interp_cspline, last);
  gsl_interp_init(spl, tauIV, tauInteg, last);
  res = gsl_interp_eval_integ(spl, tauIV, tauInteg,
                               tauIV[0], tauIV[last-1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free (acc);
#else
# error computation of modulation() without GSL is not implemented
#endif


//  printf("Posle integracije u INTENSITIJU: %li\n", rnn);


  return res;
}






/* \fcnfh
   Calculate the emeregent intensity at each wavenumber

   Return: emergent intensity  */

int
emergent_intens(struct transit *tr){ /* Main structure */
  static struct outputray st_out;
  tr->ds.out = &st_out;


  /* Initial variables:  */
  long w;
  prop_samp *rad = &tr->rads;     
  prop_samp *wn = &tr->wns;       
  transit_ray_solution *sol = tr->sol;


  /* Allocate the modulation array: */
  PREC_RES *out = st_out.o = (PREC_RES *)calloc(wn->n, sizeof(PREC_RES));
  struct optdepth *tau = tr->ds.tau;


  /* Integrate for each wavelength: */
  transitprint(1, verblevel, "\nIntegrating over wavelength.\n");

  /* This is printing process variable: */
  int nextw = wn->n/10;

  /* Calculate the intensity integral at each wavenumber: */
  for(w=0; w<wn->n; w++){
    out[w] = sol->eclIntenWn(tr, tau->t[w], wn->v[w], tau->last[w], tau->toomuch, rad);  /* JB do I need to convert this to units?  wn->v[w]*wn->fct: */


    /* Print to screen the progress status: */
    if(w==nextw){
      nextw += wn->n/10;
      transitprint(2, verblevel, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }
  transitprint(1, verblevel, "\nDone.\n");

  /* Free no longer needed memory: */
  freemem_idexrefrac(tr->ds.ir, &tr->pi);
  freemem_extinction(tr->ds.ex, &tr->pi);
  freemem_tau(tr->ds.tau,       &tr->pi);

  /* Set progress indicator, and print output: */
  tr->pi &= TRPI_MODULATION;
  printmod(tr);  
  return 0;
}


/* To the previous slantpath transit_ray_solution I just added my new function */
const transit_ray_solution slantpath = {
       "Slant Path",     /* Name of the solution                     */
       "slantpath.c",    /* Source code file name                    */
       1,                /* Equispaced impact parameter requested?   */
       NULL,        /* Per impact parameter and per wavenumber 
                            value computation                        */
/*       &totaltau,         Per impact parameter and per wavenumber 
                            value computation                        */
       &totaltau_eclipse,      /* Per per wavenumber value computation, see line 95 JB added*/
/*     &modulationperwn,              Per wavenumber value computation         */
       NULL,             /* Per wavenumber value computation         */
       &eclipse_intens,  /* Eclipse intensity per wavenumber value computation, see line 102 JB added */
       1                 /* Number of levels of modulation detail as
                            it can be handled by the above fucntion  */
       };



