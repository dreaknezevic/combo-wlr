/******************************************************************************
*   Program:  combo_wlr.sas
*   Author:   Andrea Knezevic, Memorial Sloan Kettering Cancer Center
*  				+ Pete Laud, SSU, University of Sheffield / AstraZeneca
*   Created:  February 2020 (SASGF 5062-2020)
*   Updated:  March 2021
*   Purpose:  Used to perform a combination weighted log-rank test
*				(with optional stratification)
*
*   Parameters (required):
*      data=    input dataset
*      group=   variable for assigned group (character format required)
*      time=    variable for observation time until event or censor
*      event=   variable indicator for event=1 or censor=0
*      estim=   method to be used for p-value estimation: 
*               ("GB" = Genz-Bretz, "sim" = simulation)
*
*   Parameters (optional):
*      strata=   stratification variable(s)
*      weights=  list of G(rho,gamma) weights, where rho>0, gamma>0
*                passed as: %str(rho1,gamma1 rho2,gamma2 ... )
*                default combination is: G(0,0), G(1,0), G(0,1), G(1,1)
*.     abseps=   Genz-Bretz algorithm precision paramater
*      maxpts=   Genz-Bretz algorithm precision paramater
*.     nsim=	 number of simulations to use for p-value estimate
*
*   Output datasets:
*      counts    Number of events per group (to check for over-stratification)
*      FH        Results of the selected weighted log-rank tests
*.     corr      Correlation matrix
*.     maxcombo	 MaxCombo test p-value
*
*   EXAMPLES: %combo_wlr(data = patients, 
*                        group = treatment, 
*                        time = months, 
*                        event = death, 
*                        weights = %str(0,0 0,1 1,1));
*
*             %combo_wlr(data = sashelp.bmt(where=(group in ("ALL","AML-Low Risk"))),
*                        group = group,
*                        time = T,
*                        event = status,
*                        weights = %str(0,0 2,0 0,2 2,2));
*
*             %combo_wlr(data = patients, 
*                        group = treatment, 
*                        time = months, 
*                        event = death, 
*                        weights = %str(0,0 0,1 1,1),
*                        strata = severity agegroup);
*
*   Acknowledgements:
*      Code for computing the nearest correlation matrix from The DO Loop blog by Rick Wicklin: 
*      https://blogs.sas.com/content/iml/2012/11/28/computing-the-nearest-correlation-matrix.html
*.     Code for Genz-Bretz algorithm for estimating the multivariate normal probability 
*      by Katherine Cai adapted from Genz and Bretz, from:
*      http://www.public.asu.edu/~jeffreyw/GMM/%25MVINTEGRATION%20Macro.sas
*
******************************************************************************/

%macro combo_wlr(data=, group=, time=, event=, weights=, strata=, 
					estim="GB", abseps=.00000005, maxpts=100000, nsim=50000000); 

options nosource nonotes nofmterr validvarname=V7;
ods exclude all;
%let _timer_start = %sysfunc(datetime());

data one; 
 set &data;
 keep &strata. &group &time &event; 
run;

*******************************************************;
*** Checks ;

* Check that data and variables exist ;
%let dsid = %sysfunc(open(one));
%if &dsid %then %do;
  %if %sysfunc(varnum(&dsid,&group))=0 |
      %sysfunc(varnum(&dsid,&time))=0 | 
      %sysfunc(varnum(&dsid,&event))=0 %then %do;
	    %let rc=%sysfunc(close(&dsid));
        %put ERROR: One or more variables do not exist.;
		ods exclude none; options source notes;
        %return; %end; %end;
%else %do;
  %let rc=%sysfunc(close(&dsid));
  %put ERROR: Data set does not exist.;
  ods exclude none; options source notes;
  %return; 
%end; 
%let rc=%sysfunc(close(&dsid));

* Assign default weights ;
%if &weights= %then %let weights=%str(0,0 1,0 0,1 1,1);

* Check if weights are correctly entered ;
%let count=%sysfunc(countw(&weights,%str( )));

%do i=1 %to &count;
 %let r = %scan(%qscan(&weights,&i,%str( )),1,%str(,));
 %let g = %scan(%qscan(&weights,&i,%str( )),2,%str(,));

 %if &r<0|&g<0 %then %do;
  %put ERROR: Weights are not entered correctly.;
  ods exclude none; options source notes;
  %return; %end;
%end;

* Check for duplicate weights ;
%let weights2 = %qscan(&weights,1,%str( ));
%do i=2 %to &count;
  %let w = %qscan(&weights,&i,%str( ));
  %if %sysfunc(find(&weights2,&w))=0 %then %let weights2 = &weights2 &w;
%end;

%if &weights ~= &weights2 %then %do;
 %put WARNING: Duplicate sets of weights will be used only once.;
 %let weights = &weights2;
 %let count=%sysfunc(countw(&weights,%str( )));
%end;

* Check to make sure group variable has only 2 levels ;
ods output nlevels=nlevels;
proc freq data=one nlevels;
 tables &group;
run; 

proc sql; 
 select nlevels into:nlevels from nlevels;
quit;

%if &nlevels~=2 %then %do;
  %put ERROR: Group variable must have 2 levels.;
  ods exclude none; options source notes;
%return; %end;

*******************************************************;
*******************************************************;
%put Calculating weighted log-rank statistics...;

data fh; stop; 
 length test $50;
run;

%do i=1 %to &count;

 %let r = %scan(%qscan(&weights,&i,%str( )),1,%str(,));
 %let g = %scan(%qscan(&weights,&i,%str( )),2,%str(,));

 ods output HomTests=chisq HomStats=Q FlemingHomCov=var;
 proc lifetest data=one;
  time &time*&event(0);
  strata &strata. / test=FH(&r,&g) group=&group.;
 run;

 proc contents data=var(drop=&group) out=name(keep=name); run;
 data _null_; set name(obs=1); call symput('name',name); run;

 data stat; merge
  chisq q(obs=1) var(obs=1);
  var=&name;
  z=fleming/sqrt(var);
  keep test var z probchisq;
 run;

 %if &count=1 %then %do;
  %put WARNING: Only one set of weights given. No combination test performed.;
  data fh; set stat; run;
  %goto exit; 
 %end;

 data fh; set fh stat; run;

%end;

*******************************************************;
%put Calculating correlation matrix...;

data covar; stop; 
 length test $50; 
run;

* Put weights in data set ;
data weights; retain i;
length w $50;
%do i=1 %to &count;
  i=&i;
  w="%qscan(&weights,&i,%str( ))";
  output;
%end;
run;

* List out combination indices and merge weights back in ; 
proc iml;
idx = allcomb(&count, 2);
create comb from idx[colname={i1 i2}];
 append from idx;
close comb;
quit;

proc sort data=comb; by i2; run;

data comb; merge 
 comb(in=a rename=(i2=i)) 
 weights(rename=(w=w2));
  by i;
  if a;
run;

proc sort data=comb; by i1; run;

data comb; merge 
 comb(in=a)
 weights(rename=(i=i1 w=w1));
  by i1;
  if a;
run;

data comb; set comb;
 r1 = scan(w1,1,',');
 g1 = scan(w1,2,',');
 r2 = scan(w2,1,',');
 g2 = scan(w2,2,',');
 r = (r1+r2)/2;
 g = (g1+g2)/2;
 keep r g; 
run;

proc sql noprint; 
 select count(*) into:ncomb from comb; 
quit;

%do i=1 %to &ncomb;

 data _null_; set comb(obs=&i);
  call symput("r",r);
  call symput("g",g);
 run; 

 ods output HomTests=chisq FlemingHomCov=var;
 proc lifetest data=one;
  time &time*&event(0);
  strata &strata. / test=FH(&r,&g) group=&group.; 
 run;

 data stat; merge
  chisq var(obs=1);
  var=&name;
  keep test var;
 run;

 data covar; set covar stat; run;

%end;

*******************************************************;

proc iml;
 **************************************************************;
 /* Project symmetric X onto S={positive semidefinite matrices}.
    Replace any negative eigenvalues of X with zero */
 start ProjS(X);
   call eigen(D, Q, X); /* note X = Q*D*Q` */
   V = choose(D>.0001, D, .0001);
   W = Q#sqrt(V`);      /* form Q*diag(V)*Q` */
   return( W*W` );      /* W*W` = Q*diag(V)*Q` */
 finish;
 
 /* project square X onto U={matrices with unit diagonal}.
    Return X with the diagonal elements replaced by ones. */
 start ProjU(X);
   n = nrow(X);
   Y = X;
   diagIdx = do(1, n*n, n+1);
   Y[diagIdx] = 1;      /* set diagonal elements to 1 */
   return ( Y );
 finish;
 
 /* Helper function: the infinity norm is the max abs value of the row sums */
 start MatInfNorm(A);
   return( max(abs(A[,+])) );
 finish;
 
 /* Given a symmetric correlation matrix, A, 
    project A onto the space of positive semidefinite matrices.
    The function uses the algorithm of Higham (2002) to return
    the matrix X that is closest to A in the Frobenius norm. */
 start NearestCorr(A);
   maxIter = 100; tol = 1e-8;  /* parameters...you can change these */
   iter = 1;      maxd   = 1;  /* initial values */ 
   Yold = A;  Xold = A;  dS = 0;
 
   do while( (iter <= maxIter) & (maxd > tol) );
     R = Yold - dS; /* dS is Dykstra's correction */
     X = ProjS(R);  /* project onto S={positive semidefinite} */
     dS = X - R;
     Y = ProjU(X);  /* project onto U={matrices with unit diagonal} */
 
     /* How much has X changed? (Eqn 4.1) */
     dx = MatInfNorm(X-Xold) / MatInfNorm(X);
     dy = MatInfNorm(Y-Yold) / MatInfNorm(Y);
     dxy = MatInfNorm(Y - X) / MatInfNorm(Y);
     maxd = max(dx,dy,dxy);
 
     iter = iter + 1;
     Xold = X; Yold = Y; /* update matrices */
   end;
   return( X ); /* X is positive semidefinite */
 finish;

 /*Genz-Bretz algorithm functions and constants*/
   START MVN_DIST( N, LOWER, UPPER, INFIN, COVAR, MAXPTS, ABSEPS, RELEPS,   ERROR, VALUE, NEVALS, INFORM );
      NEVALS = 0;
      RUN MVNDNT( N, COVAR, LOWER, UPPER, INFIN,   INFIS, VALUE, ERROR, INFORM );
      IF ( INFORM = 0 ) THEN DO;
         IF ( N-INFIS > 2 ) THEN RUN DKBVRC( N-INFIS-1, 0, MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, NEVALS, INFORM );
      END;
   FINISH MVN_DIST;

   START MVN_DFN( N, W ) GLOBAL( COVARS, DONE, EONE, INFI, A, B, Y );
      VALUE = 1;
      INFA = 0;
      INFB = 0;
      IK = 1;
      DO I = 1 TO N+1;
         VSUM = 0;
         IF ( IK > 1 ) THEN VSUM = VSUM + SUM( COVARS[I,1:IK-1]#T(Y[1:IK-1]));
         IF ( INFI[I] ^= 0 ) THEN DO;
            IF ( INFA = 1 ) THEN AI = MAX( AI, A[I] - VSUM );
            ELSE DO;
               AI = A[I] - VSUM;
               INFA = 1;
            END;
         END;
         IF ( INFI[I] ^= 1 ) THEN DO;
            IF ( INFB = 1 ) THEN BI = MIN( BI, B[I] - VSUM );
            ELSE DO;
               BI = B[I] - VSUM;
               INFB = 1;
            END;
         END;
         IF ( I < N + 1 & IK < N + 1 ) THEN AAA = COVARS[I+1,IK+1];
         ELSE AAA = 0;
         IF ( I = N+1 | AAA > 0 ) THEN DO;
            IF ( I = 1 ) THEN DO;
               DI = DONE;
               EI = EONE;
            END;
            ELSE DO;
               DI = 0;
               EI = 1;
               J = 2*INFA+INFB-1;
               IF ( J >= 0 ) THEN DO;
                  IF ( J ^= 0 ) THEN DI = PROBNORM(AI);
                  IF ( J ^= 1 ) THEN EI = PROBNORM(BI);
               END;
               EI = MAX( EI, DI );
            END;
            IF ( DI >= EI ) THEN DO;
               VALUE = 0;
               I = N+1;
            END;
            ELSE DO;
               VALUE = VALUE*( EI - DI );
               IF ( I <= N ) THEN Y[IK] = PROBIT( DI + W[IK]*( EI - DI ) );
               IK = IK + 1;
               INFA = 0;
               INFB = 0;
            END;
         END;
      END;
      MVNDFN = VALUE;
      RETURN( MVNDFN );
   FINISH MVN_DFN;

   START MVNDNT( N, COVAR, LOWER, UPPER, INFIN,   INFIS, VALUE, ERROR, INFORM ) GLOBAL( COVARS, DONE, EONE, INFI, A, B, NL );
      INFORM = 0;
      IF ( N > NL | N < 1 ) THEN INFORM = 2;
      ELSE DO  I = 1 TO N;
         IF ( INFIN[I] > 2 ) THEN INFORM = 3;
         ELSE IF ( INFIN[I] = 2 & LOWER[I] > UPPER[I] ) THEN INFORM = 3;
      END;
      IF ( INFORM = 0 ) THEN RUN COVSRT( N, LOWER, UPPER, COVAR, INFIN, INFIS, INFORM );
      IF ( INFORM = 0 ) THEN DO;
         IF ( N - INFIS = 0 ) THEN DO;
            VALUE = 1;
            ERROR = 0;
         END;
         ELSE DO;
            IF ( N - INFIS = 1 ) THEN DO;
               VALUE = EONE - DONE;
               ERROR = 2E-15;
            END;
            ELSE DO;
               IF ( N - INFIS = 2 ) THEN DO;
                  IF ( ABS( COVARS[2,2] ) > 0 ) THEN DO;
                     D = SQRT( 1 + COVARS[2,1]**2 );
                     IF ( INFI[2] ^= 0 ) THEN A[2] = A[2]/D;
                     IF ( INFI[2] ^= 1 ) THEN B[2] = B[2]/D;
                     VALUE = PROBBVN( A, B, INFI, COVARS[2,1]/D );
                  END;
                  ELSE DO;
                     IF ( INFI[1] ^= 0 ) THEN DO;
                        IF ( INFI[2] ^= 0 ) THEN A[1] = MAX( A[1], A[2] );
                     END;
                     ELSE DO;
                        IF ( INFI[2] ^= 0 ) THEN A[1] = A[2];
                     END;

                     IF ( INFI[1] ^= 1 ) THEN DO;
                        IF ( INFI[2] ^= 1 ) THEN B[1] = MIN( B[1], B[2] );
                     END;
                     ELSE DO;
                        IF ( INFI[2] ^= 1 ) THEN B[1] = B[2];
                     END;

                     IF ( INFI[1] ^= INFI[2] ) THEN INFI[1] = 2;
                     RUN MVNLMS( A[1], B[1], INFI[1],   D, E );
                     VALUE = E - D;
                  END;
                  ERROR = 2E-15;
               END;
               ELSE DO;
                  VALUE = 0;
                  ERROR = 1;
               END;
            END;
         END;
      END;
      ELSE DO;
         VALUE = 0;
         ERROR = 1;
      END;
   FINISH MVNDNT;

   START MVNLMS( A, B, INFIN,   LOWER, UPPER );
      LOWER = 0;
      UPPER = 1;
      IF ( INFIN >= 0 ) THEN DO;
         IF ( INFIN ^= 0 ) THEN LOWER = PROBNORM(A);
         IF ( INFIN ^= 1 ) THEN UPPER = PROBNORM(B);
      END;
      UPPER = MAX( UPPER, LOWER );
   FINISH MVNLMS;

   START COVSRT( N, LOWER, UPPER, COVAR, INFIN,   INFIS, INFORM )
      GLOBAL( EPS, SQTWPI, COVARS, DONE, EONE, INFI, A, B, Y );

      Y = J( 1, N, . );
      INFI = INFIN;
      A = J( 1, N, 0 );
      B = J( 1, N, 0 );
      COVARS = COVAR;
      INFIS = N - SUM( SIGN( SIGN( INFI ) + 1 ) );
      DO I = 1 TO N;
         IF ( INFI[I] >= 0 ) THEN DO;
            IF ( INFI[I] ^= 0 ) THEN A[I] = LOWER[I];
            IF ( INFI[I] ^= 1 ) THEN B[I] = UPPER[I];
         END;
      END;

      IF ( INFIS < N ) THEN DO;
         DO I = N TO N-INFIS+1 BY -1;
            IF ( INFI[I] >= 0 ) THEN DO;
               DO J = 1 TO I-1;
                  IF ( INFI[J] < 0 ) THEN DO;
                     RUN RCSWP( J, I, A, B, INFI, N, COVARS );
                     J = I-1;
                  END;
               END;
            END;
         END;

         DO I = 1 TO N-INFIS;
            DEMIN = 1;
            JMIN = I;
            CVDIAG = 0;
            EPSI = I*I*EPS;
            DO J = I TO N-INFIS;
               IF ( COVARS[J,J] > EPSI ) THEN DO;
                  SUMSQ = SQRT( COVARS[J,J] );
                  VSUM = 0;
                  IF ( I > 1 ) THEN VSUM = SUM( COVARS[J,1:I-1] # T(Y[1:I-1]) );
                  AJ = ( A[J] - VSUM )/SUMSQ;
                  BJ = ( B[J] - VSUM )/SUMSQ;
                  RUN MVNLMS( AJ, BJ, INFI[J],   DD, EE );
                  IF ( DEMIN >= EE - DD ) THEN DO;
                     JMIN = J;
                     AMIN = AJ;
                     BMIN = BJ;
                     DEMIN = EE - DD;
                     CVDIAG = SUMSQ;
                  END;
               END;
            END;

            IF ( JMIN > I ) THEN RUN RCSWP( I, JMIN, A, B, INFI, N, COVARS );

            IF ( CVDIAG > 0 ) THEN DO;
               COVARS[I,I] = CVDIAG;
               DO L = I+1 TO N-INFIS;
                  COVARS[L,I] = COVARS[L,I]/CVDIAG;
                  COVARS[L,I+1:L] = COVARS[L,I+1:L] - COVARS[L,I] # T(COVARS[I+1:L,I]);
               END;

               IF ( DEMIN > EPSI ) THEN DO;
                  YL = 0;
                  YU = 0;
                  IF ( INFI[I] ^= 0 ) THEN YL = -EXP( -AMIN**2/2 )/SQTWPI;
                  IF ( INFI[I] ^= 1 ) THEN YU = -EXP( -BMIN**2/2 )/SQTWPI;
                  Y[I] = ( YU - YL )/DEMIN;
               END;
               ELSE DO;
                  IF ( INFI[I] = 0 ) THEN Y[I] = BMIN;
                  IF ( INFI[I] = 1 ) THEN Y[I] = AMIN;
                  IF ( INFI[I] = 2 ) THEN Y[I] = ( AMIN + BMIN )/2;
               END;

               COVARS[I,1:I] = COVARS[I,1:I]/CVDIAG;
               A[I] = A[I]/CVDIAG;
               B[I] = B[I]/CVDIAG;
            END;
            ELSE DO;
               IF ( COVARS[I,I] > -EPSI ) THEN DO;
                  COVARS[I:N-INFIS,I] = 0;

                  AAA = 0;
                  DO J = I-1 TO 1 BY -1;
                     IF ( ABS( COVARS[I,J] ) > EPSI ) THEN DO;
                        A[I] = A[I]/COVARS[I,J];
                        B[I] = B[I]/COVARS[I,J];
                        IF ( COVARS[I,J] < 0 ) THEN DO;
                           AA = A[I];
                           A[I] = B[I];
                           B[I] = AA;
                           IF ( INFI[I] ^= 2 ) THEN INFI[I] = 1 - INFI[I];
                        END;
                        COVARS[I,1:J] = COVARS[I,1:J]/COVARS[I,J];
                        DO L = J+1 TO I-1;
                           IF( COVARS[L,J+1] > 0 ) THEN DO;
                              DO K = I-1 TO L BY -1;

                                 AA = COVARS[K,1:K];
                                 COVARS[K,1:K] = COVARS[K+1,1:K];
                                 COVARS[K+1,1:K] = AA;

                                 AA = A[K];
                                 A[K] = A[K+1];
                                 A[K+1] = AA;

                                 AA = B[K];
                                 B[K] = B[K+1];
                                 B[K+1] = AA;

                                 M = INFI[K];
                                 INFI[K] = INFI[K+1];
                                 INFI[K+1] = M;
                              END;
                              L = I-1;
                           END;
                        END;
                        J = 1;
                        AAA = 1;
                     END;
                     IF AAA = 1 THEN;
                     ELSE COVARS[I,J] = 0;
                  END;
                  Y[I] = 0;
               END;
               ELSE DO;
                 INFORM = 4;
                 I = N-INFIS;
               END;
            END;
         END;
         IF (INFORM = 0 ) THEN RUN MVNLMS( A[1], B[1], INFI[1], DONE, EONE );
      END;
   FINISH COVSRT;

   START RCSWP( P, Q, A, B, INFIN, N, C );
      AA = A[P];
      A[P] = A[Q];
      A[Q] = AA;

      AA = B[P];
      B[P] = B[Q];
      B[Q] = AA;

      I = INFIN[P];
      INFIN[P] = INFIN[Q];
      INFIN[Q] = I;

      AA = C[P,P];
      C[P,P] = C[Q,Q];
      C[Q,Q] = AA;

      IF (P>1) THEN DO;
         AA = C[Q,1:P-1];
         C[Q,1:P-1] = C[P,1:P-1];
         C[P,1:P-1] = AA;
      END;

      DO I = P+1 TO Q-1;
         AA = C[I,P];
         C[I,P] = C[Q,I];
         C[Q,I] = AA;
      END;

      IF (Q<N) THEN DO;
         AA = C[Q+1:N,P];
         C[Q+1:N,P] = C[Q+1:N,Q];
         C[Q+1:N,Q] = AA;
      END;
   FINISH RCSWP;

   START DKBVRC( NDIM, MINVLS, MAXVLS, ABSEPS, RELEPS,   ABSERR, FINEST, INTVLS, INFORM )
      GLOBAL( PLIM, KLIM, P, C, MINSMP );

      VK     = J( 1, KLIM, . );
      INFORM = 1;
      INTVLS = 0;
      KLIMI  = KLIM;
      IF ( MINVLS >= 0 ) THEN DO;
         FINEST = 0;
         VAREST = 0;
         SAMPLS = MINSMP;
         DO I = 1 TO PLIM;
            NP = I;
            IF ( MINVLS < 2*SAMPLS*P[I] ) THEN I = PLIM;
         END;
         IF ( MINVLS >= 2*SAMPLS*P[PLIM] ) THEN SAMPLS = MINVLS/( 2*P[PLIM] );
      END;
      VALUE = J( 1, SAMPLS, . );
      EXIT  = 0;
      DO UNTIL( EXIT = 1);
         VK[1] = 1/P[NP];
         DO I = 2 TO MIN( NDIM, KLIM );
            VK[I] = MOD( C[NP, MIN(NDIM-1,KLIM-1)]*VK[I-1], 1 );
         END;
         FINVAL = 0;
         VARSQR = 0;
         DO I = 1 TO SAMPLS;
            VALUE[I] = DKSMRC( NDIM, KLIMI, P[NP], VK );
         END;
         FINVAL = VALUE[:];
         VARSQR = (VALUE[##] - VALUE[+]##2/SAMPLS) / (SAMPLS # (SAMPLS-1));
         INTVLS = INTVLS + 2*SAMPLS*P[NP];
         VARPRD = VAREST*VARSQR;
         FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD );
         IF ( VARSQR > 0 ) THEN VAREST = ( 1 + VARPRD )/VARSQR;
         ABSERR = 3*SQRT( VARSQR/( 1 + VARPRD ) );
         IF ( ABSERR > MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN DO;
            IF ( NP < PLIM ) THEN NP = NP + 1;
            ELSE DO;
               SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P[NP] ) );
               SAMPLS = MAX( MINSMP, SAMPLS );
            END;
            IF ( INTVLS + 2*SAMPLS*P[NP] > MAXVLS ) THEN EXIT = 1;
         END;
         ELSE DO;
            INFORM = 0;
            EXIT = 1;
         END;
      END;
   FINISH DKBVRC;

   START DKSMRC( NDIM, KLIM, PRIME, VK );
      X = J( 1, NDIM, . );
      SUMKRO = 0;
      NK = MIN( NDIM, KLIM );
      DO J = 1 TO NK-1;
         JP = J + RANUNI(141071)*( NK + 1 - J );
         XT = VK[J];
         VK[J] = VK[INT(JP)];
         VK[INT(JP)] = XT;
      END;
      XP = RANUNI( J(1, NDIM, 141071) );
      DIFF = NDIM - KLIM;
      DO K = 1 TO PRIME;
         IF ( DIFF>0 ) THEN DO;
            RUN DKRCHT( DIFF, X );
            DO J = 1 TO NDIM-KLIM;
               X[NK+J] = X[J];
            END;
         END;
         X[1:NK] = MOD( K*VK[1:NK], 1 );
         DO J = 1 TO NDIM;
            XT = X[J] + XP[J];
            IF ( XT > 1 ) THEN  XT = XT - 1;
            X[J] = ABS( 2*XT - 1 );
         END;
         MVNDFN = MVN_DFN( NDIM, X );
         SUMKRO = SUMKRO + ( MVNDFN - SUMKRO )/( 2*K - 1 );
         X = 1 - X;
         MVNDFN = MVN_DFN( NDIM, X );
         SUMKRO = SUMKRO + ( MVNDFN - SUMKRO )/( 2*K );
      END;
      RETURN (SUMKRO);
   FINISH DKSMRC;

   START DKRCHT( S, QUASI ) GLOBAL( NN, PSQT, HISUM, OLDS, MXDIM, MXHSUM, BB );

      IF ( S ^= OLDS | S < 1 ) THEN DO;
         OLDS = S;
         NN[1] = 0;
         HISUM = 0;
      END;

      I = 0;
      CRIT = 0;
      DO UNTIL( CRIT = 1 | I = HISUM + 1 );
         NN[I + 1] = NN[I + 1] + 1;
         IF ( NN[I + 1] < BB ) THEN DO;
           CRIT = 1;
           I = I - 1;
         END;
         ELSE NN[I + 1] = 0;
         I = I + 1;
      END;

      IF ( I > HISUM ) THEN DO;
         HISUM = HISUM + 1;
         IF ( HISUM > MXHSUM ) THEN HISUM = 0;
         NN[HISUM + 1] = 1;
      END;

      RN = 0;
      DO I = HISUM TO 0 BY -1;
         RN = NN[I + 1] + BB * RN;
      END;
      QUASI[1:S] = MOD( RN # PSQT[1:S], 1 );
   FINISH DKRCHT;

   START PROBBVN( LOWER, UPPER, INFIN, CORREL );
      IF ( INFIN[1] = 2  & INFIN[2] = 2 ) THEN
         BVN =  PROBBNRM ( LOWER[1], LOWER[2], CORREL ) - PROBBNRM ( UPPER[1], LOWER[2], CORREL )
              - PROBBNRM ( LOWER[1], UPPER[2], CORREL ) + PROBBNRM ( UPPER[1], UPPER[2], CORREL );
      ELSE IF ( INFIN[1] = 2  & INFIN[2] = 1 ) THEN
         BVN =  PROBBNRM ( -LOWER[1], -LOWER[2], CORREL ) - PROBBNRM ( -UPPER[1], -LOWER[2], CORREL );
      ELSE IF ( INFIN[1] = 1  & INFIN[2] = 2 ) THEN
         BVN =  PROBBNRM ( -LOWER[1], -LOWER[2], CORREL ) - PROBBNRM ( -LOWER[1], -UPPER[2], CORREL );
      ELSE IF ( INFIN[1] = 2  & INFIN[2] = 0 ) THEN
         BVN =  PROBBNRM ( UPPER[1], UPPER[2], CORREL ) - PROBBNRM ( LOWER[1], UPPER[2], CORREL );
      ELSE IF ( INFIN[1] = 0  & INFIN[2] = 2 ) THEN
         BVN =  PROBBNRM ( UPPER[1], UPPER[2], CORREL ) - PROBBNRM ( UPPER[1], LOWER[2], CORREL );
      ELSE IF ( INFIN[1] = 1  & INFIN[2] = 0 ) THEN BVN = PROBBNRM ( -LOWER[1],  UPPER[2], -CORREL );
      ELSE IF ( INFIN[1] = 0  & INFIN[2] = 1 ) THEN BVN = PROBBNRM (  UPPER[1], -LOWER[2], -CORREL );
      ELSE IF ( INFIN[1] = 1  & INFIN[2] = 1 ) THEN BVN = PROBBNRM ( -LOWER[1], -LOWER[2],  CORREL );
      ELSE IF ( INFIN[1] = 0  & INFIN[2] = 0 ) THEN BVN = PROBBNRM (  UPPER[1],  UPPER[2],  CORREL );
      RETURN ( BVN );
   FINISH PROBBVN;

   HISUM = .;
   OLDS  = 0;
   MXDIM = 80;
   MXHSUM = 50;
   BB = 2;
   PSQT={1.414213562373 1.732050807569 2.236067977500 2.645751311065 3.316624790355 3.605551275464
      4.123105625618 4.358898943541 4.795831523313 5.385164807135 5.567764362830 6.082762530298
      6.403124237433 6.557438524302 6.855654600401 7.280109889281 7.681145747869 7.810249675907
      8.185352771872 8.426149773176 8.544003745318 8.888194417316 9.110433579144 9.433981132057
      9.848857801796 10.04987562112 10.14889156509 10.34408043279 10.44030650891 10.63014581273
      11.26942766958 11.44552314226 11.70469991072 11.78982612255 12.20655561573 12.28820572744
      12.52996408614 12.76714533480 12.92284798332 13.15294643797 13.37908816026 13.45362404707
      13.82027496109 13.89244398945 14.03566884762 14.10673597967 14.52583904633 14.93318452307
      15.06651917332 15.13274595042 15.26433752247 15.45962483374 15.52417469626 15.84297951775
      16.03121954188 16.21727474023 16.40121946686 16.46207763315 16.64331697709 16.76305461424
      16.82260384126 17.11724276862 17.52141546794 17.63519208855 17.69180601295 17.80449381476
      18.19340539866 18.35755975069 18.62793601020 18.68154169227 18.78829422806 18.94729532150
      19.15724406067 19.31320791583 19.46792233393 19.57038579078 19.72308292332 19.92485884517
      20.02498439450 20.22374841616};

   NL = 100;
   EPS = 1E-10;
   SQTWPI = 2.506628274631000502415765284811045253;

   PLIM = 25;
   KLIM = 20;
   MINSMP = 8;
   P = { 31 47 73 113 173 263 397 593 907 1361 2053 3079 4621 6947 10427 15641
      23473 35221 52837 79259 118891 178349 267523 401287 601942};
   C = { 12 9 9 13 12 12 12 12 12 12 12 12 3 3 3 12 7 7 12,
      13 11 17 10 15 15 15 15 15 15 22 15 15 6 6 6 15 15 9 ,
      27 28 10 11 11 20 11 11 28 13 13 28 13 13 13 14 14 14 14 ,
      35 27 27 36 22 29 29 20 45 5 5 5 21 21 21 21 21 21 21 ,
      64 66 28 28 44 44 55 67 10 10 10 10 10 10 38 38 10 10 10 ,
      111 42 54 118 20 31 31 72 17 94 14 14 11 14 14 14 94 10 10 ,
      163 154  83 43 82 92 150 59 76 76 47 11 11 100 131 116 116 116 116 ,
      246 189 242 102 250 250 102 250 280 118 196 118 191 215 121 121 49 49 49 ,
      347 402 322 418 215 220 339 339 339 337 218 315 315 315 315 167 167 167 167 ,
      505 220 601 644 612 160 206 206 206 422 134 518 134 134 518 652 382 206 158 ,
      794 325 960 528 247 247 338 366 847 753 753 236 334 334 461 711 652 381 381 ,
      1189 888 259 1082 725 811 636 965 497 497 1490 1490 392 1291 508 508 1291 1291 508 ,
      1763 1018 1500 432 1332 2203 126 2240 1719 1284 878 1983 266 266 266 266 747 747 127 ,
      2872 3233 1534 2941 2910 393 1796 919 446 919 919 1117 103 103 103 103 103 103 103 ,
      4309 3758 4034 1963 730 642 1502 2246 3834 1511 1102 1102 1522 1522 3427 3427 3928 915 915 ,
      6610 6977 1686 3819 2314 5647 3953 3614 5115 423 423 5408 7426 423 423 487 6227 2660 6227 ,
      9861 3647 4073 2535 3430 9865 2830 9328 4320 5913 10365 8272 3706 6186 7806 7806 7806 8610 2563 ,
      10327 7582 7124 8214 9600 10271 10193 10800 9086 2365 4409 13812 5661 9344 9344 10362 9344 9344 8585 ,
      19540 19926 11582 11113 24585 8726 17218 419 4918 4918 4918 15701 17710 4037 4037 15808 11401 19398 25950 ,
      34566 9579 12654 26856 37873 38806 29501 17271 3663 10763 18955 1298 26560 17132 17132 4753 4753 8713 18624 ,
      31929  49367 10982 3527 27066 13226 56010 18911 40574 20767 20767 9686 47603 47603 11736 11736 41601 12888 32948 ,
      40701  69087 77576 64590 39397 33179 10858 38935 43129 35468 35468 2196 61518 61518 27945 70975 70975 86478 86478 ,
      103650 125480 59978 46875 77172 83021 126904 14541 56299 43636 11655 52680 88549 29804 101894 113675 48040 113675 34987 ,
      165843 90647 59925 189541 67647 74795 68365 167485 143918 74912 167289 75517 8148 172106 126159 35867 35867 35867 121694 ,
      130365 236711 110235 125699 56483 93735 234469 60549 1291 93937 245291 196061 258647 162489 176631 204895 73353 172319 28881};

 **************************************************************;

 use fh; read all var{"z"} into z;
         read all var{"var"} into var; close fh;
 use covar; read all var{"var"} into covar; close covar;

 cov = j(nrow(var), nrow(var), 0);
 cov[loc(row(cov) < col(cov))] = covar;
 
 cov = cov + cov`;
 cov[loc(row(cov) = col(cov))] = var; 

 zmax = max(abs(z)); 
 corr = cov2corr(cov); print corr;

 create zmax from zmax; 
  append from zmax; 
 close zmax;

 create corr from corr;
  append from corr;
 close corr;

*******************************************************;
 %put Calculating p-value...;

 k = nrow(z);
 %if &estim. = "sim" %then %do;
   * Calculate nearest positive 
     semidefinite correlation matrix ;
   corr2 = nearestcorr(corr);
   nsim = &nsim.; 
   mean = repeat(0,k)`;
   call randseed(1);
   x = RANDNORMAL(nsim, mean, corr2);
   create test from x;
    append from x;
   close test;
 %end;
 %else %if &estim. = "GB" %then %do;
   LOWER = J(1,k,-zmax);
   UPPER = J(1,k,zmax);
   INFIN = J(1,k,2);
   MAXPTS = &maxpts.;
   ABSEPS = &abseps.;
   RELEPS = 0;
   RUN MVN_DIST( k, LOWER, UPPER, INFIN, corr, MAXPTS, ABSEPS, RELEPS,  ERROR, VALUE, NEVALS, INFORM );
   pval = (1 - value);
   create test var {"pval" "error" "nevals"};
    append;
   close test;
 %end;

quit;

proc sql noprint; 
 select col1 into:zmax from zmax; 
quit;

%if &estim. = "sim" %then %do;

  data test; set test;
   zmax = max(%do i=1 %to &count;
 			 %if &i ne &count %then abs(col&i), ;
 			  %else abs(col&i); %end; );
   if zmax > &zmax then x=1; 
   else x=0; 
  run; 

  proc sql; 
   create table maxcombo as
   select mean(x) as p label="pval",
		&zmax as zmax label="Z max", 
        &estim. as method label="Estimation method"
   from test;
  quit;
%end;
%else %if &estim.="GB" %then %do;
  proc sql; 
   create table maxcombo as
   select pval as p label="pval",
		&zmax as zmax label="Z max", 
        &estim. as method label="Estimation method",
		error,
		nevals
   from test;
  quit;
%end;

*******************************************************;
%exit:
*******************************************************;

ods exclude none;
ods noproctitle;
ods graphics on;

*ods html options(pagebreak='yes');
title h=11pt "Combination weighted log-rank tests";
ods select SurvivalPlot;
ods output CensoredSummary=counts;
proc lifetest data=one plots=s(atrisk cl);
  time &time*&event(0);
  strata &strata. / group=&group.;
run;

%if &strata. ne %then %do;
  title h=10pt "Number of events per stratum";
%end;
%else %do; 
  title h=10pt "Number of events";
%end;
proc print data=counts(where=(stratum > 0)) noobs;
 var &strata. &group. total failed;
run;

*ods html options(pagebreak='no');
title h=10pt "Weighted log-rank tests"; 
%if &strata. ne %then %do;
  title2 h=10pt "Stratified by &strata.";
%end;
proc print data=fh noobs label; 
 var test z probchisq; 
 label test="Test" z="Z statistic" probchisq="P"; 
run;

%if &count>1 %then %do;
  title h=10pt "Correlation between weighted log-rank tests"; 
  proc print data=corr noobs;* label; 
*   var test z probchisq; 
*   label test="Test" z="Z statistic" probchisq="P"; 
  run;

*ods html options(pagebreak='no');
  title h=10pt "Combination test";
  %if &strata. ne %then %do;
    title2 h=10pt "Stratified by &strata.";
  %end;
  proc print data=maxcombo noobs label; 
   format p pvalue10.5;
  run;
%end;

ods proctitle;
*ods html options(pagebreak='yes');
title;

*******************************************************;

* Clean up ;
proc datasets noprint;
 delete nlevels name one stat chisq q var covar
         test zmax weights comb;
run; quit;

* Print runtime to log ;
data _null_;
  dur = datetime() - &_timer_start;
  put 24*'-' / ' Run time:' dur time13.2 / 24*'-';
run;

options source notes;
ods graphics off;

%mend;
