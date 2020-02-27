/******************************************************************************
*   Program:  combo_wlr.sas
*   Author:   Andrea Knezevic, Memorial Sloan Kettering Cancer Center
*   Created:  February 2020 (SASGF 5062-2020)
*   Purpose:  Used to perform a combination weighted log-rank test
*
*   Parameters (required):
*      data=    input dataset
*      group=   variable for assigned group
*      time=    variable for observation time until event or censor
*      event=   variable indicator for event=1 or censor=0
*
*   Parameters (optional):
*      weights=  list of G(rho,gamma) weights, where rho>0, gamma>0
*                passed as: %str(rho1,gamma1 rho2,gamma2 ... )
*                default combination is: G(0,0), G(1,0), G(0,1), G(1,1)
*
*   EXAMPLES: %combo_wlr(data = patients, 
*                        group = treatment, 
*                        time = months, 
*                        event = death, 
*                        weights = %str(0,0 1,0 0,1));
*
*             %combo_wlr(data = sashelp.bmt(where=(group in ("ALL","AML-Low Risk"))),
*                        group = group,
*                        time = T,
*                        event = status,
*                        weights = %str(0,0 2,0 0,2 2,2));
*
*   Acknowledgement:
*      Code for computing the nearest correlation matrix from The DO Loop blog by Rick Wicklin: 
*      https://blogs.sas.com/content/iml/2012/11/28/computing-the-nearest-correlation-matrix.html
*
******************************************************************************/

%macro combo_wlr(data=,group=,time=,event=,weights=);

options nosource nonotes nofmterr;
ods exclude all;
%let _timer_start = %sysfunc(datetime());

data one; 
 set &data;
 keep &group &time &event; 
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
  strata &group / test=FH(&r,&g);
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
  strata &group / test=FH(&r,&g);
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

 * Calculate nearest positive 
   semidefinite correlation matrix ;
 corr2 = nearestcorr(corr);

*******************************************************;
 %put Calculating p-value...;

 n = 5000000; 
 k = nrow(z);
 mean = repeat(0,k)`;
 call randseed(1);
 x = RANDNORMAL(n, mean, corr2);

 create test from x;
  append from x;
 close test;

 create zmax from zmax; 
  append from zmax; 
 close zmax;
quit;

proc sql noprint; 
 select col1 into:zmax from zmax; 
quit;

data test; set test;
 zmax = max(%do i=1 %to &count;
 			 %if &i ne &count %then abs(col&i), ;
 			  %else abs(col&i); %end; );
 if zmax > &zmax then x=1; 
  else x=0; 
run; 

proc sql; 
 create table maxcombo as
 select &zmax as zmax label="Z max", 
        mean(x) as p label="P" from test;
quit;

*******************************************************;
%exit:
*******************************************************;

ods exclude none;
ods noproctitle;

ods html options(pagebreak='yes');
title h=11pt "Combination weighted log-rank tests";
ods select SurvivalPlot;
proc lifetest data=one plots=s(atrisk cl);
  time &time*&event(0);
  strata &group;
run;

ods html options(pagebreak='no');
title h=10pt "Weighted log-rank tests"; 
proc print data=fh noobs label; 
 var test z probchisq; 
 label test="Test" z="Z statistic" probchisq="P"; 
run;

%if &count>1 %then %do;

ods html options(pagebreak='no');
title h=10pt "Combination test";
proc print data=maxcombo noobs label; 
 format p pvalue10.4;
run;

%end;

ods proctitle;
ods html options(pagebreak='yes');
title;

*******************************************************;

* Clean up ;
proc datasets noprint;
 delete nlevels name one stat chisq q var covar
        fh maxcombo test zmax weights comb;
run; quit;

* Print runtime to log ;
data _null_;
  dur = datetime() - &_timer_start;
  put 24*'-' / ' Run time:' dur time13.2 / 24*'-';
run;

options source notes;

%mend;
