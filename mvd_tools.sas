/*******************************************************************************
| Name       : mvd_tools.sas
| Purpose    : Creates FCMP Functions to calculate PDFs, CDFs
|              and simulate values from n-variate Gauss and T distributions
|              these can be used for joint simulation and likelihood 
|              creation.
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 28NOV18 
|--------------------------------------------------------------------------------
| FCMP Functions and Call Routine List:
|--------------------------------------------------------------------------------
|
| Multivariate Normal Functions and Call Routines
|------------------------------------------------
| Name     : rand_mvn(y,m,v) 
| Purpose  : Simulates multivariate normal values  
| Arguments: y [REQ] = 1 Dim array with mv variables to populate
|            m [REQ] = 1 Dim array with mv mean values
|            v [REQ] = 2 Dim array with variance covariance matrix  
|
| Name     : pdf_mvn(y,m,v) 
| Purpose  : Returns multivariate normal density values  
| Arguments: y [REQ] = 1 Dim array with mv data values
|            m [REQ] = 1 Dim array with mv mean values
|            v [REQ] = 2 Dim array with variance covariance matrix  
|
|
| Multivariate T Functions and Call Routines
|-------------------------------------------
| Name     : rand_mvt(y,m,v,d) 
| Purpose  : Simulates multivariate normal values  
| Arguments: y [REQ] = 1 Dim array with mv variables to populate
|            m [REQ] = 1 Dim array with mv mean values
|            v [REQ] = 2 Dim array with variance covariance matrix  
|            d [REQ] = Degrees of freedom for MVT distribution
|
|
***********************************************************************************/;

proc fcmp outlib = work.functions.mvd_tools;

  ********************************************************************************;
  *** MVN DISTRIBUTION RELATED FUNCTIONS AND CALL ROUTINES                     ***;
  ********************************************************************************;

  subroutine rand_mvn(y[*],m[*],v[*,*]);   
     outargs y;
     ydim1 = dim1(y); 
     mdim1 = dim1(m);
     vdim1 = dim1(v);
     vdim2 = dim2(v);
     if not ( ydim1 = mdim1 = vdim1 = vdim2 ) then do;
       msg1 = "ER"||upcase("ror:(FCMP):")||"The Function SIM_MVN does not have matching array sizes.";
       msg2 = "ER"||upcase("ror:(FCMP): Dimensions:");
       put msg1;
       put msg2 ydim1= mdim1= vdim1= vdim2=;
     end;
     else do;
       n = ydim1;
       array zvec [1,1] / nosymbols;     
       array mvec [1,1] / nosymbols;
       array vmat [1,1] / nosymbols;
       call dynamic_array(zvec,n,1);     
       call dynamic_array(mvec,n,1);
       call dynamic_array(vmat,n,n);
       do ii = 1 to n;
         zvec[ii,1] = rand("NORMAL",0,1);
         mvec[ii,1] = m[ii];
         do jj = 1 to n;
           vmat[ii,jj] = v[ii,jj];
         end;
       end;
       array smat[1,1] / nosymbols;
       array sxz [1,1] / nosymbols;
       array mvn [1,1] / nosymbols;
       call dynamic_array(smat,n,n);     
       call dynamic_array(sxz,n,1);
       call dynamic_array(mvn,n,1);
       call chol(vmat,smat);
       call mult(smat,zvec,sxz);
       call addmatrix(mvec,sxz,mvn);
       do ii = 1 to n;
         y[ii] = mvn[ii,1];
       end;
     end;
  endsub;
 

  function pdf_mvn(y[*],m[*],v[*,*]);
    ydim1 = dim1(y); 
    mdim1 = dim1(m); 
    vdim1 = dim1(v);
    vdim2 = dim2(v);
    if not ( ydim1 = mdim1 = vdim1 = vdim2 ) then do;
      msg1 = "ER"||upcase("ror:(FCMP):")||"The Function MVN_PDF does not have matching array sizes.";
      msg2 = "ER"||upcase("ror:(FCMP): Dimensions:");
      put msg1;
      put msg2 ydim1= mdim1= vdim1= vdim2=;
      pdf = .e;
    end;
    else do;

      n = 0;
      do ii = 1 to ydim1;
        if y[ii] ne . then n = n + 1;
      end;

      if (n < ydim1) then do;
        msg1 = "NO"||upcase("te:(FCMP):")||"Partially missing vectors supplied to function MVN_PDF.";
        msg2 = "NO"||upcase("te:(FCMP):")||"PDF_MVN will return pdf values using vectors constructed with non missing values.";
        put msg1;
        put msg2;
      end;

      array loc[1,1] / nosymbols;
      call dynamic_array(loc,n,1);
      jj = 0;
      do ii = 1 to ydim1;
        if y[ii] ne . then do;
          jj = jj + 1;
          loc[jj,1] = ii;
	    end;
	  end;

      array dvec[1,1] / nosymbols;     
      array vmat[1,1] / nosymbols;
      array dtrs[1,1] / nosymbols;
      array vinv[1,1] / nosymbols;
      array vec1[1,1] / nosymbols;
      array val1[1,1] / nosymbols;
      call dynamic_array(dvec,n,1);
      call dynamic_array(dtrs,1,n);
      call dynamic_array(vmat,n,n);
      call dynamic_array(vinv,n,n);   
      call dynamic_array(vec1,n,1);  

      do ii = 1 to n;
        dvec[ii,1] = y[loc[ii,1]] - m[loc[ii,1]];
        dtrs[1,ii] = dvec[ii,1];
        do jj = 1 to n;
          vmat[ii,jj] = v[loc[ii,1],loc[jj,1]];
        end;
      end;

      call det(vmat,detv);
      call inv(vmat,vinv);
      call mult(vinv,dvec,vec1);
      call mult(dtrs,vec1,val1);

      pi  = arcos(-1);
      pdf = ((2*pi)**(-n/2)) * ((detv)**(-1/2)) * exp(-0.5*val1[1,1]);

    end;
    return(pdf);
  endsub; 


  ********************************************************************************;
  *** MVT DISTRIBUTION RELATED FUNCTIONS AND CALL ROUTINES                     ***;
  ********************************************************************************;

  subroutine rand_mvt(y[*],m[*],v[*,*],d);   
     outargs y;

     ydim1 = dim1(y); 
     mdim1 = dim1(m);
     vdim1 = dim1(v);
     vdim2 = dim2(v);
     if not ( ydim1 = mdim1 = vdim1 = vdim2 ) then do;
       msg1 = "ER"||upcase("ror:(FCMP):")||"The Function SIM_MVT does not have matching array sizes.";
       msg2 = "ER"||upcase("ror:(FCMP): Dimensions:");
       put msg1;
       put msg2 ydim1= mdim1= vdim1= vdim2=;
     end;
     else do;

       n = ydim1;
       array zvec [1,1] / nosymbols;     
       array mvec [1,1] / nosymbols;
       array vmat [1,1] / nosymbols;
       array wval [1,1] / nosymbols;
       call dynamic_array(zvec,n,1);     
       call dynamic_array(mvec,n,1);
       call dynamic_array(vmat,n,n);

       do ii = 1 to n;
         zvec[ii,1] = rand("NORMAL",0,1);
         mvec[ii,1] = m[ii];
         do jj = 1 to n;
           vmat[ii,jj] = v[ii,jj];
         end;
       end;

       wval[1,1] = sqrt(d/rand("CHISQ",d)); 

       array smat[1,1] / nosymbols;
       array sxz [1,1] / nosymbols;
       array szw [1,1] / nosymbols;
       array mvt [1,1] / nosymbols;
       call dynamic_array(smat,n,n);     
       call dynamic_array(sxz,n,1);
       call dynamic_array(szw,n,1);
       call dynamic_array(mvt,n,1);

       call det(vmat,detv);
       if detv le 0 then do;
         msg1 = "ER"||upcase("ror:(FCMP):")||"The input array v[.,.] is not positive definite.";
         msg2 = "ER"||upcase("ror:(FCMP): Array:");
         put msg1;
         put msg2 v[*]=;
       end;
       else do;
         call chol(vmat,smat);
         call mult(smat,zvec,sxz);
         call mult(sxz,wval,szw);
         call addmatrix(mvec,szw,mvt);
         do ii = 1 to n;
           y[ii] = mvt[ii,1];
         end;
       end;

     end;
  endsub;

run;
quit;

options cmplib = work.functions; 
  



