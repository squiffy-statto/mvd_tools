/*******************************************************************************
| Name       : mvd_tools.sas
| Purpose    : Creates macros and FCMP Functions to calculate PDFs, CDFs
|              and simulate values from n-variate gauss and T distributions
|              these can be used for joint simulation and likelihood 
|              creation.
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 28NOV18 
|--------------------------------------------------------------------------------
| Macros List: See corresponding FCMP Function
|--------------------------------------------------------------------------------
|--------------------------------------------------------------------------------
| FCMP Functions and Call Routine List:
|--------------------------------------------------------------------------------
| Name     : mvn_sim(y,m,v) 
| Purpose  : Simulates multivariate normal values  
| Arguments: y [REQ] = 1 Dim array with mv variables to populate
|            m [REQ] = 1 Dim array with mv mean values
|            v [REQ] = 2 Dim array with variance covariance matrix  
***********************************************************************************/;

proc fcmp outlib = work.functions.sim_tools;

  subroutine sim_mvn(y[*],m[*],v[*,*]);   
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
run;
quit;

options cmplib = work.functions; 
  



