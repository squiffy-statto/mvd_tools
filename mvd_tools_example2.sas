/*******************************************************************************
| Name       : mvd_tools_example2.sas
| Purpose    : Examples of using FCMP Functions from mvd_tools.
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 28NOV18 
********************************************************************************/;

*** INCLUDE TOOLS CODE ***;
options source2;
filename mvd url "https://raw.githubusercontent.com/squiffy-statto/mvd_tools/master/mvd_tools.sas";
%include mvd;

*** NON POS DEF ERROR ***;
data mvt_data1;
   call streaminit(123456);
   array y[2] y1-y2;
   array m[2]   _temporary_ (1 2);
   array r[2,2] _temporary_ (1.00 2.00
                             2.00 1.00);
   df = 3;

   do sim = 1 to 10;
     call sim_mvt(y,m,r,df);
     output;
   end;
run;


*** STANDARD EXAMPLE ***;
data mvt_data2;
   call streaminit(123456);
   array y[8] y1-y8;
   array m[8]   _temporary_ (1 2 3 4 5 6 7 8);
   array r[8,8] _temporary_ (1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30
                             0.90 1.00 0.90 0.80 0.70 0.60 0.50 0.40 
                             0.80 0.90 1.00 0.90 0.80 0.70 0.60 0.50  
                             0.70 0.80 0.90 1.00 0.90 0.80 0.70 0.60 
                             0.60 0.70 0.80 0.90 1.00 0.90 0.80 0.70  
                             0.50 0.60 0.70 0.80 0.90 1.00 0.90 0.80 
                             0.40 0.50 0.60 0.70 0.80 0.90 1.00 0.90  
                             0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00);
   df = 3;

   do sim = 1 to 10000;
     call sim_mvt(y,m,r,df);
     output;
   end;
run;


*** CHECK DEPENDENCE ***;
ods html;
proc corr data = mvt_data2;
  var y1-y8;
run;
ods html close;


proc datasets lib = work noprint;
  delete mvt_data:;
run;
quit;






