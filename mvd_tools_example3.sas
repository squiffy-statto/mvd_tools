/*******************************************************************************
| Name       : mvd_tools_example3.sas
| Purpose    : Examples of using FCMP Functions from mvd_tools.
| SAS Version: 9.4
| Created By : Thomas Drury
| Date       : 20SEP19 
********************************************************************************/;

*** INCLUDE TOOLS CODE ***;
options source2;
filename mvd url "https://raw.githubusercontent.com/squiffy-statto/mvd_tools/master/mvd_tools.sas";
%include mvd;
options cmplib=work.functions;


********************************************************************************;
*** SHOW A 3-DIM VECTOR WITH MISSINGS GIVES THE SAME AS A 2-DIM VECTOR       ***;
********************************************************************************;

data tvn1;

  array x[3] x31 x32 x33;
  array m[3]   _temporary_ (0 0 0);
  array v[3,3] _temporary_ (1.0 0.8 0.5 
                            0.8 1.0 0.8
                            0.5 0.8 1.0);
  do x31 = -2 to 2 by 1;
  do x32 = -2 to 2 by 1;
  do x33 = -2 to 2 by 1;
    i = x31;
	j = x32;
    if rand("uniform") < 0.5 then x33 = .;
    pdf3 = pdf_mvn(x,m,v);
    output;
  end;
  end;
  end; 

run;

data bvn1;

  array x[2] x21 x22;
  array m[2]   _temporary_ (0 0);
  array v[2,2] _temporary_ (1.0 0.8  
                            0.8 1.0);
  do x21 = -2 to 2 by 1;
  do x22 = -2 to 2 by 1;
    i = x21;
	j = x22;
    pdf2 = pdf_mvn(x,m,v);
    output;
  end;
  end;

run;

data tvn2;
  merge tvn1 
        bvn1;
  by i j;
run;

proc print data = tvn2;
run;



********************************************************************************;
*** CREATE CONTOUR OF 2-DIM PDF                                              ***;
********************************************************************************;

data bvn2;

  array x[2] x21 x22;
  array m[2]   _temporary_ (0 0);
  array v[2,2] _temporary_ (1.0 0.8  
                            0.8 1.0);
  do x21 = -3 to 3 by 0.05;
  do x22 = -3 to 3 by 0.05;
    i = x21;
	j = x22;
    pdf = pdf_mvn(x,m,v);
    output;
  end;
  end;

run;

*** SIMPLE GTL TEMPLATE FOR CONTOUR PLOT ***;
proc template;
  define statgraph ContourPlotParm;
  dynamic _X _Y _Z _TITLE;
  begingraph;
     entrytitle _TITLE;
     layout overlay;
        contourplotparm x=_X y=_Y z=_Z / contourtype=fill nhint=12  colormodel=twocolorramp name="Contour";
        continuouslegend "Contour" / title=_Z;
     endlayout;
  endgraph;
  end;
run;


proc sgrender data=bvn2 template=ContourPlotParm;
  dynamic _X="x21" _Y="x22" _Z="pdf" _TITLE="Graph of bivariate normal distribution";
run;




