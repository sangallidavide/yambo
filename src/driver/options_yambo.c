/*
  License-Identifier: GPL
 
  Copyright (C) 2020 The Yambo Team
 
  Authors (see AUTHORS file for details): AM

*/
#include <string.h>
#include <stdio.h>
#include <kind.h>

void options_yambo(struct options_struct options[],int *i_opt)
{ 
 char *desc;
 int i_desc=0;
 desc="Self-Energy";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Hartree-Fock";
 options[*i_opt].long_opt="hf";
 options[*i_opt].short_opt='x';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="HF_and_locXC";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="GW approximation";
 strcpy(options[*i_opt].long_desc[i_desc],"<string>=(p)PA/(m)PA/(c)HOSEX/(r)eal-axis");
#if defined _ELPH
 i_desc=i_desc+1;
 strcpy(options[*i_opt].long_desc[i_desc],"<string>=fan");
#endif
#if defined _ELPH
 i_desc=i_desc+1;
 strcpy(options[*i_opt].long_desc[i_desc],"<string>=X");
#endif
 options[*i_opt].long_opt="gw0";
 options[*i_opt].short_opt='p';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="gw0";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Dyson Equation solver";
 strcpy(options[*i_opt].long_desc[0],"<string>=(g)reen [any scattering]");
 strcpy(options[*i_opt].long_desc[1],"<string>=(n)ewton [order 1]/(s)ecant [e-e scattering]");
#if defined _PHEL
 strcpy(options[*i_opt].long_desc[2],"<string>=(n)ewton [order 2] [p-e scattering]");
#endif
 options[*i_opt].long_opt="dyson";
 options[*i_opt].short_opt='g';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="dyson";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="GoWo Quasiparticle lifetimes";
 options[*i_opt].long_opt="lifetimes";
 options[*i_opt].short_opt='l';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="life";
 options[*i_opt].section=desc;

 desc="Initializations";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Initialization";
 options[*i_opt].short_opt='i';
 options[*i_opt].long_opt="setup";
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="setup";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Coulomb potential";
 options[*i_opt].long_opt="coulomb";
 options[*i_opt].short_opt='r';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="rim_cut";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Screened coulomb potential";
 options[*i_opt].long_opt="rw";
 options[*i_opt].short_opt='w';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="rim_w";
 options[*i_opt].section=desc;


 desc="Response Functions";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Linear Response optical properties";
 strcpy(options[*i_opt].long_desc[0],"<string>=c Reciprocal-Space");
 strcpy(options[*i_opt].long_desc[1],"<string>=b for Transition-Space Bethe-Salpeter");
 options[*i_opt].long_opt="optics";
 options[*i_opt].short_opt='o';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="optics";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Inverse Dielectric/Response Matrix";
 options[*i_opt].long_opt="X";
 options[*i_opt].short_opt='d';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="screen";
 options[*i_opt].section=desc;
 strcpy(options[*i_opt].long_desc[0],"<string>=(s)static/(p)PA/m(PA)/(d)ynamical dielectric matrix");
 strcpy(options[*i_opt].long_desc[1],"<string>=(X) dynamical response matrix");
 options[*i_opt].char_var=1;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Oscillator strenghts (or dipoles)";
 options[*i_opt].long_opt="dipoles";
 options[*i_opt].short_opt='q';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="dipoles";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Kernel";
 strcpy(options[*i_opt].long_desc[0],"<string>=hartree/alda/lrc/hf/sex/bsfxc");
 strcpy(options[*i_opt].long_desc[1],"hf/sex only eh-space; lrc only G-space");
 options[*i_opt].long_opt="kernel";
 options[*i_opt].short_opt='k';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="kernel";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;

 desc="Bethe-Salpeter Equation";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="BSE solver";
#if defined _SLEPC && !defined _NL
 strcpy(options[*i_opt].long_desc[0],"<string>=h/d/s/(p/f)i");
#else
 strcpy(options[*i_opt].long_desc[0],"<string>=h/d/(p/f)i");
#endif
 strcpy(options[*i_opt].long_desc[1],"(h)aydock/(d)iagonalization");
 strcpy(options[*i_opt].long_desc[2],"(pi) perturbative inversion/ (fi) full inversion");
#if defined _SLEPC && !defined _NL
 strcpy(options[*i_opt].long_desc[2],"(s)lepc partial diagonalization");
#endif
 options[*i_opt].long_opt="Ksolver";
 options[*i_opt].short_opt='y';
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="bss";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;

 desc="Total Energy";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="ACFDT Total Energy";
 options[*i_opt].long_opt="acfdt";
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="acfdt";
 options[*i_opt].section=desc;

};
