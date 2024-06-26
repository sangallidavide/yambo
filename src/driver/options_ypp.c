/*
  License-Identifier: GPL
 
  Copyright (C) 2020 The Yambo Team
 
  Authors (see AUTHORS file for details): AM

*/
#include <string.h>
#include <stdio.h>
#include <kind.h>

void options_ypp(struct options_struct options[],int *i_opt)
{
 /* 
  Brillouin Zone 
 */
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="BZ Grid generator";
 strcpy(options[*i_opt].long_desc[0],"<string>=(k)pt,(q)pt,(s)hifted,(h)igh symmetry,(r)andom,r(e)gular");
 options[*i_opt].long_opt="grid";
 options[*i_opt].short_opt='k';
 options[*i_opt].bin="ypp";
 options[*i_opt].char_var=1;
 options[*i_opt].yambo_string="bzgrids";
 options[*i_opt].section="Brillouin Zone";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Fine to coarse grid Map";
 options[*i_opt].long_opt="map";
 options[*i_opt].short_opt='m';
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="kpts_map";
 options[*i_opt].section="Brillouin Zone";
 /* 
  Wannier
 */
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Wannier 90 interface";
 options[*i_opt].long_opt=  "wannier";
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="wannier";
 options[*i_opt].section="Wannier";
 /* 
  SOC
 */
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Perturbative SOC mapping";
 options[*i_opt].long_opt="soc";
 options[*i_opt].short_opt='w';
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="WFs_SOC_map";
 options[*i_opt].section="SOC";
 /* 
  Convertions
 */
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Remove symmetries not consistent with an external perturbation";
 options[*i_opt].long_opt=  "fixsym";
 options[*i_opt].short_opt='y';
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="fixsyms";
 options[*i_opt].section="Convertions";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Quasiparticle Databases";
 strcpy(options[*i_opt].long_desc[0],"<string>=(g)enerate-modify/(m)erge/(e)xpand");
 strcpy(options[*i_opt].long_desc[1],"         (e)xpand uses the symmetries to generate a BZ-expanded QP database");
 options[*i_opt].long_opt="qpdb";
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="QPDBs";
 options[*i_opt].section="Convertions";
 options[*i_opt].char_var=1;
#if !defined _YPP_RT  
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="gkkp databases";
 options[*i_opt].long_opt= "gkkp";
 options[*i_opt].short_opt='g';
 options[*i_opt].bin="ypp_ph";
 options[*i_opt].yambo_string="gkkp"  ;
 options[*i_opt].section="Convertions";
 strcpy(options[*i_opt].long_desc[0],"<string>=(g)kkp,(d)ouble grid,(p)lot gkkp");
 options[*i_opt].char_var=1;
#endif
 /* 
  Plots
 */
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Excitonic State Momentum";
 options[*i_opt].long_opt= "BSiq";
 options[*i_opt].short_opt='b';
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="BSiq";
 options[*i_opt].section="Plots";
 options[*i_opt].int_var=1;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Mean Potential";
 options[*i_opt].long_opt= "potential";
 options[*i_opt].short_opt='v';
 options[*i_opt].bin="ypp_sc";
 options[*i_opt].yambo_string="MeanPot";
 options[*i_opt].section="Plots";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Electronic properties";
 strcpy(options[*i_opt].long_desc[0],"<string>=(h)artree,(f)ock,(coh),(sex),(cohsex),(exx),(exxc),(srpa),(d)ef,(ip)");
#if defined _ELPH  
 strcpy(options[*i_opt].long_desc[0],"<string>=(w)ave,(d)ensity,(m)ag,do(s),(b)ands,(c)urrent,(e)lias");
#elif defined _YPP_MAGNETIC 
 strcpy(options[*i_opt].long_desc[0],"<string>=(w)ave,(d)ensity,(m)ag,do(s),(b)ands,(c)urrent,angu(l)ar,(p)osition");
#else
 strcpy(options[*i_opt].long_desc[0],"<string>=(w)ave,(d)ensity,(m)ag,do(s),(b)ands,(c)urrent");
#endif
 options[*i_opt].long_opt="electron";
 options[*i_opt].short_opt='s';
 options[*i_opt].bin="ypp"; 
 options[*i_opt].yambo_string="electrons";
 options[*i_opt].section="Plots";
 options[*i_opt].char_var=1;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Excitonic properties";
 options[*i_opt].long_opt="exciton";
 options[*i_opt].short_opt='e';
 options[*i_opt].bin="ypp ypp_ph";
#if defined _ELPH  
 strcpy(options[*i_opt].long_desc[0],"<string>=(s)ort,(sp)in,(a)mplitude,(w)ave,(i)nterpolate,");
 strcpy(options[*i_opt].long_desc[1],"         (e)lias,(g)kkp,(p)h-assisted dos");
#else
 strcpy(options[*i_opt].long_desc[0],"<string>=(s)ort,(sp)in,(a)mplitude,(w)ave,(i)nterpolate");
#endif
 options[*i_opt].yambo_string="excitons";
 options[*i_opt].section="Plots";
 options[*i_opt].char_var=1;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Dipole properties";
 options[*i_opt].long_opt="dipoles";
 options[*i_opt].bin="ypp";
 strcpy(options[*i_opt].long_desc[0],"<string>=(exc)itonic,(ip)independent-particle");
#if defined _YPP_RT  
 strcpy(options[*i_opt].long_desc[0],"<string>=(exc)itonic,(ip)independent-particle,(m)ask");
#endif
 options[*i_opt].yambo_string="dipoles";
 options[*i_opt].section="Plots";
 options[*i_opt].char_var=1;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Free hole position in the excitonic plot";
 options[*i_opt].long_opt="freehole";
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="freehole";
 options[*i_opt].section="Plots";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Average hole/electron wavefunction";
 options[*i_opt].long_opt="avehole";
 options[*i_opt].bin="ypp";
 options[*i_opt].yambo_string="avehole";
 options[*i_opt].section="Plots";
#if !defined _YPP_RT  
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Phononic properties";
 strcpy(options[*i_opt].long_desc[0],"<string>=(d)os,(e)lias,(a)mplitude");
 options[*i_opt].long_opt=  "phonon";
 options[*i_opt].short_opt='p';
 options[*i_opt].bin="ypp_ph";
 options[*i_opt].yambo_string="phonons";
 options[*i_opt].section="Plots";
 options[*i_opt].char_var=1;
#endif
/*
  Real-Time
*/
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Non-linear response analysis";
 options[*i_opt].long_opt="nl";
 options[*i_opt].short_opt='u';
 options[*i_opt].bin="ypp_nl";
 options[*i_opt].yambo_string="nonlinear";
 options[*i_opt].section="Real-Time";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Carriers database generation";
 strcpy(options[*i_opt].long_desc[0],"<string>=(e)nergy,(f)ermi");
 options[*i_opt].long_opt=  "rtdb";
 options[*i_opt].short_opt='c';
 options[*i_opt].bin="ypp_rt";
 options[*i_opt].yambo_string="RTDBs";
 options[*i_opt].section="Real-Time";
 options[*i_opt].char_var=1;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="TD observables plot";
 options[*i_opt].long_opt="rtplot";
 options[*i_opt].short_opt='n';
 options[*i_opt].char_var=1;
 options[*i_opt].bin="ypp_rt";
 options[*i_opt].yambo_string="TDplots"; /* TDplots */
 options[*i_opt].section="Real-Time";
 strcpy(options[*i_opt].long_desc[0],"<string>=(X)response,(a)bsorption,(o)ccupations,(l)ifetimes,(d)ensity,(p)olariazion,(g)reen-function");
 strcpy(options[*i_opt].long_desc[1]," ");
 strcpy(options[*i_opt].long_desc[2],"(X) response calculates the response via the time-resolved polarization");
 strcpy(options[*i_opt].long_desc[3],"(a) absorption amends the Kubo expression with the time-dependent occupations");
 strcpy(options[*i_opt].long_desc[4],"(p) polarization evaluates the k-resolved components of the time-dependent polarization");
 strcpy(options[*i_opt].long_desc[5],"(g) evaluates the two-times Green`s function");
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="TD plot control";
 options[*i_opt].char_var=1;
 options[*i_opt].long_opt= "rtmode";
 options[*i_opt].short_opt='t';
 options[*i_opt].bin="ypp_rt";
 options[*i_opt].yambo_string="TDplotmode"; /* TDpol */
 options[*i_opt].section="Real-Time";
 strcpy(options[*i_opt].long_desc[0],"rtplot=X/a => <string>=(t)ime");
 strcpy(options[*i_opt].long_desc[1],"rtplot=o   => <string>=(b)ands,(t)ime,(e)nergy,(d)os");
 strcpy(options[*i_opt].long_desc[2],"rtplot=l   => <string>=(b)ands,(t)ime,(e)nergy");
 strcpy(options[*i_opt].long_desc[3],"rtplot=d   => <string>=(t)ime");
 strcpy(options[*i_opt].long_desc[4],"rtplot=p   => <string>=(t)ime");
};
