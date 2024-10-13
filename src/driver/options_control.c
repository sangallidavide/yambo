/*
  License-Identifier: GPL
 
  Copyright (C) 2020 The Yambo Team
 
  Authors (see AUTHORS file for details): AM

*/
#include <string.h>
#include <stdio.h>
#include <kind.h>

void options_control(struct options_struct options[],int *i_opt)
{ 
 char *desc;
 /* 
  Input file 
 */
 desc="Input file & Directories";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Input file";
#if defined _a2y || defined _e2y
 options[*i_opt].short_desc="KSS or WFK file";
#endif
 options[*i_opt].short_opt='F';
 options[*i_opt].long_opt="Input"; 
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 options[*i_opt].bin="yambo ypp a2y e2y";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Input file variables verbosity";
 strcpy(options[*i_opt].long_desc[0],"<string> can be RL,kpt,sc,qp,io,gen,resp/X,ph,rt,par,nl,all");
 options[*i_opt].short_opt='V';
 options[*i_opt].long_opt="Verbosity"; 
 options[*i_opt].char_var=1;
 options[*i_opt].serial_var=1;
 options[*i_opt].yambo_string="infver";
 options[*i_opt].bin="yambo ypp";
 options[*i_opt].section=desc;
 /* 
  Utils
 */
 desc="Utilites";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Quiet input file creation";
 options[*i_opt].short_opt='Q';
 options[*i_opt].long_opt="Quiet"; 
 options[*i_opt].serial_var=1;
 options[*i_opt].bin="yambo ypp";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Verbose (fatter) log(s)";
 options[*i_opt].long_opt="fatlog";
 options[*i_opt].yambo_string="fatlog";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Assume experienced user";
 options[*i_opt].long_opt="expuser";
 options[*i_opt].yambo_string="expuser";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Databases properties";
 options[*i_opt].short_opt='D';
 options[*i_opt].long_opt="DBlist";
 options[*i_opt].serial_var=1;
 options[*i_opt].yambo_string="dbpr";
 options[*i_opt].bin="yambo";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Walltime";
 strcpy(options[*i_opt].long_desc[0],"Format is DdHhMm with D=days, H=hours and M=minutes");
 options[*i_opt].long_opt="walltime";
 options[*i_opt].int_var=1;
 options[*i_opt].yambo_string="wallt";
 options[*i_opt].bin="yambo";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Memory";
 strcpy(options[*i_opt].long_desc[0],"The value is per processor. The value can be provided in Mb/Gb. E.g. 1Gb");
 options[*i_opt].long_opt="memory";
 options[*i_opt].int_var=1;
 options[*i_opt].yambo_string="memorylimit";
 options[*i_opt].bin="yambo";
 options[*i_opt].section=desc;
#if defined _SCALAPACK
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="ScaLapacK test";
 options[*i_opt].long_opt="slktest";
 options[*i_opt].bin="yambo";
 options[*i_opt].yambo_string="slktest";
 options[*i_opt].section=desc;
#endif
 /* 
  Job control 
 */
 desc="Input file & Directories";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Job string";
 options[*i_opt].short_opt='J';
 options[*i_opt].long_opt="Job"; 
 options[*i_opt].char_var=1;
 options[*i_opt].yambo_string="jobstr";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Input directory";
 options[*i_opt].short_opt='I';
 options[*i_opt].long_opt="Idir"; 
 options[*i_opt].char_var=1;
 options[*i_opt].bin="yambo ypp p2y";
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="I/O directory";
 options[*i_opt].short_opt='O';
 options[*i_opt].long_opt="Odir";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Communication directory";
 options[*i_opt].short_opt='C';
 options[*i_opt].long_opt="Cdir";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 /* 
  Parallel
 */
 desc="Parallel Control";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Environment Parallel Variables file";
 options[*i_opt].short_opt='E';
 options[*i_opt].long_opt="parenv"; 
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
#if defined _MPI
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Switch off MPI support";
 options[*i_opt].long_opt="nompi";
 options[*i_opt].serial_var=1;
 options[*i_opt].section=desc;
#endif
#if defined _OPENMP
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Switch off OPENMP support";
 options[*i_opt].long_opt="noopenmp";
 options[*i_opt].yambo_string="noopenmp";
 options[*i_opt].section=desc;
#endif
};
