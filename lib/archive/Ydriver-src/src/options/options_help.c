/*
  License-Identifier: GPL
 
  Copyright (C) 2020 The Yambo Team
 
  Authors (see AUTHORS file for details): AM

*/
#include <stdio.h>
#include <kind.h>
#include <string.h>

void options_help(struct options_struct options[],int *i_opt)
{
 char *desc="Help & version";
 /*
  Help(s) 
 */
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="<string> can be an option (e.g. -h optics)";
 options[*i_opt].short_opt='h';
 options[*i_opt].long_opt="help";
 options[*i_opt].serial_var=1;
 options[*i_opt].optional_var=1;
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Code version & libraries";
 options[*i_opt].long_opt="version"; 
 options[*i_opt].serial_var=1;
 options[*i_opt].section=desc;

#if defined _example_driver
 desc="Self-Energy";
 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Hartree-Fock";
 options[*i_opt].long_opt="hf";
 options[*i_opt].short_opt='x';
 options[*i_opt].bin="example";
 options[*i_opt].yambo_string="HF_and_locXC";
 options[*i_opt].section=desc;

 *i_opt=*i_opt+1;
 options[*i_opt].short_desc="Self-Energy";
 strcpy(options[*i_opt].long_desc[0],"<string>=gw/(c)HOSEX (GW-kind self-energies)");
 strcpy(options[*i_opt].long_desc[1],"<string>=fan (Fan+DW self-energy)");
 strcpy(options[*i_opt].long_desc[2],"<string>=lamb (Lamb self-energy)");
 options[*i_opt].long_opt="se";
 options[*i_opt].short_opt='p';
 options[*i_opt].bin="example";
 options[*i_opt].yambo_string="self_energy";
 options[*i_opt].char_var=1;
 options[*i_opt].section=desc;
#endif
};
