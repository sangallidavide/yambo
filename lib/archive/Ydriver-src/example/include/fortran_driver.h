/*
  License-Identifier: GPL
 
  Copyright (C) 2020 The Yambo Team
 
  Authors (see AUTHORS file for details): DS
*/
/* 
 Tool drivers 
*/
/* MAIN 
*/
#if defined _example_driver
 #if defined _FORTRAN_US
  int main_
 #else
  int main
 #endif
#endif
 (int *, int *,int *,int *,int *,int *,int *,int *,
  char *string, char *in_file, char *in_dir, char *out_dir, char *com_dir, char *job,
  int string_N, int in_file_N, int in_dir_N, int out_dir_N, int com_dir_N, int job_N);
