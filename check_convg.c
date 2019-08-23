   /********* Print out if converge or diverge *********/
   double **w = (app->w);
   double check = w[1][0];
   if(isinf(check)||isnan(check))
   {
    check=0;
   }
   else
   {
    check=1;
   }

   char    filename[255];
   FILE   *file;
   sprintf(filename, "%s.%d.%f", "out/advec-diff-imp.conv", ntime,nu);
   file = fopen(filename, "w");
   fprintf(file, "%f", check);
   fflush(file);
   fclose(file);
   /****************************************************/ 