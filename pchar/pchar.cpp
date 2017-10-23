/****************************************************************************
This program determines the burn rate of a propellant given a pressure 
versus time input file. 

Adapted from Ballisti.c by Jeroen Louwers which is based on 
"Determination of the Interior Ballistic Properties of Solid Propellants 
from Progressive Burning Grains Data." FROM: J. Spacecraft Vol. 5, No.5, 
623-624 (1968)

   Usage: balli [-f infile]

   Author: Mark Spiegl

***************************************************************************/

// Remove this line for non-Windows compile
#define WINDOWS_COMPILE

// To keep Visual Studio quiet about deprecated functions
#define _CRT_SECURE_NO_WARNINGS   

#ifdef WINDOWS_COMPILE
#include "stdafx.h"
#include <stdlib.h> 
#include <ctype.h> 
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAXDATA 2000
#define MAXSTR  80

void iterate();
void input_raw_data(char *);
void get_output_file(char *);
double a_mean(double *, int); 
double slope_bf_line(double *, double *, int, int); 
double y_intercept (double *, double *, int, int); 


double pressure[MAXDATA];
double core_radius;
double grain_radius;
double grain_length;
double begin_time;
double end_time;
double delta_t;
FILE *outfp;


int
main(int argc, char argv[])
{
   int cc;
   char *infile = (char *) 0;
   extern char *optarg;
   int err = 0;
   extern int getopt (int,char**, char*);

   /*
    * Parse input parameters
    */
   while ((cc = getopt(argc, (char **)argv, "f:")) != EOF) {
      switch (cc) {
      case 'f':
         infile = optarg;
         break;
      default:
         err++;
         break;
      }
   }

   if (err) {
      (void) fprintf(stderr, "Usage: %s [-f infile] \n", argv[0]);
      return(-1);
   }
   /*
    * Get the data for simulation run
    */
   input_raw_data(infile);
   get_output_file(infile);

   /*
    * Execute the simulation
    */
   iterate();

   return (0);
}


/**********************************************************************/
/* Actual Simulation                                                  */
/**********************************************************************/
void
iterate()
{
   double cur_radius, cur_time, b, rate1, pres1, pres2;
   static double rate[MAXDATA];
   static double log_rate[MAXDATA]; 
   static double log_pressure[MAXDATA]; 
   double tol;
   int i;

   tol = 0.005;

   (void) printf("Iterating for rate (t=%lf,%lf) ....\n", begin_time, end_time);
   rate[(int) (1 / delta_t * begin_time + 0.5)] = 0.00;
   for (;;) {
      cur_time = begin_time;
      cur_radius = core_radius;
      while (cur_time <= end_time + 0.5 * delta_t) {
         b = delta_t / (2 * cur_radius);
         rate1 = rate[(int) (1 / delta_t * cur_time + 0.5)];
         pres1 = pressure[(int) (1 / delta_t * cur_time + 0.5)];
         pres2 = pressure[(int) (1 / delta_t * cur_time + 1.5)];
         rate[(int) (1 / delta_t * cur_time + 1.5)] =
            (-(1 + b * rate1) + sqrt((1 + b * rate1) * (1 + b * rate1) +
                                  4 * b * rate1 * pres2 / pres1)) / (2 * b);
         cur_radius = cur_radius + delta_t *
            (rate1 + rate[(int) (1 / delta_t * cur_time + 1.5)]) / 2;
         cur_time += delta_t;
      }
      if (cur_radius > grain_radius + tol / 2.0) {
         rate[(int) (1 / delta_t * begin_time + 0.5)] =
            rate[(int) (1 / delta_t * begin_time + 0.5)] -
            (cur_radius - grain_radius) / (end_time - begin_time);
      } else {
         if (cur_radius < grain_radius - tol / 2.0) {
            rate[(int) (1 / delta_t * begin_time + 0.5)] =
               rate[(int) (1 / delta_t * begin_time + 0.5)] -
               (cur_radius - grain_radius) / (end_time - begin_time);
         } else {
            break;
         }
      }
   }
   (void) printf("Done\n");

   // Run least square regression of Pressure v Log(burnRate) to get "a" and "n"
   cur_time = begin_time;
   while (cur_time < end_time + 0.5 * delta_t) {
      i = (int) (1 / delta_t * cur_time + 0.5);
      log_rate[i] = log (rate[i]); 
      log_pressure[i] = log (pressure[i]);
      cur_time += delta_t;
   }
   double slope = slope_bf_line(log_pressure, log_rate, i, i); 
   double inctp = y_intercept(log_pressure, log_rate, i, i); 
   fprintf (outfp, "# Rate=a*Pres^n\n"); 
   fprintf (outfp, "# a=%6.5f\n", exp(inctp)); 
   fprintf (outfp, "# n=%6.5f\n", slope); 


   // Now output a matrix of Pressure, burnRate
   fprintf (outfp, "#\n#Column1=pressure...Column2=burnRate\n#\n"); 
   cur_time = begin_time;
   while (cur_time < end_time + 0.5 * delta_t) {
      i = (int) (1 / delta_t * cur_time + 0.5);
      (void) fprintf(outfp, "%6.5f %6.5f\n", pressure[i], rate[i]);
      cur_time += delta_t;
   }
   (void) fclose(outfp);
}


/**********************************************************************/
/* read the input parmaters                                           */
/**********************************************************************/
void
input_raw_data(char *fname)
{
   char str[MAXSTR];
   void get_from_file(char *, char *, char *);
   FILE *fp;
   char *tmp;
   double t;
   int no_datapoints = 0;
   int line_no = 0; 

   if (fname == (char *) 0) {
      (void) printf("Outside Grain Diam. (inches): ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "GRAIN_DIAMETER", str);
   (void) sscanf(str, "%lf", &grain_radius);
   grain_radius /= 2.0;

   if (fname == (char *) 0) {
      (void) printf("Grain Core Diameter (inches): ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "CORE_DIAMETER", str);
   (void) sscanf(str, "%lf", &core_radius);
   core_radius /= 2.0;

   if (fname == (char *) 0) {
      (void) printf("Grain Length (inches): ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "GRAIN_LENGTH", str);
   (void) sscanf(str, "%lf", &grain_length);

   if (fname == (char *) 0) {
      (void) printf("Begin Time: ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "BEGIN_TIME", str);
   (void) sscanf(str, "%lf", &begin_time);

   if (fname == (char *) 0) {
      (void) printf("End Time: ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "END_TIME", str);
   (void) sscanf(str, "%lf", &end_time);

   if (fname == (char *) 0) {
      (void) printf("Delta t (sec): ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "DELTA_T", str);
   (void) sscanf(str, "%lf", &delta_t);

   if (fname == (char *) 0) {
      (void) printf("Time vs Pressure datafile: ");
      (void) fgets(str, MAXSTR, stdin);
   } else
      get_from_file(fname, "DATAFILE", str);
   for (tmp = str; (*tmp != '\n') && (*tmp != '\0'); tmp++);
   *tmp = (char) 0;
   if ((fp = fopen(str, "r")) == NULL) {
      (void) fprintf(stderr, "Cannot open datafile %s\n", str);
      exit(-1);
   }
   for (;;) {
      line_no++; 
      if (fgets(str, MAXSTR, fp) == (char *) 0)
         break; 
      if (*str == '#')
         continue; 
      if (sscanf(str, "%lf %lf", &t, &pressure[no_datapoints++]) != 2)
         break;
   }
   if (!feof(fp)) {
      (void) fprintf(stderr, "Invalid entry in datafile at line %d\n",
                     line_no);
      exit(-1);
   }
   (void) fclose(fp);
}

/**********************************************************************/
/* Get the paramaters from an input file instead of stdin             */
/**********************************************************************/
void
get_from_file(char *fname, char *keyword, char *ret)
{
   char str[MAXSTR];
   FILE *fp;
   char *strcpy(char *, const char *);
   char *t1, *t2;
   int len = (int) strlen(keyword);

   if ((fp = fopen(fname, "r")) == (FILE *) NULL) {
      (void) fprintf(stderr, "Cannot read input file %s\n", fname);
      exit(-1);
   }
   for (;;) {
      if (fgets(str, MAXSTR, fp) == (char *) 0) {
         (void) fprintf(stderr, "Could not find keyword %s\n", keyword);
         exit(-1);
      }
      if (strncmp(str, keyword, len) == 0) {
         for (t1 = str; !isspace(*t1); t1++);
         for (; isspace(*t1); t1++);
         for (t2 = t1; !isspace(*t2); t2++);
         *t2 = (char) 0;
         (void) strcpy(ret, t1);
         (void) fclose(fp);
         return;
      }
   }
}


/**********************************************************************/
/* Setup the output file                                              */
/**********************************************************************/
void
get_output_file(char *fname)
{
   char outfile[MAXSTR];
   void get_from_file(char *, char *, char *);

   if (fname == (char *) 0) {
      (void) printf("Output file [return for stdout]: ");
      (void) gets(outfile);
   } else
      get_from_file(fname, "OUTFILE", outfile);
   if ((strncmp(outfile, "stdout", strlen("stdout")) == 0) || (*outfile == (char) 0))
      outfp = stdout;
   else if ((outfp = fopen(outfile, "w")) == (FILE *) NULL) {
      (void) fprintf(stderr, "Cannot open output file %s\n", outfile);
      exit(-1);
   }
}


/**********************************************************************/
/* Return the slope of the line of                                    */
/* best fit of a set of x,y data. (output OK)                         */
/* Copied from Sgt Pepper's Statistics Library
/**********************************************************************/

double 
slope_bf_line(double *xlist, double *ylist, int xn, int yn)
{
    int idx;
    double num, den, xmean, ymean;

    if(xn != yn)
    {
        fprintf(stderr,"error: data sets unequal\n"
                       " in slope_bf_line()\n");
        exit(-1);
    }

    xmean = a_mean(xlist, xn);
    ymean = a_mean(ylist, yn);
    num = 0.0;
    den = 0.0;
    for(idx = 0; idx < xn; idx++)
    {
        num += (xlist[idx] - xmean) * (ylist[idx] - ymean);
        den += (xlist[idx] - xmean) * (xlist[idx] - xmean);
    }
    
    return num / den;       
}

/**********************************************************************/
/* Return the y-intercept of the line                                 */
/* of best fit of a set of x,y data. (output OK)                      */
/* Copied from Sgt Pepper's Statistics Library
/**********************************************************************/
double 
y_intercept(double *xlist, double *ylist, int xn, int yn)
{
    double xmean, ymean, slope;

    if(xn != yn)
    {
        fprintf(stderr,"error: data sets unequal\n"
                       " in y_intercept()\n");
        exit(-1);
    }

    xmean = a_mean(xlist, xn);
    ymean = a_mean(ylist, yn);
    slope = slope_bf_line(xlist, ylist, xn, yn);
    
    return ymean - (slope * xmean);
}



/**********************************************************************/
/* Return the arithmetic mean                                         */
/* of values in datalist. (output OK)                                 */  
/* Copied from Sgt Pepper's Statistics Library
/**********************************************************************/
double 
a_mean(double *datalist, int listsize)
{
    int idx;
    double sum;
  
    sum = 0;
    for(idx = 0; idx < listsize; ++idx)
    {
        sum += datalist[idx];
    }

    return sum / listsize;
}

