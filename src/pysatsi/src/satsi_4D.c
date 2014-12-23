 //-------------------------------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//-------------------------------------------------------------------------------------------------

/* COORDINATES ARE EAST,NORTH,UP */

//-------------------------------------------------------------------------------------------------
// SATSI_4D
//
// Original code:
//   Jeanne Hardebeck <jhardebeck@usgs.gov>
//   available at: http://earthquake.usgs.gov/research/software/
// 
// Corrections to the original code:
//   Grzegorz Kwiatek [GK] <kwiatek@gfz-potsdam.de> <http://www.sejsmologia-gornicza.pl/about>
//   Patricia Martinez-Garzon [PM] <patricia@gfz-potsdam.de>
// 
//   Code updated to C99 standard. 
//
// $Last revision: 1.0 $  $Date: 2012/07/11  $  
//-------------------------------------------------------------------------------------------------

#define TODEG 57.29577951
//-------------------------------------------------------------------------------------------------
// [GK 2013.05.15] Original satsi setup.
//#define MAXDATA 7200
//#define MAXX 14
//#define MAXY 14
//#define MAXZ 14
//#define MAXT 14
//#define MAXBOX 1764
//#define SRMAXBOX 42
// [GK 2013.05.15] Standard version.
#define MAXDATA 70000
#define MAXX 50
#define MAXY 50
#define MAXZ 50
#define MAXT 50
#define MAXBOX 10000
#define SRMAXBOX 100
// [GK 2013.05.15] Extended SATSI setup.
//#define MAXDATA 70000
//#define MAXX 200
//#define MAXY 200
//#define MAXZ 200
//#define MAXT 200
//#define MAXBOX 10000
//#define SRMAXBOX 100
//-------------------------------------------------------------------------------------------------

// [GK 2013.03.03] Additional declarations to suppress warning messages;
void sprsax(double sa[], int ija[], double x[], double b[], int m, int n);
void leasq_sparse(int a_ija[], double a_sa[], int d_ija[], double d_sa[], int m,
    int n, int p, double x[], double b[]);

//-------------------------------------------------------------------------------------------------
int main(argc, argv)
  /* slickenside inversion program */
  int argc; /* argument count */
  char **argv; /* argument string */
  {
  double *ddir, *dip, *rake; /* focal mechanism data */
  int *x, *y, *dep, *time; /* focal mechanism bins */
  int nobs, nloc, nrows, nlocfill; /* number of observations, bins, rows */
  double cwt, twt; /* damping parameter, time/space damping ratio */
  double *diag_sa, *amat_sa, *d_sa; /* inversion matrices in sparse matrix, */
  int *diag_ija, *amat_ija, *d_ija; /*      row-indexed form */
  int *lindex, nx, ny, nz, nt, index; /* book-keeping, which bins where in matrix */
  double stress[5 * MAXBOX]; /* stress tensor in vector form xx,xy,xz,yy,yz,zz */
  double *slick; /* slickenside vector elements vector */
  double n1, n2, n3; /* normal vector elements */
  //double norm[MAXDATA][3]; /* storage of n1,n2,n3 */
  double angavg, angstd; /* average and standard deviation of fit angle */
  double magavg, magstd; /* same for tangential stress size */
  double *slick_pre; /* predicted slip vector */
  char line[80]; /* character line */
  FILE *fpin; /* input file pointer */
  FILE *fpout; /* output file pointer */
  int i, j, k, k2, ip, m, n, p, i2; /* dummy variables */
  double z, z2, z3, temp[5], ls1, ls2, ls3; /* more dummy variables */

  ddir = (double *) malloc(MAXDATA * sizeof(double));
  dip = (double *) malloc(MAXDATA * sizeof(double));
  rake = (double *) malloc(MAXDATA * sizeof(double));
  x = (int *) malloc(MAXDATA * sizeof(int));
  y = (int *) malloc(MAXDATA * sizeof(int));
  dep = (int *) malloc(MAXDATA * sizeof(int));
  time = (int *) malloc(MAXDATA * sizeof(int));
  n = 3 * MAXDATA + 25 * MAXBOX;
  diag_sa = (double *) malloc(n * sizeof(double));
  diag_ija = (int *) malloc(n * sizeof(int));
  n = 18 * MAXDATA + 5 * MAXBOX + 1;
  amat_sa = (double *) malloc(n * sizeof(double));
  amat_ija = (int *) malloc(n * sizeof(int));
  n = 60 * (MAXBOX - SRMAXBOX) + 1;
  d_sa = (double *) malloc(n * sizeof(double));
  d_ija = (int *) malloc(n * sizeof(int));
  n = 3 * MAXDATA;
  slick = (double *) malloc(n * sizeof(double));
  slick_pre = (double *) malloc(n * sizeof(double));
  lindex = (int *) malloc(MAXBOX * sizeof(int));

  /* get file pointers */
  --argc;
  ++argv;
  if (argc != 4) /* [PM 11.04.2013] Corrected from 3 to 4, from < to != */
    {
    printf("usage: satsi_4D.exe data_file outfile damping time/space_damping\n");
    return -4001;/* [PM 11.04.2013] Changed from -1 to -4001*/
    }
  fpin = fopen(*argv, "r");
  if (fpin == NULL )
    {
    printf("unable to open %s.\n", *argv);
    return -4002; /* [PM 11.04.2013] Changed from -2 to -4002*/
    }
  ++argv;
  fpout = fopen(*argv, "a");
  if (fpout == NULL )
    {
    printf("unable to open %s.\n", *argv);
    return -4003; /* [PM 11.04.2013] Changed from -3 to -4003*/
    }
  ++argv;
  sscanf(*argv, "%lf", &cwt);
  ++argv;
  sscanf(*argv, "%lf", &twt);

  fgets(line, 80, fpin);

  for (i = 0; i < 3 * MAXDATA; i++)
    slick[i] = 0;
  for (i = 0; i < 3 * MAXDATA + 25 * MAXBOX; i++)
    {
    diag_sa[i] = 0;
    diag_ija[i] = 0;
    }

  /* loop to get data and make up equation */
  nobs = 0;
  nloc = 0;
  index = 0;
  while (fscanf(fpin, "%d %d %d %d %lf %lf %lf", &x[nobs], &y[nobs], &dep[nobs],
      &time[nobs], &ddir[nobs], &dip[nobs], &rake[nobs]) != EOF)
    {
    j = 3 * nobs;
    z = ddir[nobs] / TODEG;
    z2 = dip[nobs] / TODEG;
    z3 = rake[nobs] / TODEG;

    n1 = sin(z) * sin(z2); /* normal vector to fault plane */
    n2 = cos(z) * sin(z2);
    n3 = cos(z2);

    // [GK 2013.03.08] Commented as of no use.
    //norm[nobs][0] = n1;
    //norm[nobs][1] = n2;
    //norm[nobs][2] = n3;

    /* slickenside vector calculation */
    slick[j] = -cos(z3) * cos(z) - sin(z3) * sin(z) * cos(z2);
    slick[j + 1] = cos(z3) * sin(z) - sin(z3) * cos(z) * cos(z2);
    slick[j + 2] = sin(z3) * sin(z2);

    /* find the matrix elements */
    k = -1;
    for (i = 0; i < nloc; i++)
      if (lindex[i]
          == MAXY * MAXZ * MAXT * x[nobs] + MAXZ * MAXT * y[nobs]
              + MAXT * dep[nobs] + time[nobs])
        k = 5 * i;
    if (k == -1)
      {
      lindex[nloc] = MAXY * MAXZ * MAXT * x[nobs] + MAXZ * MAXT * y[nobs]
          + MAXT * dep[nobs] + time[nobs];
      k = 5 * nloc;
      nloc++;
      }

    temp[0] = n1 - n1 * n1 * n1 + n1 * n3 * n3;
    temp[1] = n2 - 2. * n1 * n1 * n2;
    temp[2] = n3 - 2. * n1 * n1 * n3;
    temp[3] = -n1 * n2 * n2 + n1 * n3 * n3;
    temp[4] = -2. * n1 * n2 * n3;
    diag_ija[j] = index;
    for (i = 0; i < 5; i++)
      {
      if ((k + i) == j)
        diag_sa[j] = temp[i];
      else
        {
        amat_ija[index] = k + i;
        amat_sa[index] = temp[i];
        index++;
        }
      }

    temp[0] = -n2 * n1 * n1 + n2 * n3 * n3;
    temp[1] = n1 - 2. * n1 * n2 * n2;
    temp[2] = -2. * n1 * n2 * n3;
    temp[3] = n2 - n2 * n2 * n2 + n2 * n3 * n3;
    temp[4] = n3 - 2. * n2 * n2 * n3;
    diag_ija[j + 1] = index;
    for (i = 0; i < 5; i++)
      {
      if ((k + i) == (j + 1))
        diag_sa[j + 1] = temp[i];
      else
        {
        amat_ija[index] = k + i;
        amat_sa[index] = temp[i];
        index++;
        }
      }

    temp[0] = -n3 * n1 * n1 - n3 + n3 * n3 * n3;
    temp[1] = -2. * n1 * n2 * n3;
    temp[2] = n1 - 2. * n1 * n3 * n3;
    temp[3] = -n3 * n2 * n2 - n3 + n3 * n3 * n3;
    temp[4] = n2 - 2. * n2 * n3 * n3;
    diag_ija[j + 2] = index;
    for (i = 0; i < 5; i++)
      {
      if ((k + i) == (j + 2))
        diag_sa[j + 2] = temp[i];
      else
        {
        amat_ija[index] = k + i;
        amat_sa[index] = temp[i];
        index++;
        }
      }

    ++nobs;
    /* check to see if all possible data has been read */
    if (nobs == MAXDATA)
      {
      printf("NOT ALL DATA COULD BE READ.\n");
      break;
      }
    } /* end of data read loop */

  nlocfill = nloc;

  /* fill in holes in grid */
  for (nz = 0; nz < MAXZ; nz++)
    for (nt = 0; nt < MAXT; nt++)
      for (nx = 0; nx < MAXX; nx++)
        {
        i = MAXY;
        j = 0;
        for (ny = 0; ny < MAXY; ny++)
          {
          for (i2 = 0; i2 < nloc; i2++)
            {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (ny < i))
              i = ny;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (ny > j))
              j = ny;
            }
          }
        if (i < j)
          for (ny = i + 1; ny < j; ny++)
            {
            ip = 0;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
                      + nt)
                ip = 1;
            if (ip == 0)
              {
              lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                  + MAXT * nz + nt;
              nloc++;
              }
            }
        }
  for (nz = 0; nz < MAXZ; nz++)
    for (nt = 0; nt < MAXT; nt++)
      for (ny = 0; ny < MAXY; ny++)
        {
        i = MAXX;
        j = 0;
        for (nx = 0; nx < MAXX; nx++)
          {
          for (i2 = 0; i2 < nloc; i2++)
            {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nx < i))
              i = nx;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nx > j))
              j = nx;
            }
          }
        if (i < j)
          for (nx = i + 1; nx < j; nx++)
            {
            ip = 0;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
                      + nt)
                ip = 1;
            if (ip == 0)
              {
              lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                  + MAXT * nz + nt;
              nloc++;
              }
            }
        }
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      for (nt = 0; nt < MAXT; nt++)
        {
        i = MAXZ;
        j = 0;
        for (nz = 0; nz < MAXZ; nz++)
          {
          for (i2 = 0; i2 < nloc; i2++)
            {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nz < i))
              i = nz;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nz > j))
              j = nz;
            }
          }
        if (i < j)
          for (nz = i + 1; nz < j; nz++)
            {
            ip = 0;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
                      + nt)
                ip = 1;
            if (ip == 0)
              {
              lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                  + MAXT * nz + nt;
              nloc++;
              }
            }
        }
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      for (nz = 0; nz < MAXZ; nz++)
        {
        i = MAXT;
        j = 0;
        for (nt = 0; nt < MAXT; nt++)
          {
          for (i2 = 0; i2 < nloc; i2++)
            {
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nt < i))
              i = nt;
            if ((lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
                && (nt > j))
              j = nt;
            }
          }
        if (i < j)
          for (nt = i + 1; nt < j; nt++)
            {
            ip = 0;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz
                      + nt)
                ip = 1;
            if (ip == 0)
              {
              lindex[nloc] = MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                  + MAXT * nz + nt;
              nloc++;
              }
            }
        }

  /* fill in diagonal */
  nrows = 3 * nobs;
  m = 5 * nloc;
  if (nrows < m)
    {
    for (i = nrows; i < m; i++)
      diag_ija[i] = index;
    nrows = m;
    }
  for (i = index - 1; i >= 0; i--)
    {
    amat_ija[i + nrows + 1] = amat_ija[i];
    amat_sa[i + nrows + 1] = amat_sa[i];
    }
  for (i = 0; i < nrows; i++)
    {
    amat_ija[i] = diag_ija[i] + nrows + 1;
    amat_sa[i] = diag_sa[i];
    }
  amat_ija[nrows] = index + nrows + 1;
  amat_sa[nrows] = 0;
  n = nrows;

  /* set up smoothing constraints */
  for (i = 0; i < 3 * MAXDATA + 25 * MAXBOX; i++)
    {
    diag_sa[i] = 0;
    diag_ija[i] = 0;
    }
  index = 0;
  j = 0;
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      for (nz = 0; nz < MAXZ; nz++)
        for (nt = 0; nt < MAXT; nt++)
          {
          k = -1;
          for (i2 = 0; i2 < nloc; i2++)
            if (lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
              k = 5 * i2;
          if (k > -1)
            {
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * (nx + 1) + MAXZ * MAXT * ny
                      + MAXT * nz + nt)
                k2 = 5 * i2;
            if ((nx < MAXX - 1) && (k2 > -1))
              {
              for (i = 0; i < 5; i++)
                {
                diag_ija[j + i] = index;
                if ((k + i) == (j + i))
                  diag_sa[j + i] = 1;
                else
                  {
                  d_ija[index] = k + i;
                  d_sa[index] = 1;
                  index++;
                  }
                if ((k2 + i) == (j + i))
                  diag_sa[j + i] = -1;
                else
                  {
                  d_ija[index] = k2 + i;
                  d_sa[index] = -1;
                  index++;
                  }
                }
              j += 5;
              }
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * (ny + 1)
                      + MAXT * nz + nt)
                k2 = 5 * i2;
            if ((ny < MAXY - 1) && (k2 > -1))
              {
              for (i = 0; i < 5; i++)
                {
                diag_ija[j + i] = index;
                if ((k + i) == (j + i))
                  diag_sa[j + i] = 1;
                else
                  {
                  d_ija[index] = k + i;
                  d_sa[index] = 1;
                  index++;
                  }
                if ((k2 + i) == (j + i))
                  diag_sa[j + i] = -1;
                else
                  {
                  d_ija[index] = k2 + i;
                  d_sa[index] = -1;
                  index++;
                  }
                }
              j += 5;
              }
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny
                      + MAXT * (nz + 1) + nt)
                k2 = 5 * i2;
            if ((nz < MAXZ - 1) && (k2 > -1))
              {
              for (i = 0; i < 5; i++)
                {
                diag_ija[j + i] = index;
                if ((k + i) == (j + i))
                  diag_sa[j + i] = 1;
                else
                  {
                  d_ija[index] = k + i;
                  d_sa[index] = 1;
                  index++;
                  }
                if ((k2 + i) == (j + i))
                  diag_sa[j + i] = -1;
                else
                  {
                  d_ija[index] = k2 + i;
                  d_sa[index] = -1;
                  index++;
                  }
                }
              j += 5;
              }
            k2 = -1;
            for (i2 = 0; i2 < nloc; i2++)
              if (lindex[i2]
                  == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt
                      + 1)
                k2 = 5 * i2;
            if ((nt < MAXT - 1) && (k2 > -1))
              {
              for (i = 0; i < 5; i++)
                {
                diag_ija[j + i] = index;
                if ((k + i) == (j + i))
                  diag_sa[j + i] = twt;
                else
                  {
                  d_ija[index] = k + i;
                  d_sa[index] = twt;
                  index++;
                  }
                if ((k2 + i) == (j + i))
                  diag_sa[j + i] = -twt;
                else
                  {
                  d_ija[index] = k2 + i;
                  d_sa[index] = -twt;
                  index++;
                  }
                }
              j += 5;
              }
            }
          }
  nrows = j;
  if (5 * nloc > nrows)
    {
    for (i = nrows; i < 5 * nloc; i++)
      diag_ija[i] = index;
    nrows = 5 * nloc;
    }
  for (i = index - 1; i >= 0; i--)
    {
    d_ija[i + nrows + 1] = d_ija[i];
    d_sa[i + nrows + 1] = d_sa[i];
    }
  for (i = 0; i < nrows; i++)
    {
    d_ija[i] = diag_ija[i] + nrows + 1;
    d_sa[i] = diag_sa[i];
    }
  d_ija[nrows] = index + nrows + 1;
  d_sa[nrows] = 0;
  for (i = 0; i < d_ija[d_ija[0] - 1]; i++)
    d_sa[i] = cwt * cwt * d_sa[i];
  p = nrows;

  /* solve equations via linear least squares */
  leasq_sparse(amat_ija, amat_sa, d_ija, d_sa, m, n, p, stress, slick);

  fprintf(fpout, "\nCOORDINATES ARE EAST,NORTH,UP.\n");
  fprintf(fpout, "stress tensors are:\n");
  fprintf(fpout, "X Y Z T See Sen Seu Snn Snu Suu\n");

  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      for (nz = 0; nz < MAXZ; nz++)
        for (nt = 0; nt < MAXT; nt++)
          {
          k = -1;
          for (i2 = 0; i2 < nlocfill; i2++)
            if (lindex[i2]
                == MAXY * MAXZ * MAXT * nx + MAXZ * MAXT * ny + MAXT * nz + nt)
              k = 5 * i2;
          if (k > -1)
            {
            fprintf(fpout, "%3d %3d %3d %3d ", nx, ny, nz, nt);
            fprintf(fpout, "%9.6f %9.6f %9.6f ", stress[k], stress[k + 1],
                stress[k + 2]);
            fprintf(fpout, "%9.6f %9.6f %9.6f\n", stress[k + 3], stress[k + 4],
                -(stress[k] + stress[k + 3]));
            }
          }
  fprintf(fpout, "\n");

  fprintf(fpout,
      "\ndip direction, dip, rake, fit angle, mag tau, X, Y, Z, T\n");
  sprsax(amat_sa, amat_ija, stress, slick_pre, m, n);
  angavg = 0.;
  angstd = 0.;
  magavg = 0.;
  magstd = 0.;
  for (i = 0; i < nobs; i++)
    {
    ls1 = sqrt(
        slick[3 * i] * slick[3 * i] + slick[3 * i + 1] * slick[3 * i + 1]
            + slick[3 * i + 2] * slick[3 * i + 2]);
    ls2 = sqrt(
        slick_pre[3 * i] * slick_pre[3 * i]
            + slick_pre[3 * i + 1] * slick_pre[3 * i + 1]
            + slick_pre[3 * i + 2] * slick_pre[3 * i + 2]);
    ls3 = slick[3 * i] * slick_pre[3 * i]
        + slick[3 * i + 1] * slick_pre[3 * i + 1]
        + slick[3 * i + 2] * slick_pre[3 * i + 2];
    z = TODEG * acos(ls3 / (ls1 * ls2));
    angavg += z;
    angstd += z * z;
    magavg += ls2;
    magstd += ls2 * ls2;
    fprintf(fpout, "%7.1f  %7.1f  %7.1f  %7.1f %7.2f %4d %4d %4d %4d\n",
        ddir[i], dip[i], rake[i], z, ls2, x[i], y[i], dep[i], time[i]);
    }
  z3 = (double) nobs - 1;
  angstd = angstd - (angavg * angavg / nobs);
  angstd = angstd / z3;
  angstd = sqrt(angstd);
  angavg = angavg / nobs;
  magstd = magstd - (magavg * magavg / nobs);
  magstd = magstd / z3;
  magstd = sqrt(magstd);
  magavg = magavg / nobs;
  fprintf(fpout, "\nfit angle mean= %f standard deviation= %f\n", angavg,
      angstd);
  fprintf(fpout, "avg tau= %f , std. dev.= %f\n", magavg, magstd);

  free(ddir);
  free(dip);
  free(rake);
  free(x);
  free(y);
  free(dep);
  free(time);
  free(diag_sa);
  free(diag_ija);
  free(amat_sa);
  free(amat_ija);
  free(d_sa);
  free(d_ija);
  free(slick);
  free(slick_pre);
  free(lindex);

  return 0; // [GK 2013.03.08] Default return value.
  }
