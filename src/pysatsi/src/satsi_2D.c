//-------------------------------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//-------------------------------------------------------------------------------------------------

/* COORDINATES ARE EAST,NORTH,UP */

//-------------------------------------------------------------------------------------------------
// SATSI_2D
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
// Lion Krischer <krischer@geophysik.uni-muenchen.de>, Dez. 2014:
//   * Refactoring to allow usage as a library.
//   * Removed separate reading of first line which most likely was a bug.
//
// $Last revision: 1.0 $  $Date: 2012/07/11  $
//-------------------------------------------------------------------------------------------------

#define TODEG 57.29577951
//-------------------------------------------------------------------------------------------------
// [GK 2013.05.15] Original SATSI setup.
//#define MAXDATA 7000
//#define MAXX 56
//#define MAXY 56
//#define MAXBOX 900
//#define SRMAXBOX 30
// [GK 2013.05.15] Standard version.
//#define MAXDATA 70000
//#define MAXX 50
//#define MAXY 50
//#define MAXBOX 10000
//#define SRMAXBOX 100
// [GK 2013.05.15] Extended SATSI setup.
#define MAXDATA 70000
#define MAXX 200
#define MAXY 200
#define MAXBOX 10000
#define SRMAXBOX 100
//-------------------------------------------------------------------------------------------------

// [GK 2013.03.03] Additional declarations to suppress warning messages.
void sprsax(double sa[], int ija[], double x[], double b[], int m, int n);
void leasq_sparse(int a_ija[], double a_sa[], int d_ija[], double d_sa[], int m,
    int n, int p, double x[], double b[]);


// 2D stress tensor.
struct stress_tensor_2D {
  int X;
  int Y;
  double See;
  double Sen;
  double Seu;
  double Snn;
  double Snu;
  double Suu;
};


// for each input focal mechanism this struct represents the magnitude of
// shear traction vector and deviation between predicted and calculated slip.
struct results_for_input_mt {
  double ddir;
  double dip;
  double rake;
  double fit_angle;
  double mag_tau;
  int X;
  int Y;
};


// Output structure.
struct satsi_2D_result {
  double fit_angle_mean; // Fit angle mean
  double fit_angle_std; // Fit angle standard deviation
  double mag_mean; // Tangential stress mean.
  double mag_std; // Tangential stress standard deviation.
  struct stress_tensor_2D *tensors;
  int tensor_count;
  struct results_for_input_mt *results;
  int result_count;
};


// Function to free the result struct.
void free_satsi_2D_result(struct satsi_2D_result *result) {
  free(result->tensors);
  free(result->results);
}


//-------------------------------------------------------------------------------------------------
  // Parameters:
  //   x_in: Integer array with the x coordinate of the values.
  //   y_in: Integer array with the y coordinate of the values.
  //   ddir_in: Double array containing the dip direction.
  //   dip_in: Double array containing the dip angle.
  //   rake_in: Double arrays containing the rake angles.
  //   input_length: The length of the arrays. All must have the same length!
  //   cwt: The damping parameter.
struct satsi_2D_result satsi_2D(
  int *x_in, int *y_in, double *ddir_in, double *dip_in, double *rake_in,
  int input_length, double cwt)
{
  double *ddir, *dip, *rake; /* focal mechanism data */
  int *x, *y; /* focal mechanism bins */
  int nobs, nloc, nrows, nlocfill; /* number of observations, bins, rows */
  double *diag_sa, *amat_sa, *d_sa; /* inversion matrices in sparse matrix, */
  int *diag_ija, *amat_ija, *d_ija; /*      row-indexed form */
  int *loclist, lindex, nx, ny, index; /* book-keeping, which bins where in matrix*/
  double stress[5 * MAXBOX]; /* stress tensor in vector form xx,xy,xz,yy,yz,zz */
  double *slick; /* slickenside vector elements vector */
  double n1, n2, n3; /* normal vector elements */
  //double norm[MAXDATA][3]; /* storage of n1,n2,n3 */ // GK: 2013.03.03 Variable not used.
  double angavg, angstd; /* average and standard deviation of fit angle */
  double magavg, magstd; /* same for tangential stress size */
  double *slick_pre; /* predicted slip vector */
  int i, j, k, k2, m, n, p; /* dummy variables */
  double z, z2, z3, temp[5], ls1, ls2, ls3; /* more dummy variables */

  // Current array index.
  int cur_idx = 0;
  struct satsi_2D_result result;
  int output_stress_tensor_count;

  ddir = (double *) malloc(MAXDATA * sizeof(double));
  dip = (double *) malloc(MAXDATA * sizeof(double));
  rake = (double *) malloc(MAXDATA * sizeof(double));
  x = (int *) malloc(MAXDATA * sizeof(int));
  y = (int *) malloc(MAXDATA * sizeof(int));
  n = 3 * MAXDATA + 25 * MAXBOX;
  diag_sa = (double *) malloc(n * sizeof(double));
  diag_ija = (int *) malloc(n * sizeof(int));
  n = 18 * MAXDATA + 5 * MAXBOX + 1;
  amat_sa = (double *) malloc(n * sizeof(double));
  amat_ija = (int *) malloc(n * sizeof(int));
  n = 30 * (MAXBOX - SRMAXBOX) + 1;
  d_sa = (double *) malloc(n * sizeof(double));
  d_ija = (int *) malloc(n * sizeof(int));
  n = 3 * MAXDATA;
  slick = (double *) malloc(n * sizeof(double));
  slick_pre = (double *) malloc(n * sizeof(double));
  loclist = (int *) malloc(MAXX * MAXY * sizeof(int));

  for (i = 0; i < 3 * MAXDATA; i++)
    slick[i] = 0;
  for (i = 0; i < 3 * MAXDATA + 25 * MAXBOX; i++)
    {
    diag_sa[i] = 0;
    diag_ija[i] = 0;
    }
  for (i = 0; i < MAXX; i++)
    for (j = 0; j < MAXY; j++)
      loclist[MAXY * i + j] = -999;

  /* loop to get data and make up equation */
  nobs = 0;
  nloc = 0;
  index = 0;

  while (cur_idx < input_length) {
    x[nobs] = x_in[cur_idx];
    y[nobs] = y_in[cur_idx];
    ddir[nobs] = ddir_in[cur_idx];
    dip[nobs] = dip_in[cur_idx];
    rake[nobs] = rake_in[cur_idx];
    cur_idx++;

    j = 3 * nobs;
    z = ddir[nobs] / TODEG;
    z2 = dip[nobs] / TODEG;
    z3 = rake[nobs] / TODEG;

    n1 = sin(z) * sin(z2); /* normal vector to fault plane */
    n2 = cos(z) * sin(z2);
    n3 = cos(z2);

    //norm[nobs][0] = n1; [GK 2012.07.10] Suppressed as it is not used.
    //norm[nobs][1] = n2;
    //norm[nobs][2] = n3;

    /* slickenside vector calculation */
    slick[j] = -cos(z3) * cos(z) - sin(z3) * sin(z) * cos(z2);
    slick[j + 1] = cos(z3) * sin(z) - sin(z3) * cos(z) * cos(z2);
    slick[j + 2] = sin(z3) * sin(z2);

    /* find the matrix elements */
    lindex = MAXY * x[nobs] + y[nobs];
    if (loclist[lindex] > -999)
      k = 5 * loclist[lindex];
    else
      {
      loclist[lindex] = nloc;
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
  for (nx = 0; nx < MAXX; nx++)
    {
    i = 0;
    while ((i < MAXY) && (loclist[MAXY * nx + i] == -999))
      i++;
    j = MAXY - 1;
    while ((j > -1) && (loclist[MAXY * nx + j] == -999))
      j--;
    if (i < j)
      for (ny = i + 1; ny < j; ny++)
        if (loclist[MAXY * nx + ny] == -999)
          {
          loclist[MAXY * nx + ny] = nloc;
          nloc++;
          }
    }
  for (ny = 0; ny < MAXY; ny++)
    {
    i = 0;
    while ((i < MAXX) && (loclist[MAXY * i + ny] == -999))
      i++;
    j = MAXX - 1;
    while ((j > -1) && (loclist[MAXY * j + ny] == -999))
      j--;
    if (i < j)
      for (nx = i + 1; nx < j; nx++)
        if (loclist[MAXY * nx + ny] == -999)
          {
          loclist[MAXY * nx + ny] = nloc;
          nloc++;
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
      if (loclist[MAXY * nx + ny] > -999)
        {
        k = 5 * loclist[MAXY * nx + ny];
        if ((nx < MAXX - 1) && (loclist[MAXY * (nx + 1) + ny] > -999))
          {
          k2 = 5 * loclist[MAXY * (nx + 1) + ny];
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
        if ((ny < MAXY - 1) && (loclist[MAXY * nx + ny + 1] > -999))
          {
          k2 = 5 * loclist[MAXY * nx + ny + 1];
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

  output_stress_tensor_count = 0;
  /* Count the amount of output stress tensors */
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      if ((loclist[MAXY * nx + ny] > -999)
          && (loclist[MAXY * nx + ny] < nlocfill))
        {
          output_stress_tensor_count++;
        }

  // Allocate space for the output tensors.
  result.tensors = (struct stress_tensor_2D *) malloc (output_stress_tensor_count * sizeof(struct stress_tensor_2D));
  result.tensor_count = output_stress_tensor_count;


  // Write the output tensors.
  cur_idx = 0;
  for (nx = 0; nx < MAXX; nx++)
    for (ny = 0; ny < MAXY; ny++)
      if ((loclist[MAXY * nx + ny] > -999)
          && (loclist[MAXY * nx + ny] < nlocfill))
        {
        k = 5 * loclist[MAXY * nx + ny];

        result.tensors[cur_idx].X = nx;
        result.tensors[cur_idx].Y = ny;
        result.tensors[cur_idx].See = stress[k];
        result.tensors[cur_idx].Sen = stress[k + 1];
        result.tensors[cur_idx].Seu = stress[k + 2];
        result.tensors[cur_idx].Snn = stress[k + 3];
        result.tensors[cur_idx].Snu = stress[k + 4];
        result.tensors[cur_idx].Suu = -(stress[k] + stress[k + 3]);

        cur_idx++;
        }

  sprsax(amat_sa, amat_ija, stress, slick_pre, m, n);
  angavg = 0.;
  angstd = 0.;
  magavg = 0.;
  magstd = 0.;


  result.results = (struct results_for_input_mt *) malloc (nobs * sizeof(struct results_for_input_mt));
  result.result_count = nobs;

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

    result.results[i].ddir = ddir[i];
    result.results[i].dip = dip[i];
    result.results[i].rake = rake[i];
    result.results[i].fit_angle = z;
    result.results[i].mag_tau = ls2;
    result.results[i].X = x[i];
    result.results[i].Y = y[i];
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

  result.fit_angle_mean = angavg;
  result.fit_angle_std = angstd;
  result.mag_mean = magavg;
  result.mag_std = magstd;

  free(ddir);
  free(dip);
  free(rake);
  free(x);
  free(y);
  free(diag_sa);
  free(diag_ija);
  free(amat_sa);
  free(amat_ija);
  free(d_sa);
  free(d_ija);
  free(slick);
  free(slick_pre);
  free(loclist);

  return result;
  }
