#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1986-2004 by the University of Washington
  and by Joseph Felsenstein.  Written by Joseph Felsenstein.  Permission is
  granted to copy and use this program provided no fee is charged for it
  and provided that this copyright notice is not removed. */

#define over            60       /* Maximum xcoord of tip nodes */

/* These are redefined from phylip.h */
#undef epsilon
#undef smoothings
#define epsilon         0.000005   /* used during branch length
                                    * optimization. Smaller values
                                    * yield better estimates, but
                                    * take longer. */
#define smoothings      4

/* DEBUG: The following will print many details about node tyme optimization */
/* #define MAKENEWV_DEBUG */

/* DEBUG: This does extra consistency checks in nuview() */
/* #define NUVIEW_DEBUG */

const double MIN_BRANCH_LENGTH = epsilon/4.0;
const double MIN_ROOT_TYME = -10;

typedef struct valrec {
  double rat, ratxi, ratxv, orig_zz, z1, y1, z1zz, z1yy, xiz1,
         xiy1xv;
  double *ww, *zz, *wwzz, *vvzz; 
} valrec;

typedef double contribarr[maxcategs];

valrec ***tbl;

#ifndef OLDC
/* function prototypes */
void    getoptions(void);
void    allocrest(void);
void    doinit(void);
void    inputoptions(void);
void    makeweights(void);
void    getinput(void);
void    inittable_for_usertree (FILE *);
void    inittable(void);
void    exmake(double, long);
void    alloc_nvd(long, nuview_data *);
void    free_nvd(nuview_data *);
node   *invalid_descendant_view(node *);
boolean nuview(node *);
double  evaluate(node *);
node   *pnode(tree *, node *);
double  min_child_tyme(node *);
double  parent_tyme(node *);
boolean all_tymes_valid(node *, double, boolean);
boolean valid_tyme(node *, double);
void    invalidate_tyme(node *);
void    invalidate_traverse(node *);
double  set_tyme_evaluate(node *, double);
void    makenewv(node *);
void    update(node *);
void    smooth(node *);
void    restoradd(node *, node *, node *, double);
void    dnamlk_add(node *, node *, node *);
void    dnamlk_re_move(node **, node **, boolean);
void    tryadd(node *, node **, node **);
void    addpreorder(node *, node *, node *, boolean, boolean);
boolean tryrearr(node *);
boolean repreorder(node *);
void    rearrange(node **);
void    initdnamlnode(node **, node **, node *, long, long, long *, long *,
                initops, pointarray, pointarray, Char *, Char *, FILE *);
void    tymetrav(node *, double *);
void    dnamlk_coordinates(node *, long *);
void    dnamlk_drawline(long, double);
void    dnamlk_printree(void);
void    describe(node *);
void    reconstr(node *, long);
void    rectrav(node *, long, long);
void    summarize(void);
void    dnamlk_treeout(node *);
void    nodeinit(node *);
void    initrav(node *);
void    travinit(node *);
void    travsp(node *);
void    treevaluate(void);
void    maketree(void);
void    reallocsites(void);
void    save_tree_tyme(tree *, double[]);
void    restore_saved_tyme(tree *, double[]);
void    freetable(void);
void    setnodetymes(node *, double);
/* function prototypes */
#endif


Char infilename[FNMLNGTH], outfilename[FNMLNGTH], intreename[FNMLNGTH], outtreename[FNMLNGTH],
     catfilename[FNMLNGTH], weightfilename[FNMLNGTH];
double *rrate;
long sites, weightsum, categs, datasets, ith, njumble, jumb, numtrees, shimotrees;
/*  sites = number of sites in actual sequences
  numtrees = number of user-defined trees */
long inseed, inseed0, mx, mx0, mx1;
boolean freqsfrom, global, global2=0, jumble, trout, usertree, weights, 
          rctgry, ctgry, ttr, auto_, progress, mulsets, firstset, hypstate, 
          smoothit, polishing, justwts, gama, invar;
boolean lengthsopt = false;             /* Use lengths in user tree option */
boolean lngths = false;                 /* Actually use lengths (depends on
                                           each input tree) */
tree curtree, bestree, bestree2, priortree;
node *qwhere, *grbg;
double *tymes;
double xi, xv, ttratio, ttratio0, freqa, freqc, freqg, freqt, freqr,
              freqy, freqar, freqcy, freqgr, freqty, fracchange, sumrates,
              cv, alpha, lambda, lambda1, invarfrac;
long *enterorder;
steptr aliasweight;
double *rate;
longer seed;
double *probcat;
contribarr *contribution;
char *progname;
long rcategs, nonodes2;
long **mp;
char basechar[16]="acmgrsvtwyhkdbn";


/* Local variables for maketree, propagated globally for C version: */
long    k, maxwhich, col;
double  like, bestyet, maxlogl;
boolean lastsp;

boolean smoothed; /* set true before each smoothing run, and set false
                     each time a branch cannot be completely optimized. */
double  *l0gl;
double  expon1i[maxcategs], expon1v[maxcategs],
        expon2i[maxcategs], expon2v[maxcategs];
node   *there;
double  **l0gf;
Char ch, ch2;


void save_tree_tyme(tree* save_tree, double tymes[])
{
  int i;
  assert( all_tymes_valid(curtree.root, 0.0, false) );
  for ( i = spp ; i < nonodes ; i++) {
    tymes[i - spp] = save_tree->nodep[i]->tyme;
  }
}

void restore_saved_tyme(tree *load_tree, double tymes[])
{
  int i;
  for ( i = spp ; i < nonodes ; i++) {
    setnodetymes(load_tree->nodep[i], tymes[i-spp]);
  }
  /* Fix invalid tymes */
  all_tymes_valid(curtree.root, 0.0, true);
}

void getoptions()
{
  /* interactively set options */
  long i, loopcount, loopcount2;
  Char ch;
  boolean done;
  boolean didchangecat, didchangercat;
  double probsum;

  fprintf(outfile, "\nNucleic acid sequence\n");
  fprintf(outfile, "   Maximum Likelihood method with molecular ");
  fprintf(outfile, "clock, version %s\n\n", VERSION);
  putchar('\n');
  auto_ = false;
  ctgry = false;
  didchangecat = false;
  rctgry = false;
  didchangercat = false;
  categs = 1;
  rcategs = 1;
  freqsfrom = true;
  gama = false;
  invar = false;
  global = false;
  hypstate = false;
  jumble = false;
  njumble = 1;
  lambda = 1.0;
  lambda1 = 0.0;
  lengthsopt = false;
  trout = true;
  ttratio = 2.0;
  ttr = false;
  usertree = false;
  weights = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  treeprint = true;
  interleaved = true;
  loopcount = 0;
  do {
    cleerhome();
    printf("\nNucleic acid sequence\n");
    printf("   Maximum Likelihood method with molecular clock, version %s\n\n",
           VERSION);
    printf("Settings for this run:\n");
    printf("  U                 Search for best tree?");
    if (usertree)
     printf("  No, use user trees in input file\n");
    else
      printf("  Yes\n");
    if (usertree) {
      printf("  L           Use lengths from user tree?");
      if (lengthsopt)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  T        Transition/transversion ratio:");
    if (!ttr)
      printf("  2.0\n");
    else
      printf("  %8.4f\n", ttratio);
    printf("  F       Use empirical base frequencies?");
    if (freqsfrom)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  C   One category of substitution rates?");
    if (!ctgry)
      printf("  Yes\n");
    else
      printf("  %ld categories\n", categs);
    printf("  R           Rate variation among sites?");
    if (!rctgry)
      printf("  constant rate\n");
    else {
      if (gama)
        printf("  Gamma distributed rates\n");
      else {
        if (invar)
          printf("  Gamma+Invariant sites\n");
        else
          printf("  user-defined HMM of rates\n");
      }
      printf("  A   Rates at adjacent sites correlated?");
      if (!auto_)
        printf("  No, they are independent\n");
      else
        printf("  Yes, mean block length =%6.1f\n", 1.0 / lambda);
    }
    if (!usertree) {
      printf("  G                Global rearrangements?");
      if (global)
        printf("  Yes\n");
      else
        printf("  No\n");
    }
    printf("  W                       Sites weighted?  %s\n",
           (weights ? "Yes" : "No"));
    if (!usertree) {
      printf("  J   Randomize input order of sequences?");
      if (jumble)
        printf("  Yes (seed = %8ld, %3ld times)\n", inseed0, njumble);
      else
        printf("  No. Use input order\n");
    }
    printf("  M           Analyze multiple data sets?");
    if (mulsets)
      printf("  Yes, %2ld %s\n", datasets,
               (justwts ? "sets of weights" : "data sets"));
    else
      printf("  No\n");
    printf("  I          Input sequences interleaved?");
    if (interleaved)
      printf("  Yes\n");
    else
      printf("  No, sequential\n");
    printf("  0   Terminal type (IBM PC, ANSI, none)?");
    if (ibmpc)
      printf("  IBM PC\n");
    if (ansi)
      printf("  ANSI\n");
    if (!(ibmpc || ansi))
      printf("  (none)\n");
    printf("  1    Print out the data at start of run");
    if (printdata)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  2  Print indications of progress of run");
    if (progress)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  3                        Print out tree");
    if (treeprint)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  4       Write out trees onto tree file?");
    if (trout)
      printf("  Yes\n");
    else
      printf("  No\n");
    printf("  5   Reconstruct hypothetical sequences?  %s\n",
           (hypstate ? "Yes" : "No"));
    printf("\nAre these settings correct? (type Y or the letter for one to change)\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    if (ch == '\n')
      ch = ' ';
    uppercase(&ch);
    done = (ch == 'Y');
    if (!done) {
      uppercase(&ch);
      if (((!usertree) && (strchr("JUCRAFWGTMI012345", ch) != NULL))
          || (usertree && ((strchr("UCRAFWLTMI012345", ch) != NULL)))){
        switch (ch) {

        case 'C':
          ctgry = !ctgry;
          if (ctgry) {
            printf("\nSitewise user-assigned categories:\n\n");
            initcatn(&categs);
            if (rate){
              free(rate);
            }
            rate    = (double *) Malloc( categs * sizeof(double));
            didchangecat = true;
            initcategs(categs, rate);
          }
          break;

        case 'R':
          if (!rctgry) {
            rctgry = true;
            gama = true;
          } else {
            if (gama) {
              gama = false;
              invar = true;
            } else {
              if (invar)
                invar = false;
              else
                rctgry = false;
            }
          }
          break;

        case 'A':
          auto_ = !auto_;
          if (auto_) {
            initlambda(&lambda);
            lambda1 = 1.0 - lambda;
          }
          break;

        case 'F':
          freqsfrom = !freqsfrom;
          if (!freqsfrom)
            initfreqs(&freqa, &freqc, &freqg, &freqt);
          break;

        case 'G':
          global = !global;
          break;

        case 'W':
          weights = !weights;
          break;

        case 'J':
          jumble = !jumble;
          if (jumble)
            initjumble(&inseed, &inseed0, seed, &njumble);
          else njumble = 1;
          break;

        case 'L':
          lengthsopt = !lengthsopt;
          break;

        case 'T':
          ttr = !ttr;
          if (ttr)
            initratio(&ttratio);
          break;

        case 'U':
          usertree = !usertree;
          break;

        case 'M':
          mulsets = !mulsets;
          if (mulsets) {
            printf("Multiple data sets or multiple weights?");
            loopcount2 = 0;
            do {
              printf(" (type D or W)\n");
#ifdef WIN32
              phyFillScreenColor();
#endif
              scanf("%c%*[^\n]", &ch2);
              getchar();
              if (ch2 == '\n')
                ch2 = ' ';
              uppercase(&ch2);
              countup(&loopcount2, 10);
            } while ((ch2 != 'W') && (ch2 != 'D'));
            justwts = (ch2 == 'W');
            if (justwts)
              justweights(&datasets);
            else
              initdatasets(&datasets);
            if (!jumble) {
              jumble = true;
              initjumble(&inseed, &inseed0, seed, &njumble);
            }
          }
          break;

        case 'I':
          interleaved = !interleaved;
          break;

        case '0':
          initterminal(&ibmpc, &ansi);
          break;

        case '1':
          printdata = !printdata;
          break;

        case '2':
          progress = !progress;
          break;

        case '3':
          treeprint = !treeprint;
          break;

        case '4':
          trout = !trout;
          break;

        case '5':
          hypstate = !hypstate;
          break;
        }
      } else
        printf("Not a possible option!\n");
    }
    countup(&loopcount, 100);
  } while (!done);
  if (gama || invar) {
    loopcount = 0;
    do {
      printf(
"\nCoefficient of variation of substitution rate among sites (must be positive)\n");
      printf(
  " In gamma distribution parameters, this is 1/(square root of alpha)\n");
#ifdef WIN32
      phyFillScreenColor();
#endif
      scanf("%lf%*[^\n]", &cv);
      getchar();
      countup(&loopcount, 10);
    } while (cv <= 0.0);
    alpha = 1.0 / (cv * cv);
  }
  if (!rctgry)
    auto_ = false;
  if (rctgry) {
    printf("\nRates in HMM");
    if (invar)
      printf(" (including one for invariant sites)");
    printf(":\n");
    initcatn(&rcategs);
    if (probcat){
      free(probcat);
      free(rrate);
    }
    probcat = (double *) Malloc(rcategs * sizeof(double));
    rrate   = (double *) Malloc(rcategs * sizeof(double));
    didchangercat = true;
    if (gama)
      initgammacat(rcategs, alpha, rrate, probcat); 
    else {
      if (invar) {
        loopcount = 0;
        do {
          printf("Fraction of invariant sites?\n");
          scanf("%lf%*[^\n]", &invarfrac);
          getchar();
          countup(&loopcount, 10);
        } while ((invarfrac <= 0.0) || (invarfrac >= 1.0));
        initgammacat(rcategs-1, alpha, rrate, probcat); 
        for (i = 0; i < rcategs-1; i++)
          probcat[i] = probcat[i]*(1.0-invarfrac);
        probcat[rcategs-1] = invarfrac;
        rrate[rcategs-1] = 0.0;
      } else {
        initcategs(rcategs, rrate);
        initprobcat(rcategs, &probsum, probcat);
      }
    }
  }
  if (!didchangercat){
    rrate      = Malloc( rcategs*sizeof(double));
    probcat    = Malloc( rcategs*sizeof(double));
    rrate[0]   = 1.0;
    probcat[0] = 1.0;
  }
  if (!didchangecat){
    rate       = Malloc( categs*sizeof(double));
    rate[0]    = 1.0;
  }
}  /* getoptions */

void reallocsites(void) 
{
  long i;

  for (i = 0; i < spp; i++) {
    free(y[i]);
    y[i] = (char *)Malloc(sites * sizeof(char));
  }
  free(weight);
  free(category);
  free(alias);
  free(aliasweight);
  free(ally);
  free(location);
  
  weight      = (long *)Malloc(sites*sizeof(long));
  category    = (long *)Malloc(sites*sizeof(long));
  alias       = (long *)Malloc(sites*sizeof(long));
  aliasweight = (long *)Malloc(sites*sizeof(long));
  ally        = (long *)Malloc(sites*sizeof(long));
  location    = (long *)Malloc(sites*sizeof(long));
} /* reallocsites */


void allocrest()
{
  long i;

  y     = (Char **)Malloc(spp*sizeof(Char *));
  nayme  = (naym *)Malloc(spp*sizeof(naym));
  for (i = 0; i < spp; i++)
    y[i] = (char *)Malloc(sites * sizeof(char));
  enterorder  = (long *)Malloc(spp*sizeof(long));
  weight      = (long *)Malloc(sites*sizeof(long));
  category    = (long *)Malloc(sites*sizeof(long));
  alias       = (long *)Malloc(sites*sizeof(long));
  aliasweight = (long *)Malloc(sites*sizeof(long));
  ally        = (long *)Malloc(sites*sizeof(long));
  location    = (long *)Malloc(sites*sizeof(long));
  tymes       = (double *)Malloc((nonodes - spp) * sizeof(double));
}  /* allocrest */


void doinit()
{
  /* initializes variables */

  inputnumbers(&spp, &sites, &nonodes, 1);
  getoptions();
  if (printdata)
    fprintf(outfile, "%2ld species, %3ld  sites\n", spp, sites);
  alloctree(&curtree.nodep, nonodes, usertree);
  allocrest();
  if (usertree)
    return;
  alloctree(&bestree.nodep, nonodes, 0);
  if (njumble <= 1)
    return;
  alloctree(&bestree2.nodep, nonodes, 0);
}  /* doinit */


void inputoptions()
{
  long i;

  if (!firstset && !justwts) {
    samenumsp(&sites, ith);
    reallocsites();
  }

  for (i = 0; i < sites; i++)
    category[i] = 1;
  for (i = 0; i < sites; i++)
    weight[i] = 1;
  
  if (justwts || weights)
    inputweights(sites, weight, &weights);
  weightsum = 0;
  for (i = 0; i < sites; i++)
    weightsum += weight[i];
  if (ctgry && categs > 1) {
    inputcategs(0, sites, category, categs, "DnaMLK");
    if (printdata)
      printcategs(outfile, sites, category, "Site categories");
  }
  if (weights && printdata)
    printweights(outfile, 0, sites, weight, "Sites");
}  /* inputoptions */


void makeweights()
{
  /* make up weights vector to avoid duplicate computations */
  long i;

   for (i = 1; i <= sites; i++) {
    alias[i - 1] = i;
    ally[i - 1] = 0;
    aliasweight[i - 1] = weight[i - 1];
    location[i - 1] = 0;
  }
  sitesort2(sites, aliasweight);
  sitecombine2(sites, aliasweight);
  sitescrunch2(sites, 1, 2, aliasweight);
  for (i = 1; i <= sites; i++) {
    if (aliasweight[i - 1] > 0)
      endsite = i;
  }
  for (i = 1; i <= endsite; i++) {
    ally[alias[i - 1] - 1] = alias[i - 1];
    location[alias[i - 1] - 1] = i;
  }
  contribution = (contribarr *) Malloc( endsite*sizeof(contribarr));
}  /* makeweights */


void getinput()
{
  
  /* reads the input data */
  inputoptions();
  if (!freqsfrom)
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
                   &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
                   freqsfrom, true);
  if (!justwts || firstset)
    inputdata(sites);
  makeweights();
  setuptree2(&curtree);
  if (!usertree) {
    setuptree2(&bestree);
    if (njumble > 1)
      setuptree2(&bestree2);
  }
  allocx(nonodes, rcategs, curtree.nodep, usertree);
  if (!usertree) {
    allocx(nonodes, rcategs, bestree.nodep, 0);
    if (njumble > 1)
      allocx(nonodes, rcategs, bestree2.nodep, 0);
  }
  makevalues2(rcategs, curtree.nodep, endsite, spp, y, alias);
  if (freqsfrom) {
    empiricalfreqs(&freqa, &freqc, &freqg, &freqt, aliasweight, curtree.nodep);
    getbasefreqs(freqa, freqc, freqg, freqt, &freqr, &freqy, &freqar, &freqcy,
                   &freqgr, &freqty, &ttratio, &xi, &xv, &fracchange,
                   freqsfrom, true);
  }
  if (!justwts || firstset)
    fprintf(outfile, "\nTransition/transversion ratio = %10.6f\n\n", ttratio);
}  /* getinput */


void inittable_for_usertree (FILE *intree)
{
  /* If there's a user tree, then the ww/zz/wwzz/vvzz elements need
     to be allocated appropriately. */
  long num_comma;
  long i, j;

  /* First, figure out the largest possible furcation, i.e. the number
     of commas plus one */
  countcomma (&intree, &num_comma);
  num_comma++;
  
  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      /* Free the stuff allocated assuming bifurcations */
      free (tbl[i][j]->ww);
      free (tbl[i][j]->zz);
      free (tbl[i][j]->wwzz);
      free (tbl[i][j]->vvzz);

      /* Then allocate for worst-case multifurcations */
      tbl[i][j]->ww   = (double *) Malloc( num_comma * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc( num_comma * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc( num_comma * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc( num_comma * sizeof (double));
    }
  }
}  /* inittable_for_usertree */


void freetable()
{
  long i, j;
 
  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      free(tbl[i][j]->ww);
      free(tbl[i][j]->zz);
      free(tbl[i][j]->wwzz);
      free(tbl[i][j]->vvzz);
    }
  }
  for (i = 0; i < rcategs; i++)  {
    for (j = 0; j < categs; j++) 
      free(tbl[i][j]);
    free(tbl[i]);
  }
  free(tbl);
}

void inittable()
{
  /* Define a lookup table. Precompute values and print them out in tables */
  long i, j;
  double sumrates;
  
  tbl = (valrec ***) Malloc( rcategs * sizeof(valrec **));
  for (i = 0; i < rcategs; i++) {
    tbl[i] = (valrec **) Malloc( categs*sizeof(valrec *));
    for (j = 0; j < categs; j++)
      tbl[i][j] = (valrec *) Malloc( sizeof(valrec));
  }

  for (i = 0; i < rcategs; i++) {
    for (j = 0; j < categs; j++) {
      tbl[i][j]->rat = rrate[i]*rate[j];
      tbl[i][j]->ratxi = tbl[i][j]->rat * xi;
      tbl[i][j]->ratxv = tbl[i][j]->rat * xv;

      /* Allocate assuming bifurcations, will be changed later if
         neccesarry (i.e. there's a user tree) */
      tbl[i][j]->ww   = (double *) Malloc( 2 * sizeof (double));
      tbl[i][j]->zz   = (double *) Malloc( 2 * sizeof (double));
      tbl[i][j]->wwzz = (double *) Malloc( 2 * sizeof (double));
      tbl[i][j]->vvzz = (double *) Malloc( 2 * sizeof (double));
    }
  }
  sumrates = 0.0;
  for (i = 0; i < endsite; i++) {
    for (j = 0; j < rcategs; j++)
      sumrates += aliasweight[i] * probcat[j] 
        * tbl[j][category[alias[i] - 1] - 1]->rat;
  }
  sumrates /= (double)sites;
  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
      tbl[i][j]->rat /= sumrates;
      tbl[i][j]->ratxi /= sumrates;
      tbl[i][j]->ratxv /= sumrates;
    }

  if(jumb > 1)
    return;

  if (gama || invar) {
    fprintf(outfile, "\nDiscrete approximation to gamma distributed rates\n");
    fprintf(outfile,
    " Coefficient of variation of rates = %f  (alpha = %f)\n", cv, alpha);
  }
  if (rcategs > 1) {
    fprintf(outfile, "\nState in HMM    Rate of change    Probability\n\n");
    for (i = 0; i < rcategs; i++)
      if (probcat[i] < 0.0001)
        fprintf(outfile, "%9ld%16.3f%20.6f\n", i+1, rrate[i], probcat[i]);
      else if (probcat[i] < 0.001)
          fprintf(outfile, "%9ld%16.3f%19.5f\n", i+1, rrate[i], probcat[i]);
        else if (probcat[i] < 0.01)
            fprintf(outfile, "%9ld%16.3f%18.4f\n", i+1, rrate[i], probcat[i]);
          else
            fprintf(outfile, "%9ld%16.3f%17.3f\n", i+1, rrate[i], probcat[i]);
    putc('\n', outfile);
    if (auto_) {
      fprintf(outfile,
     "Expected length of a patch of sites having the same rate = %8.3f\n",
             1/lambda);
      putc('\n', outfile);
    }
  }
  if (categs > 1) {
    fprintf(outfile, "\nSite category   Rate of change\n\n");
    for (i = 0; i < categs; i++)
      fprintf(outfile, "%9ld%16.3f\n", i+1, rate[i]);
    fprintf(outfile, "\n\n");
  }
}  /* inittable */


void exmake(double lz, long n)
{
  /* pretabulate tables of exponentials so need not do for each site */
  long i;
  double rat;

  for (i = 0; i < categs; i++) {
    rat = rate[i];
    switch (n) {

    case 1:
      expon1i[i] = exp(rat * xi * lz);
      expon1v[i] = exp(rat * xv * lz);
      break;

    case 2:
      expon2i[i] = exp(rat * xi * lz);
      expon2v[i] = exp(rat * xv * lz);
      break;
    }
  }
}  /* exmake */


void alloc_nvd(long num_sibs, nuview_data *local_nvd)
{
  /* Allocate blocks of memory appropriate for the number of siblings
     a given node has */
  local_nvd->yy     = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->wwzz   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->vvzz   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->vzsumr = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->vzsumy = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->sum    = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->sumr   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->sumy   = (double *) Malloc( num_sibs * sizeof (double));
  local_nvd->xx     = (sitelike *) Malloc( num_sibs * sizeof (sitelike));
}  /* alloc_nvd */


void free_nvd(nuview_data *local_nvd)
{
  /* The natural complement to the alloc version */
  free (local_nvd->yy);
  free (local_nvd->wwzz);
  free (local_nvd->vvzz);
  free (local_nvd->vzsumr);
  free (local_nvd->vzsumy);
  free (local_nvd->sum);
  free (local_nvd->sumr);
  free (local_nvd->sumy);
  free (local_nvd->xx);
}  /* free_nvd */

node *invalid_descendant_view(node *p)
{  /* Traverses all descendants of p looking for uninitialized views.
    * Returns NULL if none exist, otherwise a pointer to the first one found.
    */

  node *q = NULL;
  node *s = NULL;

  if ( p == NULL || p->tip )
    return NULL;

  for ( q = p->next; q != p; q = q->next ) {
    s = invalid_descendant_view(q->back);
    if ( s != NULL )
      return s;
  }

  return s;
} /* invalid_descendant_view */


boolean nuview(node *p)
/* current (modified dnaml) nuview */
/* Returns true if view was updated. False if already initialized */
{
  /* Set true to check all descendants for invalid views.
   * False to stop at first descendant with valid view */
  const boolean NUVIEW_SANITY_CHECK = false;

  long i, j, k, l, num_sibs, sib_index;
  nuview_data *local_nvd;
  node *sib_ptr, *sib_back_ptr;
  sitelike p_xx;
  double lw;
  double correction;
  double maxx;

  
  if ( p == NULL ) return false;
  if ( p->tip ) return false; /* Tips do not need to be initialized */

#ifdef NUVIEW_DEBUG
  /* DEBUG: Update nuviews regardless */
  p->initialized = false;
#endif

  if ( !NUVIEW_SANITY_CHECK && p->initialized )
      return false;

  /* Figure out how many siblings the current node has */
  num_sibs    = count_sibs (p);
  /* Recursive calls, should be called for all children */
  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;

    if ( sib_back_ptr != NULL )
      if ( !sib_back_ptr->tip )
        if ( nuview (sib_back_ptr) )
          p->initialized = false;
  }

  if ( p->initialized )
    return false;

  if ( NUVIEW_SANITY_CHECK ) {
    /* At this point, all views downstream should be initialized.
     * If not, we have a problem. */
    assert( invalid_descendant_view(p) == NULL );
  }

  /* Allocate the structure and blocks therein for variables used in
     this function */
  local_nvd = (nuview_data *) Malloc( sizeof (nuview_data));
  alloc_nvd (num_sibs, local_nvd);

  /* Loop 1: makes assignments to tbl based on some combination of
     what's already in tbl and the children's value of v */
  sib_ptr = p;
  for (sib_index=0; sib_index < num_sibs; sib_index++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    
    if (sib_back_ptr != NULL)
      lw = -fabs(p->tyme - sib_back_ptr->tyme);
    else
      lw = 0.0;

    for (i = 0; i < rcategs; i++)
      for (j = 0; j < categs; j++) {
        tbl[i][j]->ww[sib_index]   = exp(tbl[i][j]->ratxi * lw);
        tbl[i][j]->zz[sib_index]   = exp(tbl[i][j]->ratxv * lw);
        tbl[i][j]->wwzz[sib_index] = tbl[i][j]->ww[sib_index] * tbl[i][j]->zz[sib_index];
        tbl[i][j]->vvzz[sib_index] = (1.0 - tbl[i][j]->ww[sib_index]) *
                                                                tbl[i][j]->zz[sib_index];
      }
  }

  /* Loop 2: */
  for (i = 0; i < endsite; i++) {
    correction = 0; 
    maxx = 0;
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {

      /* Loop 2.1 */
      sib_ptr = p;
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
        sib_ptr         = sib_ptr->next;
        sib_back_ptr    = sib_ptr->back;
        

        local_nvd->wwzz[sib_index] = tbl[j][k]->wwzz[sib_index];
        local_nvd->vvzz[sib_index] = tbl[j][k]->vvzz[sib_index];
        local_nvd->yy[sib_index]   = 1.0 - tbl[j][k]->zz[sib_index];
        if (sib_back_ptr != NULL) {
          memcpy(local_nvd->xx[sib_index],
               sib_back_ptr->x[i][j],
               sizeof(sitelike));
          if ( j == 0)
            correction += sib_back_ptr->underflows[i];
        }
        else {
          local_nvd->xx[sib_index][0] = 1.0;
          local_nvd->xx[sib_index][(long)C - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)G - (long)A] = 1.0;
          local_nvd->xx[sib_index][(long)T - (long)A] = 1.0;
        }
      }

      /* Loop 2.2 */
      for (sib_index=0; sib_index < num_sibs; sib_index++) {
        local_nvd->sum[sib_index] =
          local_nvd->yy[sib_index] *
          (freqa * local_nvd->xx[sib_index][(long)A] +
           freqc * local_nvd->xx[sib_index][(long)C] +
           freqg * local_nvd->xx[sib_index][(long)G] +
           freqt * local_nvd->xx[sib_index][(long)T]);
        local_nvd->sumr[sib_index] =
          freqar * local_nvd->xx[sib_index][(long)A] +
          freqgr * local_nvd->xx[sib_index][(long)G];
        local_nvd->sumy[sib_index] =
          freqcy * local_nvd->xx[sib_index][(long)C] +
          freqty * local_nvd->xx[sib_index][(long)T];
        local_nvd->vzsumr[sib_index] =
          local_nvd->vvzz[sib_index] * local_nvd->sumr[sib_index];
        local_nvd->vzsumy[sib_index] =
          local_nvd->vvzz[sib_index] * local_nvd->sumy[sib_index];
      }

      /* Initialize to one, multiply incremental values for every
         sibling a node has */
      p_xx[(long)A] = 1 ;
      p_xx[(long)C] = 1 ; 
      p_xx[(long)G] = 1 ;
      p_xx[(long)T] = 1 ;

      for (sib_index=0; sib_index < num_sibs; sib_index++) {
        p_xx[(long)A] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)A] +
          local_nvd->vzsumr[sib_index];
        p_xx[(long)C] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)C] +
          local_nvd->vzsumy[sib_index];
        p_xx[(long)G] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)G] +
          local_nvd->vzsumr[sib_index];
        p_xx[(long)T] *=
          local_nvd->sum[sib_index] +
          local_nvd->wwzz[sib_index] *
          local_nvd->xx[sib_index][(long)T] +
          local_nvd->vzsumy[sib_index];
      }

      for ( l = 0 ; l < ((long)T - (long)A + 1 ) ; l++ ) {
        if (  p_xx[l] > maxx )
          maxx = p_xx[l];
      }

      /* And the final point of this whole function: */
      memcpy(p->x[i][j], p_xx, sizeof(sitelike));
    }
    p->underflows[i] = 0;
    if ( maxx < MIN_DOUBLE)
      fix_x(p, i, maxx,rcategs);
    p->underflows[i] += correction;
  }

  p->initialized = true;

  free_nvd (local_nvd);
  free (local_nvd);

  return true;
}  /* nuview */


double evaluate(node *p)
{ /* Evaluate and return the log likelihood of the current tree
   * as seen from the branch from p to p->back. If p is the root node,
   * the first child branch is used instead. Views are updated as needed. */
  contribarr tterm;
  static contribarr like, nulike, clai;
  double sum, sum2, sumc=0, y, lz, y1, z1zz, z1yy, 
                prod12, prod1, prod2, prod3, sumterm, lterm;
  long i, j, k, lai;
  node *q, *r;
  sitelike x1, x2;
  sum = 0.0;

  if ( p == curtree.root ) {
    /* Assure that tymes are set alike for all nodes in fork. */
    setnodetymes(p, p->tyme);
    p = p->next;
  }

  r = p;
  q = p->back;
  nuview(r);
  nuview(q);
  y = fabs(r->tyme - q->tyme);
  lz = -y;

  for (i = 0; i < rcategs; i++)
    for (j = 0; j < categs; j++) {
    tbl[i][j]->orig_zz = exp(tbl[i][j]->ratxi * lz);
    tbl[i][j]->z1 = exp(tbl[i][j]->ratxv * lz);
    tbl[i][j]->z1zz = tbl[i][j]->z1 * tbl[i][j]->orig_zz;
    tbl[i][j]->z1yy = tbl[i][j]->z1 - tbl[i][j]->z1zz;
  }
  for (i = 0; i < endsite; i++) {
    k = category[alias[i]-1] - 1;
    for (j = 0; j < rcategs; j++) {
      if (y > 0.0) {
        y1 = 1.0 - tbl[j][k]->z1;
        z1zz = tbl[j][k]->z1zz;
        z1yy = tbl[j][k]->z1yy;
      } else {
        y1 = 0.0;
        z1zz = 1.0;
        z1yy = 0.0;
      }
      memcpy(x1, r->x[i][j], sizeof(sitelike));
      prod1 = freqa * x1[0] + freqc * x1[(long)C - (long)A] +
             freqg * x1[(long)G - (long)A] + freqt * x1[(long)T - (long)A];
      memcpy(x2, q->x[i][j], sizeof(sitelike));
      prod2 = freqa * x2[0] + freqc * x2[(long)C - (long)A] +
             freqg * x2[(long)G - (long)A] + freqt * x2[(long)T - (long)A];
      prod3 = (x1[0] * freqa + x1[(long)G - (long)A] * freqg) *
              (x2[0] * freqar + x2[(long)G - (long)A] * freqgr) +
          (x1[(long)C - (long)A] * freqc + x1[(long)T - (long)A] * freqt) *
         (x2[(long)C - (long)A] * freqcy + x2[(long)T - (long)A] * freqty);
      prod12 = freqa * x1[0] * x2[0] +
               freqc * x1[(long)C - (long)A] * x2[(long)C - (long)A] +
               freqg * x1[(long)G - (long)A] * x2[(long)G - (long)A] +
               freqt * x1[(long)T - (long)A] * x2[(long)T - (long)A];
      tterm[j] = z1zz * prod12 + z1yy * prod3 + y1 * prod1 * prod2;
    }
    sumterm = 0.0;
    for (j = 0; j < rcategs; j++)
      sumterm += probcat[j] * tterm[j];
    lterm = log(sumterm) + p->underflows[i] + q->underflows[i];
    for (j = 0; j < rcategs; j++)
      clai[j] = tterm[j] / sumterm;
    memcpy(contribution[i], clai, sizeof(contribarr));
    if (!auto_ && usertree && (which <= shimotrees))
      l0gf[which - 1][i] = lterm;
    sum += aliasweight[i] * lterm;
  }
  if (auto_) {
    for (j = 0; j < rcategs; j++)
      like[j] = 1.0;
    for (i = 0; i < sites; i++) {
      sumc = 0.0;
      for (k = 0; k < rcategs; k++)
        sumc += probcat[k] * like[k];
      sumc *= lambda;
      if ((ally[i] > 0) && (location[ally[i]-1] > 0)) {
        lai = location[ally[i] - 1];
        memcpy(clai, contribution[lai - 1], sizeof(contribarr));
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc) * clai[j];
      } else {
        for (j = 0; j < rcategs; j++)
          nulike[j] = ((1.0 - lambda) * like[j] + sumc);
      }
      memcpy(like, nulike, sizeof(contribarr));
    }
    sum2 = 0.0;
    for (i = 0; i < rcategs; i++)
      sum2 += probcat[i] * like[i];
    sum += log(sum2);
  }
  curtree.likelihood = sum;
  if (auto_ || !usertree)
    return sum;
  if(which <= shimotrees)
    l0gl[which - 1] = sum;
  if (which == 1) {
    maxwhich = 1;
    maxlogl = sum;
    return sum;
  }
  if (sum > maxlogl) {
    maxwhich = which;
    maxlogl = sum;
  }
  return sum;
}  /* evaluate */


node *pnode(tree *t, node *p) {
  /* Get the "parent nodelet" of p's node group */
  return t->nodep[p->index - 1];
}

double min_child_tyme(node *p) {
  /* Return the minimum tyme of all children. p must be a parent nodelet */
  double min;
  node *q;

  min = 1.0; /* Tymes are always nonpositive */
  
  for ( q = p->next; q != p; q = q->next ) {
    if ( q->back == NULL ) continue;
    if ( q->back->tyme < min )
      min = q->back->tyme;
  }
  
  return min;
} /* min_child_tyme */


double parent_tyme(node *p) {
  /* Return the tyme of the parent of node p. p must be a parent nodelet. */
  if ( p->back ) {
    return p->back->tyme;
  } else {
    return MIN_ROOT_TYME;
  }
} /* parent_tyme */

boolean all_tymes_valid(node *p, double minlength, boolean fix) {
  /* Ensures that all node tymes at node p and descending from it are
   * valid, with all branches being not less than minlength. If any
   * inconsistencies are found, returns true. If 'fix' is given,
   * adjustments are made to make the subtree consistent. Otherwise if
   * assertions are enabled, all inconsistencies are fatal. No effort is
   * made to check that the parent node tyme p->back->tyme is less than
   * p->tyme. */

  node *q;
  double max_tyme;
  boolean ret = true;;

  /* All tips should have tyme == 0.0 */
  if ( p->tip ) {
    if ( p->tyme == 0.0 )
      return true;
    else { /* this would be very bad. */
      if ( fix ) 
        p->tyme = 0.0;
      else
        assert( p->tyme == 0 );

      return false;
    }
  }      

  for ( q = p->next; q != p; q = q->next ) {
    /* All nodes in ring should have same tyme */
    if ( q->tyme != p->tyme ) {
      if ( fix )
        q->tyme = p->tyme;
      else
        assert( q->tyme == p->tyme );
      ret = false;
    }

    /* All substrees should be OK too */
    if ( all_tymes_valid(q->back, minlength, fix) == false )
      ret = false;
  }
  
  /* Tymes cannot be greater than the minimum child time, less
   * branch length */
  max_tyme = min_child_tyme(p) - minlength;
  if ( p->tyme > max_tyme ) {
    if ( fix )
      setnodetymes(p, max_tyme);
    else
      assert( p->tyme < max_tyme );
    return false;
  }

  return ret;
}

boolean valid_tyme(node *p, double tyme) {
  /* Return true if tyme is a valid tyme to assign to node p. tyme must be
   * finite, not greater than any of p's children, and not less than p's
   * parent. Also, tip nodes can only be assigned 0. Otherwise false is
   * returned. */

  /* p must be the parent nodelet of its node group. */

  assert( p->tip != true || tyme == 0.0 );
  assert( tyme <= min_child_tyme(p) );
  assert( tyme >= parent_tyme(p) );

  return true;
} /* valid_tyme */


double set_tyme_evaluate(node *p, double tyme) {
  /* Change tyme of node p and return likelihood
   * Views should be invalidated and regenerated before calling
   * evaluate() anywhere else in the tree. */
  
  /* node *sib_ptr;
     long num_sibs, i; */
 
  assert( valid_tyme(curtree.nodep[p->index - 1], tyme) );

  setnodetymes(p, tyme);
  return evaluate(p);
} /* set_tyme_evaluate */


void makenewv(node *p)
{
  /* Improve a node tyme using Newton-Raphson
   *
   * Slope and curvature are estimated at the current point and used to
   * interpolate a new point. If the curvature is positive, the next point
   * the estimations are pushed uphill by a fraction of the total range.
   * If any interpolation fails to produce a better result, the result is
   * retracted by a given factor toward the original point and tested again.
   *
   * The function terminates when max_iterations have been performed, or
   * when the likelihood improvement for any iteration is less than epsilon,
   * or when a retraction fails. If the iterations are exhausted,
   * 'smoothed' is set false, indicating that further improvement may be
   * possible by additional calls to makenewv(). Otherwise 'smoothed' is left
   * untouched.
   *
   * Define MAKENEWV_DEBUG to get lots of junk printed to stdout.
   * Each loop prints a character, as follows:
   *   '/' ->   Positive curvature, positive slope, moved +
   *   '\' ->   Positive curvature, negative slope, moved -
   *   ')' ->   Negative curvature, moved +
   *   '(' ->   Negative curvature, moved -
   *   '<' ->   Retracting back by retract_factor
   *   'X' ->   Retraction failed, keeping current point 
   */

  /* Tuning constants */
  const double likelihood_epsilon = epsilon/1000.0;
                                                /* Any change in likelihood
                                                   less than this, and we're
                                                   done. */
  const double tyme_delta = epsilon;            /* Small tyme difference used
                                                   to compute slope and
                                                   curvature. */
  const double min_tyme_delta = tyme_delta / 10.0;
                                                /* Absolute smallest
                                                 * tyme_delta */
  const double uphill_step_factor = 0.05;       /* Fraction of the current
                                                   branch length to move uphill
                                                   in positive curvature
                                                   regions. */
  const double retract_factor = 0.5;            /* This defines how far back we
                                                   go if the interpolated point
                                                   is lower than the original
                                                   */
  const double min_tdelta = epsilon;            /* Minimum to which we will
                                                   retract before giving up.
                                                   */
  const long max_iterations = 100;              /* Maximum iterations -
                                                   typically we stop much
                                                   sooner */
  
  double min_tyme, max_tyme;
  double current_tyme, new_tyme;
  double current_likelihood, new_likelihood;
  double x[3], lnl[3];                          /* tyme (x) and log likelihood
                                                   (lnl) points below, at, and
                                                   above the current tyme */
  double s21, s10, slope;
  double curv;
  double uphill_step;
  double tdelta;                /* interpolated point minus current point */
  long iteration;
  boolean done;
 
  
  node *s = curtree.nodep[p->index - 1];
 
#ifdef MAKENEWV_DEBUG
  double start_tyme = s->tyme;
  double start_likelihood = curtree.likelihood;
  long uphill_steps = 0;
#endif /* MAKENEWV_DEBUG */
  assert( valid_tyme(s, s->tyme) );
  
  /* Tyme cannot be less than parent */
  if (s == curtree.root)
    min_tyme = MIN_ROOT_TYME;
  else
    min_tyme = parent_tyme(s) + MIN_BRANCH_LENGTH;

  /* Tyme cannot be greater than any children */
  max_tyme = min_child_tyme(s) - MIN_BRANCH_LENGTH;

  /* Nothing to do if we can't move */
  if ( max_tyme < min_tyme + 2.0*min_tyme_delta ) {
    done = true;
    return;
  }
 
  current_tyme = s->tyme;
  current_likelihood = evaluate(s);

  uphill_step = (max_tyme - min_tyme) * uphill_step_factor;
  
  done = false;
  for ( iteration = 0; iteration < max_iterations; iteration++) {
    /* Evaluate three points for interpolation */
    x[0] = current_tyme - tyme_delta;
    if ( x[0] < min_tyme )
      x[0] = min_tyme;
    x[2] = current_tyme + tyme_delta;
    if ( x[2] > max_tyme )
      x[2] = max_tyme;
    x[1] = (x[0] + x[2]) / 2.0;

    lnl[0] = set_tyme_evaluate(s, x[0]);
    lnl[1] = set_tyme_evaluate(s, x[1]);
    lnl[2] = set_tyme_evaluate(s, x[2]);
   
    /* Compute slopes */
    s21 = (lnl[2] - lnl[1]) / (x[2] - x[1]);
    s10 = (lnl[1] - lnl[0]) / (x[1] - x[0]);
    slope = s21 + s10 / 2.0;

    /* Compute curvature */
    if (fabs(x[2] - x[0]) > epsilon)
      curv = (s21 - s10) / ((x[2] - x[0]) / 2);
    else
      curv = 0.0;
    
    if (curv >= 0.0) {
      /* In negative curvature regions, just move uphill by a
       * fraction of the current length */
      tdelta = copysign(uphill_step, slope);
#ifdef MAKENEWV_DEBUG
      uphill_steps++;
      if ( tdelta > 0 ) putchar('/');
      else putchar('\\');
#endif /* MAKENEWV_DEBUG */
    } else {
      /* Otherwise guess where slope is 0 */
      tdelta = -(slope / curv);
#ifdef MAKENEWV_DEBUG
      if ( tdelta > 0 ) putchar(')');
      else putchar('(');
#endif /* MAKENEWV_DEBUG */
    }
  
    new_tyme = current_tyme + tdelta;
    if ( new_tyme <= min_tyme ) {
      new_tyme = min_tyme;
      tdelta = new_tyme - current_tyme;
    } else if ( new_tyme >= max_tyme ) {
      new_tyme = max_tyme;
      tdelta = new_tyme - current_tyme;
    }
    
    new_likelihood = set_tyme_evaluate(s, new_tyme);
    
    while ( new_likelihood < current_likelihood ) {
      /* If our estimate is worse, retract until we find a better one */
#ifdef MAKENEWV_DEBUG
      putchar('<');
#endif /* MAKENEWV_DEBUG */
      tdelta *= retract_factor;
      uphill_step *= retract_factor;
      
      if ( fabs(tdelta) < min_tdelta ) {
        /* Keep the current point and quit */
        new_likelihood = set_tyme_evaluate(s, current_tyme);
        done = true;
#ifdef MAKENEWV_DEBUG
        putchar('X');
#endif /* MAKENEWV_DEBUG */
        break;
      }

      new_tyme = current_tyme + tdelta;
      new_likelihood = set_tyme_evaluate(s, new_tyme);
      
    }
    
    
    if ( new_likelihood - current_likelihood < likelihood_epsilon ) {
      done = true;
    }
   
    current_likelihood = new_likelihood;
    if ( done ) break;
  }
  
  /* invalidate views */
  invalidate_tyme(s);
  
  if ( !done ) smoothed = false;
  
#ifdef MAKENEWV_DEBUG
  fprintf(stdout, "\nmakenewv(): node %ld: %ld iterations (%f,%f) => (%f,%f)\n",
      p->index, iteration+1, start_tyme, start_likelihood, current_tyme, current_likelihood);
#endif
}  /* makenewv */


void update(node *p)
{
  /* Conditionally optimize tyme at a node */

  if (p == NULL)
    return;

  if ( (!usertree) || (usertree && !lngths) ) {
    makenewv(p);
    return;
  }
}  /* update */


void smooth(node *p)
{
  node *sib_ptr;
  long i, num_sibs;

  if (p == NULL)
    return;
  if (p->tip)
    return;

  /* optimize tyme here */
  update(p);

  smoothed = false;
  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0; i < num_sibs; i++) {
    sib_ptr = sib_ptr->next;
    if (polishing || (smoothit && !smoothed)) {
      /* smooth subtrees */
      smooth(sib_ptr->back);
      p->initialized = false;
      sib_ptr->initialized = false;
    }
    /* optimize tyme again after each subtree */
    update(p);  
  }
}  /* smooth */


void restoradd(node *below, node *newtip, node *newfork, double prevtyme)
{
/* restore "new" tip and fork to place "below".  restore tymes */
/* assumes bifurcation */
  hookup(newfork, below->back);
  hookup(newfork->next, below);
  hookup(newtip, newfork->next->next);
  curtree.nodep[newfork->index-1] = newfork;
  setnodetymes(newfork,prevtyme);
} /* restoradd */


void dnamlk_add(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its descendant, newtip, into the tree. */
  long i;
  boolean done;
  node *p;
  double newtyme;

  /* Get parent nodelets */
  below = curtree.nodep[below->index - 1];
  newfork = curtree.nodep[newfork->index-1];
  newtip = curtree.nodep[newtip->index-1];
  /* Join above node to newfork */
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  /* Join below to newfork->next->next */
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  /* Join newtip to newfork->next */
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  /* assign newfork minimum child tyme */
  if (newtip->tyme < below->tyme)
    p = newtip;
  else
    p = below;
  newtyme = p->tyme;
  setnodetymes(newfork,newtyme);
  /* Move root if inserting there */
  if (curtree.root == below)
    curtree.root = newfork;
  /* If not at root, set newfork tyme to average below/above */
  if (newfork->back != NULL) {
    if (p->tyme > newfork->back->tyme)
      newtyme = (p->tyme + newfork->back->tyme) / 2.0;
    else
      newtyme = p->tyme - INSERT_MIN_TYME;
    setnodetymes(newfork, newtyme);
    /* Now move from p to root, setting parent tymes older than children 
     * by at least INSERT_MIN_TYME */
    do {
      p = curtree.nodep[p->back->index - 1];
      done = (p == curtree.root);
      if (!done)
        done = (curtree.nodep[p->back->index - 1]->tyme < p->tyme - INSERT_MIN_TYME);
      if (!done) {
        setnodetymes(curtree.nodep[p->back->index - 1], p->tyme - INSERT_MIN_TYME);
      }
    } while (!done);
  } else { /* root == newfork */
    /* make root 2x older */
    setnodetymes(newfork, newfork->tyme - 2*INSERT_MIN_TYME);
  }

  /* Invalidate views */
  inittrav(newtip);
  inittrav(newtip->back);
  
  /* Adjust branch lengths throughout */
  for ( i = 1; i < smoothings; i++ ) {
    smoothed = true;
    smooth(newfork);
    smooth(newfork->back);
    if ( smoothed ) break;
  }
}  /* dnamlk_add */


void dnamlk_re_move(node **item, node **fork, boolean tempadd)
{
  /* removes nodes item and its ancestor, fork, from the tree.
    the new descendant of fork's ancestor is made to be
    fork's second descendant (other than item).  Also
    returns pointers to the deleted nodes, item and fork */
  node *p, *q;
  long i;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *item = curtree.nodep[(*item)->index-1];
  *fork = curtree.nodep[(*item)->back->index - 1];
  if (curtree.root == *fork) {
    if (*item == (*fork)->next->back)
      curtree.root = (*fork)->next->next->back;
    else
      curtree.root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
  inittrav(p);
  inittrav(q);
  if (tempadd)
    return;

  for ( i = 1; i <= smoothings; i++ ) {
    smoothed = true;
    smooth(q);
    if ( smoothit )
      smooth(q->back);
    if ( smoothed ) break;
  }
}  /* dnamlk_re_move */


void tryadd(node *p, node **item, node **nufork)
{  /* temporarily adds one fork and one tip to the tree.
    if the location where they are added yields greater
    likelihood than other locations tested up to that
    time, then keeps that location as there */

  dnamlk_add(p, *item, *nufork);
  like = evaluate(p);
  if (lastsp) {
      if (like >= bestree.likelihood || bestree.likelihood == UNDEFINED)  {
            copy_(&curtree, &bestree, nonodes, rcategs);
            if (global2) 
                save_tree_tyme(&curtree,tymes);
      }
  }
  if (like > bestyet || bestyet == UNDEFINED) {
    bestyet = like;
    there = p;
  }
  dnamlk_re_move(item, nufork, true);
  if ( global2 ) {
      restore_saved_tyme(&curtree,tymes);
  }
}  /* tryadd */


void addpreorder(node *p, node *item_, node *nufork_, boolean contin,
                        boolean continagain)
{
  /* traverses a binary tree, calling function tryadd
    at a node before calling tryadd at its descendants */
  node *item, *nufork;

  item = item_;
  nufork = nufork_;
  if (p == NULL)
    return;
  tryadd(p, &item, &nufork);
  contin = continagain;
  if ((!p->tip) && contin) {
    /* assumes bifurcation (OK) */
    addpreorder(p->next->back, item, nufork, contin, continagain);
    addpreorder(p->next->next->back, item, nufork, contin, continagain);
  }
}  /* addpreorder */


boolean tryrearr(node *p)
{
  /* evaluates one rearrangement of the tree.
    if the new tree has greater likelihood than the old
    keeps the new tree and returns true.
    otherwise, restores the old tree and returns false. */

  node *frombelow, *whereto, *forknode;
  double oldlike, prevtyme, like_delta;
  boolean wasonleft;

  if (p == curtree.root)
    return false;
  forknode = curtree.nodep[p->back->index - 1];
  if (forknode == curtree.root)
    return false;
  oldlike = bestyet;
  prevtyme = forknode->tyme;
  /* assumes bifurcation (OK) */
  if (forknode->next->back == p) {
    frombelow = forknode->next->next->back;
    wasonleft = true;
  }
  else {
    frombelow = forknode->next->back;
    wasonleft = false;
  }
  whereto = curtree.nodep[forknode->back->index - 1];
  dnamlk_re_move(&p, &forknode, true);
  dnamlk_add(whereto, p, forknode);
  like = evaluate(p);
  like_delta = like - oldlike;
  if ( like_delta < LIKE_EPSILON && oldlike != UNDEFINED) {
    dnamlk_re_move(&p, &forknode, true);
    restoradd(frombelow, p, forknode, prevtyme);
    if (wasonleft && (forknode->next->next->back == p)) {
       hookup (forknode->next->back, forknode->next->next);
       hookup (forknode->next, p);
    }
    curtree.likelihood = oldlike;

    /* assumes bifurcation (OK) */
    inittrav(forknode);
    inittrav(forknode->next);
    inittrav(forknode->next->next);

    return false;
  } else {
    bestyet = like;
  }
  return true;
}  /* tryrearr */


boolean repreorder(node *p)
{
  /* traverses a binary tree, calling function tryrearr
    at a node before calling tryrearr at its descendants.
    Returns true the first time rearrangement increases the
    tree's likelihood. */

  if (p == NULL)
    return false;
  if ( !tryrearr(p) )
    return false;
  if (p->tip)
    return true;
  /* assumes bifurcation */
  if ( !repreorder(p->next->back) )
    return false;
  if ( !repreorder(p->next->next->back) )
    return false;
  
  return true;
}  /* repreorder */


void rearrange(node **r)
{
  /* traverses the tree (preorder), finding any local
    rearrangement which increases the likelihood.
    if traversal succeeds in increasing the tree's
    likelihood, function rearrange runs traversal again */
  while ( repreorder(*r) )
    /* continue */;

}  /* rearrange */


void initdnamlnode(node **p, node **grbg, node *q, long len, long nodei,
                     long *ntips, long *parens, initops whichinit,
                     pointarray treenode, pointarray nodep, Char *str, Char *ch,
                     FILE *intree)
{
  /* Initializes each node as it is read from user tree by treeread().
   * whichinit specifies which type of initialization is to be done.
   */
  boolean minusread;
  double valyew, divisor;

  switch (whichinit) {
  case bottom:
    gnu(grbg, p);
    (*p)->index = nodei;
    (*p)->tip = false;
    malloc_pheno((*p), endsite, rcategs);
    nodep[(*p)->index - 1] = (*p);
    break;
  case nonbottom:
    gnu(grbg, p);
    malloc_pheno(*p, endsite, rcategs);
    (*p)->index = nodei;
    break;
  case tip:
    match_names_to_data (str, nodep, p, spp);
    break;
  case iter:
    (*p)->initialized = false;
    (*p)->v = initialv;
    (*p)->iter = true;
    if ((*p)->back != NULL)
      (*p)->back->iter = true;
    break;
  case length:
    processlength(&valyew, &divisor, ch, &minusread, intree, parens);
    (*p)->v = valyew / divisor / fracchange;
    (*p)->iter = false;
    if ((*p)->back != NULL) {
      (*p)->back->v = (*p)->v;
      (*p)->back->iter = false;
    }
    break;
  case hslength:
    break;
  case hsnolength:
    if (usertree && lengthsopt && lngths) {
      printf("Warning: one or more lengths not defined in user tree number %ld.\n", which);
      printf("DNAMLK will attempt to optimize all branch lengths.\n\n");
      lngths = false;
    }
    break;
  case treewt:
    break;
  case unittrwt:
    break;
  }
} /* initdnamlnode */


void tymetrav(node *p, double *x)
{
  /* set up times of nodes */
  node *sib_ptr;
  long i, num_sibs;
  double xmax;

  xmax = 0.0;
  if (!p->tip) {
    sib_ptr  = p;
    num_sibs = count_sibs(p);
    for (i=0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      tymetrav(sib_ptr->back, x);
      if (xmax > (*x))
        xmax = (*x);
    }
  } else
    (*x)     = 0.0;
  setnodetymes(p,xmax);
  (*x) = p->tyme - p->v;
}  /* tymetrav */


void dnamlk_coordinates(node *p, long *tipy)
{
  /* establishes coordinates of nodes */
  node *q, *first, *last, *pp1 =NULL, *pp2 =NULL;
  long num_sibs, p1, p2, i;

  if (p->tip) {
    p->xcoord = 0;
    p->ycoord = (*tipy);
    p->ymin   = (*tipy);
    p->ymax   = (*tipy);
    (*tipy)  += down;
    return;
  }
  q = p->next;
  do {
    dnamlk_coordinates(q->back, tipy);
    q = q->next;
  } while (p != q);
  num_sibs = count_sibs(p);
  p1 = (long)((num_sibs+1)/2.0);
  p2 = (long)((num_sibs+2)/2.0);
  i = 1;
  q = p->next;
  first  = q->back;
  do {
    if (i == p1) pp1 = q->back;
    if (i == p2) pp2 = q->back;
    last = q->back;
    q = q->next;
    i++;
  } while (q != p);
  p->xcoord = (long)(0.5 - over * p->tyme);
  p->ycoord = (pp1->ycoord + pp2->ycoord) / 2;
  p->ymin = first->ymin;
  p->ymax = last->ymax;
}  /* dnamlk_coordinates */


void dnamlk_drawline(long i, double scale)
{
  /* draws one row of the tree diagram by moving up tree */
  node *p, *q, *r, *first =NULL, *last =NULL;
  long n, j;
  boolean extra, done;

  p = curtree.root;
  q = curtree.root;
  extra = false;
  if ((long)(p->ycoord) == i) {
    if (p->index - spp >= 10)
      fprintf(outfile, "-%2ld", p->index - spp);
    else
      fprintf(outfile, "--%ld", p->index - spp);
    extra = true;
  } else
    fprintf(outfile, "  ");
  do {
    if (!p->tip) {
      r = p->next;
      done = false;
      do {
        if (i >= r->back->ymin && i <= r->back->ymax) {
          q = r->back;
          done = true;
        }
        r = r->next;
      } while (!(done || r == p));
      first = p->next->back;
      r = p->next;
      while (r->next != p)
        r = r->next;
      last = r->back;
    }
    done = (p == q);
    n = (long)(scale * ((long)(p->xcoord) - (long)(q->xcoord)) + 0.5);
    if (n < 3 && !q->tip)
      n = 3;
    if (extra) {
      n--;
      extra = false;
    }
    if ((long)(q->ycoord) == i && !done) {
      if (p->ycoord != q->ycoord)
        putc('+', outfile);
      else
        putc('-', outfile);
      if (!q->tip) {
        for (j = 1; j <= n - 2; j++)
          putc('-', outfile);
        if (q->index - spp >= 10)
          fprintf(outfile, "%2ld", q->index - spp);
        else
          fprintf(outfile, "-%ld", q->index - spp);
        extra = true;
      } else {
        for (j = 1; j < n; j++)
          putc('-', outfile);
      }
    } else if (!p->tip) {
      if ((long)(last->ycoord) > i && (long)(first->ycoord) < i && 
           i != (long)(p->ycoord)) {
        putc('!', outfile);
        for (j = 1; j < n; j++)
          putc(' ', outfile);
      } else {
        for (j = 1; j <= n; j++)
          putc(' ', outfile);
      }
    } else {
      for (j = 1; j <= n; j++)
        putc(' ', outfile);
    }
    if (p != q)
      p = q;
  } while (!done);
  if ((long)(p->ycoord) == i && p->tip) {
    for (j = 0; j < nmlngth; j++)
      putc(nayme[p->index - 1][j], outfile);
  }
  putc('\n', outfile);
}  /* dnamlk_drawline */


void dnamlk_printree()
{
 /* prints out diagram of the tree */
  long tipy;
  double scale;
  long i;
  node *p;

  if (!treeprint)
    return;
  putc('\n', outfile);
  tipy = 1;
  dnamlk_coordinates(curtree.root, &tipy);
  p = curtree.root;
  while (!p->tip)
    p = p->next->back;
  scale = 1.0 / (long)(p->tyme - curtree.root->tyme + 1.000);
  putc('\n', outfile);
  for (i = 1; i <= tipy - down; i++)
    dnamlk_drawline(i, scale);
  putc('\n', outfile);
}  /* dnamlk_printree */


void describe(node *p)
{
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;
  double v;

  if (p == curtree.root)
    fprintf(outfile, " root         ");
  else
    fprintf(outfile, "%4ld          ", p->back->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  if (p != curtree.root) {
    fprintf(outfile, "%11.5f", fracchange * (p->tyme - curtree.root->tyme));
    v = fracchange * (p->tyme - curtree.nodep[p->back->index - 1]->tyme);
    fprintf(outfile, "%13.5f", v);
  }
  putc('\n', outfile);
  if (!p->tip) {

    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      describe(sib_back_ptr);
    }
  }
}  /* describe */


void reconstr(node *p, long n)
{
  /* reconstruct and print out base at site n+1 at node p */
  long i, j, k, m, first, second, num_sibs;
  double f, sum, xx[4];
  node *q;

  if ((ally[n] == 0) || (location[ally[n]-1] == 0))
    putc('.', outfile);
  else {
    j = location[ally[n]-1] - 1;
    for (i = 0; i < 4; i++) {
      f = p->x[j][mx-1][i];
      num_sibs = count_sibs(p);
      q = p;
      for (k = 0; k < num_sibs; k++) {
        q = q->next;
        f *= q->x[j][mx-1][i];
      }
      f = sqrt(f);
      xx[i] = f;
    }
    xx[0] *= freqa;
    xx[1] *= freqc;
    xx[2] *= freqg;
    xx[3] *= freqt;
    sum = xx[0]+xx[1]+xx[2]+xx[3];
    for (i = 0; i < 4; i++)
      xx[i] /= sum;
    first = 0;
    for (i = 1; i < 4; i++)
      if (xx [i] > xx[first])
        first = i;
    if (first == 0)
      second = 1;
    else
      second = 0;
    for (i = 0; i < 4; i++)
      if ((i != first) && (xx[i] > xx[second]))
        second = i;
    m = 1 << first;
    if (xx[first] < 0.4999995)
      m = m + (1 << second);
    if (xx[first] > 0.95)
      putc(toupper(basechar[m - 1]), outfile);
    else
      putc(basechar[m - 1], outfile);
    if (rctgry && rcategs > 1)
      mx = mp[n][mx - 1];    
    else
      mx = 1;
  }
} /* reconstr */


void rectrav(node *p, long m, long n)
{
  /* print out segment of reconstructed sequence for one branch */
  long num_sibs, i;
  node *sib_ptr;

  putc(' ', outfile);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index-1][i], outfile);
  } else
    fprintf(outfile, "%4ld      ", p->index - spp);
  fprintf(outfile, "  ");
  mx = mx0;
  for (i = m; i <= n; i++) {
    if ((i % 10 == 0) && (i != m))
      putc(' ', outfile);
    if (p->tip)
      putc(y[p->index-1][i], outfile);
    else
      reconstr(p, i);
  }
  putc('\n', outfile);
  if (!p->tip) {
    num_sibs = count_sibs(p);
    sib_ptr = p;
    for (i = 0; i < num_sibs; i++) {
      sib_ptr = sib_ptr->next;
      rectrav(sib_ptr->back, m, n);
    }
  }
  mx1 = mx;
}  /* rectrav */


void summarize()
{
  long i, j, mm;
  double mode, sum;
  double like[maxcategs], nulike[maxcategs];
  double **marginal;

  mp = (long **)Malloc(sites * sizeof(long *));
  for (i = 0; i <= sites-1; ++i)
    mp[i] = (long *)Malloc(sizeof(long)*rcategs);
  fprintf(outfile, "\nLn Likelihood = %11.5f\n\n", curtree.likelihood);
  fprintf(outfile, " Ancestor      Node      Node Height     Length\n");
  fprintf(outfile, " --------      ----      ---- ------     ------\n");
  describe(curtree.root);
  putc('\n', outfile);
  if (rctgry && rcategs > 1) {
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        mp[i][j] = j + 1;
        for (k = 1; k <= rcategs; k++) {
          if (k != j + 1) {
            if (lambda * probcat[k - 1] * like[k - 1] > nulike[j]) {
              nulike[j] = lambda * probcat[k - 1] * like[k - 1];
              mp[i][j] = k;
            }
          }
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    mode = 0.0;
    mx = 1;
    for (i = 1; i <= rcategs; i++) {
      if (probcat[i - 1] * like[i - 1] > mode) {
        mx = i;
        mode = probcat[i - 1] * like[i - 1];
      }
    }
    mx0 = mx;
    fprintf(outfile,
 "Combination of categories that contributes the most to the likelihood:\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 1; i <= sites; i++) {
      fprintf(outfile, "%ld", mx);
      if (i % 10 == 0)
        putc(' ', outfile);
      if (i % 60 == 0 && i != sites) {
        putc('\n', outfile);
        for (j = 1; j <= nmlngth + 3; j++)
          putc(' ', outfile);
      }
      mx = mp[i - 1][mx - 1];
    }
    fprintf(outfile, "\n\n");
    marginal = (double **) Malloc( sites*sizeof(double *));
    for (i = 0; i < sites; i++)
      marginal[i] = (double *) Malloc( rcategs*sizeof(double));
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = sites - 1; i >= 0; i--) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++) {
          if (k != j + 1)
              nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        if ((ally[i] > 0) && (location[ally[i]-1] > 0))
          nulike[j] *= contribution[location[ally[i] - 1] - 1][j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++) {
        nulike[j] /= sum;
        marginal[i][j] = nulike[j];
      }
      memcpy(like, nulike, rcategs * sizeof(double));
    }
    for (i = 0; i < rcategs; i++)
      like[i] = 1.0;
    for (i = 0; i < sites; i++) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++) {
        nulike[j] = (lambda1 + lambda * probcat[j]) * like[j];
        for (k = 1; k <= rcategs; k++) {
          if (k != j + 1)
              nulike[j] += lambda * probcat[k - 1] * like[k - 1];
        }
        marginal[i][j] *= like[j] * probcat[j];
        sum += nulike[j];
      }
      for (j = 0; j < rcategs; j++)
        nulike[j] /= sum;
      memcpy(like, nulike, rcategs * sizeof(double));
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        sum += marginal[i][j];
      for (j = 0; j < rcategs; j++)
        marginal[i][j] /= sum;
    }
    fprintf(outfile, "Most probable category at each site if > 0.95 probability (\".\" otherwise)\n\n");
    for (i = 1; i <= nmlngth + 3; i++)
      putc(' ', outfile);
    for (i = 0; i < sites; i++) {
      sum = 0.0;
      for (j = 0; j < rcategs; j++)
        if (marginal[i][j] > sum) {
          sum = marginal[i][j];
          mm = j;
        }
        if (sum >= 0.95)
        fprintf(outfile, "%ld", mm+1);
      else
        putc('.', outfile);
      if ((i+1) % 60 == 0) {
        if (i != 0) {
          putc('\n', outfile);
          for (j = 1; j <= nmlngth + 3; j++)
            putc(' ', outfile);
        }
      }
      else if ((i+1) % 10 == 0)
        putc(' ', outfile);
    }
    putc('\n', outfile);
    for (i = 0; i < sites; i++)
      free(marginal[i]);
    free(marginal);
  }
  putc('\n', outfile);
  putc('\n', outfile);
  if (hypstate) {
    fprintf(outfile, "Probable sequences at interior nodes:\n\n");
    fprintf(outfile, "  node       ");
    for (i = 0; (i < 13) && (i < ((sites + (sites-1)/10 - 39) / 2)); i++)
      putc(' ', outfile);
    fprintf(outfile, "Reconstructed sequence (caps if > 0.95)\n\n");
    if (!rctgry || (rcategs == 1))
      mx0 = 1;
    for (i = 0; i < sites; i += 60) {
      k = i + 59;
      if (k >= sites)
        k = sites - 1;
      rectrav(curtree.root, i, k);
      putc('\n', outfile);
      mx0 = mx1;
    }
  }
  for (i = 0; i < sites; ++i)
    free(mp[i]);
  free(mp);
}  /* summarize */


void dnamlk_treeout(node *p)
{
  /* write out file with representation of final tree */
  node *sib_ptr;
  long i, n, w, num_sibs;
  Char c;
  double x;

  if (p->tip) {
    n = 0;
    for (i = 1; i <= nmlngth; i++) {
      if (nayme[p->index - 1][i - 1] != ' ')
        n = i;
    }
    for (i = 0; i < n; i++) {
      c = nayme[p->index - 1][i];
      if (c == ' ')
        c = '_';
      putc(c, outtree);
    }
    col += n;
  } else {
    sib_ptr = p;
    num_sibs = count_sibs(p);
    putc('(', outtree);
    col++;

    for (i=0; i < (num_sibs - 1); i++) {
      sib_ptr = sib_ptr->next;
      dnamlk_treeout(sib_ptr->back);
      putc(',', outtree);
      col++;
      if (col > 55) {
        putc('\n', outtree);
        col = 0;
      }
    }
    sib_ptr = sib_ptr->next;
    dnamlk_treeout(sib_ptr->back);
    putc(')', outtree);
    col++;
  }
  if (p == curtree.root) {
    fprintf(outtree, ";\n");
    return;
  }
  x = fracchange * (p->tyme - curtree.nodep[p->back->index - 1]->tyme);
  if (x > 0.0)
    w = (long)(0.4342944822 * log(x));
  else if (x == 0.0)
    w = 0;
  else
    w = (long)(0.4342944822 * log(-x)) + 1;
  if (w < 0)
    w = 0;
  fprintf(outtree, ":%*.5f", (int)(w + 7), x);
  col += w + 8;
}  /* dnamlk_treeout */


void nodeinit(node *p)
{
  /* set up times at one node */
  node *sib_ptr, *sib_back_ptr;
  long i, num_sibs;
  double lowertyme;

  sib_ptr = p;
  num_sibs = count_sibs(p);

  /* lowertyme = lowest of children's times */
  lowertyme = p->next->back->tyme;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    if (sib_back_ptr->tyme < lowertyme)
      lowertyme = sib_back_ptr->tyme;
  }

  setnodetymes(p, lowertyme - 0.1);

  sib_ptr = p;
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    sib_back_ptr->v = sib_back_ptr->tyme - p->tyme;
    sib_ptr->v = sib_back_ptr->v;
  }
}  /* nodeinit */

void invalidate_traverse(node *p)
{ /* Invalidates p's view and all views looking toward p from p->back
   * on out. */
  node *q;
  
  if (p == NULL)
    return;
  if (p->tip)
    return;

  p->initialized = false;

  q = p->back;
  if ( q == NULL ) return;
  if ( q->tip ) return;

  /* Call ourselves on p->back's sibs */
  for ( q = q->next ; q != p->back ; q = q->next) {
    invalidate_traverse(q);
  }
} /* invalidate_traverse */

void invalidate_tyme(node *p) {
  /* Must be called on a node after changing its tyme, and before calling
   * evaluate on any other node. */
  
  node *q;

  if ( p == NULL )
    return;
  invalidate_traverse(p);
  if ( p->tip )
    return;
  for ( q = p->next; q != p; q = q->next ) {
    invalidate_traverse(q);
  }
}

void setnodetymes(node* p, double newtyme)
{ /* Set node tyme for an entire fork. Also clears initialized flags on this 
   * fork, but not recursively. inittrav() must be called before evaluating
   * elsewhere. */
  node * q;

  p->tyme = newtyme;
  p->initialized = false;
  if ( p->tip ) return;
  for ( q = p->next; q != p; q = q->next ) {
    if ( !q ) break;
    q->tyme = newtyme;
    q->initialized = false;
  }
} /* setnodetymes */

void initrav(node *p)
{ /* traverse calling nodeinit() on each internal node */

  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up times throughout tree */
  if (p->tip)
    return;

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    initrav(sib_back_ptr);
  }

  nodeinit(p);
}  /* initrav */

void travinit(node *p)
{ /* This function is obsolete. nuview() will call itself on
   * any uninitialized children, and set initialized when
   * it is done. */
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to set up initial values */
  if (p == NULL)
    return;
  if (p->tip)
    return;
  if (p->initialized)
    return;

  sib_ptr = p;
  num_sibs = count_sibs(p);
  for (i=0 ; i < num_sibs; i++) {
    sib_ptr      = sib_ptr->next;
    sib_back_ptr = sib_ptr->back;
    travinit(sib_back_ptr);
  }

  nuview(p);
  p->initialized = true;
}  /* travinit */


void travsp(node *p)
{ /* This function is obsolete. To update all nuviews, just
   * call nuview() from all directions (i.e. from each tip
   * and from the root) */
  long i, num_sibs;
  node *sib_ptr, *sib_back_ptr;

  /* traverse to find tips */
  if (p == curtree.root)
    travinit(p);
  if (p->tip)
    travinit(p->back);
  else {
    sib_ptr = p;
    num_sibs = count_sibs(p);
    for (i=0 ; i < num_sibs; i++) {
      sib_ptr      = sib_ptr->next;
      sib_back_ptr = sib_ptr->back;
      travsp(sib_back_ptr);
    }
  }
}  /* travsp */


void treevaluate()
{
  /* evaluate likelihood of tree, after iterating branch lengths */
  long i, j,  num_sibs;
  node *sib_ptr;

  polishing = true;
  smoothit = true;
  if ( usertree && !lngths ) {
    for (i = 0; i < spp; i++)
        curtree.nodep[i]->initialized = false;
    for (i = spp; i < nonodes; i++) {
        sib_ptr = curtree.nodep[i];
        sib_ptr->initialized = false;
        num_sibs = count_sibs(sib_ptr);
        for (j=0 ; j < num_sibs; j++) {
          sib_ptr      = sib_ptr->next;
          sib_ptr->initialized = false;
        }
    }
    initrav(curtree.root);
    travsp(curtree.root);
  }
  for (i = 1; i <= smoothings * 4; i++)
    smooth(curtree.root);
  evaluate(curtree.root);
}  /* treevaluate */


void maketree()
{
  /* constructs a binary tree from the pointers in curtree.nodep,
     adds each node at location which yields highest likelihood
     then rearranges the tree for greatest likelihood */

  long i, j, numtrees;
  double x;
  node *item, *nufork, *dummy, *q, *root=NULL;
  boolean succeded, dummy_haslengths, dummy_first, goteof;
  long max_nonodes;        /* Maximum number of nodes required to
                            * express all species in a bifurcating tree
                            * */
  long nextnode;
  pointarray dummy_treenode=NULL;
  double oldbest;
  node *tmp;

  inittable();  

  if (!usertree) {
    for (i = 1; i <= spp; i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree.root = curtree.nodep[spp];
    curtree.root->back = NULL;
    for (i = 0; i < spp; i++)
       curtree.nodep[i]->back = NULL;
    for (i = spp; i < nonodes; i++) {
       q = curtree.nodep[i];
       q->back = NULL;
       while ((q = q->next) != curtree.nodep[i])
         q->back = NULL;
    }
    polishing = false;
    dnamlk_add(curtree.nodep[enterorder[0] - 1], curtree.nodep[enterorder[1] - 1],
        curtree.nodep[spp]);
    if (progress) {
      printf("\nAdding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    lastsp = false;
    smoothit = false;
    for (i = 3; i <= spp; i++) {
      bestree.likelihood = UNDEFINED;
      bestyet = UNDEFINED;
      there = curtree.root;
      item = curtree.nodep[enterorder[i - 1] - 1];
      nufork = curtree.nodep[spp + i - 2];
      lastsp = (i == spp);
      addpreorder(curtree.root, item, nufork, true, true);
      dnamlk_add(there, item, nufork);
      like = evaluate(curtree.root);
      copy_(&curtree, &bestree, nonodes, rcategs);
      rearrange(&curtree.root);
      if (curtree.likelihood > bestree.likelihood) {
        copy_(&curtree, &bestree, nonodes, rcategs);
      }
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      if (lastsp && global) {
        if (progress) {
          printf("Doing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= nonodes; j++)
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");
        }
        global2 = true;
        do {
          succeded = false;
          if (progress)
            printf("   ");
          save_tree_tyme(&curtree, tymes);
          for (j = 0; j < nonodes; j++) {
            oldbest = bestree.likelihood;
            bestyet = UNDEFINED;
            item = curtree.nodep[j];
            if (item != curtree.root) {
              nufork = curtree.nodep[curtree.nodep[j]->back->index - 1];
              
              if (nufork != curtree.root) {  
                tmp = nufork->next->back;
                if (tmp == item) 
                    tmp = nufork->next->next->back; 
                    /* can't figure out why we never get here */
              }
              else {
                  if (nufork->next->back != item)
                      tmp  = nufork->next->back;
                  else tmp = nufork->next->next->back;
              } /* if we add item at tmp we have done nothing */
              dnamlk_re_move(&item, &nufork, false);
              there = curtree.root;
              addpreorder(curtree.root, item, nufork, true, true);
              if ( tmp != there && bestree.likelihood > oldbest)
                succeded = true;
              dnamlk_add(there, item, nufork);
              restore_saved_tyme(&curtree,tymes);
            }
            if (progress) {
              if ( j % (( nonodes / 72 ) + 1 ) == 0 )
                putchar('.');
              fflush(stdout);
            }
          } 
          if (progress)
            putchar('\n');
        } while ( succeded );
      }
    }
    if (njumble > 1 && lastsp) {
      for (i = 0; i < spp; i++ )
        dnamlk_re_move(&curtree.nodep[i], &dummy, false);
      if (jumb == 1 || bestree2.likelihood < bestree.likelihood)
        copy_(&bestree, &bestree2, nonodes, rcategs);
    }
    if (jumb == njumble) {
      if (njumble > 1)
        copy_(&bestree2, &curtree, nonodes, rcategs);
      else copy_(&bestree, &curtree, nonodes, rcategs);
      fprintf(outfile, "\n\n");
      treevaluate();
      curtree.likelihood = evaluate(curtree.root);
      dnamlk_printree();
      summarize();
      if (trout) {
        col = 0;
        dnamlk_treeout(curtree.root);
      }
    } 
  } else { /* if ( usertree ) */
    openfile(&intree, INTREE, "input tree file", "r", progname, intreename);
    inittable_for_usertree (intree);
    numtrees = countsemic(&intree);
    if(numtrees > MAXSHIMOTREES)
      shimotrees = MAXSHIMOTREES;
    else
      shimotrees = numtrees;
    if (numtrees > 2)
      initseed(&inseed, &inseed0, seed);
    l0gl = (double *)Malloc(shimotrees * sizeof(double));
    l0gf = (double **)Malloc(shimotrees * sizeof(double *));
    for (i=0; i < shimotrees; ++i)
      l0gf[i] = (double *)Malloc(endsite * sizeof(double));
    if (treeprint) {
      fprintf(outfile, "User-defined tree");
      if (numtrees > 1)
        putc('s', outfile);
      fprintf(outfile, ":\n\n");
    }
    fprintf(outfile, "\n\n");
    which = 1;
    max_nonodes = nonodes;
    while (which <= numtrees) {
      
      /* These initializations required each time through the loop
         since multiple trees require re-initialization */
      dummy_haslengths = true;
      nextnode         = 0;
      dummy_first      = true;
      goteof           = false;
      lngths           = lengthsopt;
      nonodes          = max_nonodes;

      treeread(intree, &root, dummy_treenode, &goteof, &dummy_first,
               curtree.nodep, &nextnode, &dummy_haslengths, &grbg, 
               initdnamlnode, false, nonodes);

      nonodes = nextnode;
      
      root = curtree.nodep[root->index - 1];
      curtree.root = root;

      if (lngths)
        tymetrav(curtree.root, &x);

      if (goteof && (which <= numtrees)) {
        /* if we hit the end of the file prematurely */
        printf ("\n");
        printf ("ERROR: trees missing at end of file.\n");
        printf ("\tExpected number of trees:\t\t%ld\n", numtrees);
        printf ("\tNumber of trees actually in file:\t%ld.\n\n", which - 1);
        exxit(-1);
      }
      curtree.start = curtree.nodep[0]->back;
      treevaluate();
      dnamlk_printree();
      summarize();
      if (trout) {
        col = 0;
        dnamlk_treeout(curtree.root);
      }
      if(which < numtrees){
        freex_notip(nonodes, curtree.nodep);
        gdispose(curtree.root, &grbg, curtree.nodep);
      }
      which++;
    }      

    FClose(intree);
    if (!auto_ && numtrees > 1 && weightsum > 1 )
      standev2(numtrees, maxwhich, 0, endsite, maxlogl, l0gl, l0gf,
               aliasweight, seed);
  }
  
  if (jumb == njumble) {
    if (progress) {
      printf("\nOutput written to file \"%s\"\n\n", outfilename);
      if (trout)
        printf("Tree also written onto file \"%s\"\n\n", outtreename);
    }
    free(contribution);
    freex(nonodes, curtree.nodep);
    if (!usertree) {
      freex(nonodes, bestree.nodep);
      if (njumble > 1)
        freex(nonodes, bestree2.nodep);
    }
  }
  free(root);
  freetable();
} /* maketree */

/*?? Dnaml has a clean-up function for freeing memory, closing files, etc.
     Put one here too? */

int main(int argc, Char *argv[])
{  /* DNA Maximum Likelihood with molecular clock */

#ifdef MAC
  argc = 1;                /* macsetup("Dnamlk", "Dnamlk");        */
  argv[0] = "Dnamlk";
#endif
  init(argc,argv);
  progname = argv[0];
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  openfile(&outfile, OUTFILE, "output file", "w", argv[0], outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  datasets = 1;
  mulsets = false;
  firstset = true;
  
  doinit();

  ttratio0 = ttratio;

  /* Open output tree, categories, and weights files if needed */
  if ( trout )
    openfile(&outtree, OUTTREE, "output tree file", "w", argv[0], outtreename);
  if ( ctgry )
    openfile(&catfile, CATFILE, "categories file", "r", argv[0], catfilename);
  if ( weights || justwts )
   openfile(&weightfile, WEIGHTFILE, "weights file", "r", argv[0], weightfilename);

  /* Data set loop */
  for ( ith = 1; ith <= datasets; ith++ ) {
    ttratio = ttratio0;
    if ( datasets > 1 ) {
      fprintf(outfile, "Data set # %ld:\n\n", ith);
      if ( progress )
        printf("\nData set # %ld:\n", ith);
    }
    getinput();
    if ( ith == 1 )
      firstset = false;
    
    /* Jumble loop */
    for ( jumb = 1; jumb <= njumble; jumb++ )
      maketree();
  }

  /* Close files */
  FClose(infile);  
  FClose(outfile);
  if ( trout )
    FClose(outtree);
  if ( ctgry )
    FClose(catfile);
  if ( weights || justwts )
    FClose(weightfile);

#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}  /* DNA Maximum Likelihood with molecular clock */
