#include "phylip.h"
#include "seq.h"

/* version 3.6. (c) Copyright 1993-2005 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Doug Buxton.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

typedef enum {
  seqs, morphology, restsites, genefreqs
} datatype;

typedef enum {
  dna, rna, protein
} seqtype;


#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   seqboot_inputnumbers(void);
void   seqboot_inputfactors(void);
void   inputoptions(void);
void   seqboot_inputdata(void);
void   allocrest(void);
void   allocnew(void);
void   doinput(int argc, Char *argv[]);
void   bootweights(void);
void   sppermute(long);
void   charpermute(long, long);
void   writedata(void);
void   writeweights(void);
void   writecategories(void);
void   writeauxdata(steptr, FILE*);
void   writefactors(void);
void   bootwrite(void);
void   seqboot_inputaux(steptr, FILE*);
/* function prototypes */
#endif


FILE *outcatfile, *outweightfile, *outmixfile, *outancfile, *outfactfile;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], catfilename[FNMLNGTH], outcatfilename[FNMLNGTH],
  weightfilename[FNMLNGTH], outweightfilename[FNMLNGTH], mixfilename[FNMLNGTH], outmixfilename[FNMLNGTH], ancfilename[FNMLNGTH], outancfilename[FNMLNGTH],
  factfilename[FNMLNGTH], outfactfilename[FNMLNGTH];
long sites, loci, maxalleles, groups, newsites, newersites,
  newgroups, newergroups, nenzymes, reps, ws, blocksize, categs, maxnewsites;
boolean bootstrap, permute, ild, lockhart, jackknife, regular, xml, nexus,
  weights, categories, factors, enzymes, all, justwts, progress, mixture,
  firstrep, ancvar;
double fracsample;
datatype data;
seqtype seq;
steptr oldweight, where, how_many, newwhere, newhowmany,
  newerwhere, newerhowmany, factorr, newerfactor, mixdata, ancdata;
steptr *charorder;
Char *factor;
long *alleles;
Char **nodep;
double **nodef;
long **sppord;
longer seed;


void getoptions()
{
  /* interactively set options */
  long reps0;
  long inseed, inseed0, loopcount, loopcount2;
  Char ch;
  boolean done1;

  data = seqs;
  seq = dna;
  bootstrap = true;
  jackknife = false;
  permute = false;
  ild = false;
  lockhart = false;
  blocksize = 1;
  regular = true;
  fracsample = 1.0;
  all = false;
  reps = 100;
  weights = false;
  mixture = false;
  ancvar = false;
  categories = false;
  justwts = false;
  printdata = false;
  dotdiff = true;
  progress = true;
  interleaved = true;
  xml = false;
  nexus = false;
  factors = false;
  loopcount = 0;
  for (;;) {
    cleerhome();
    printf("\nBootstrapping algorithm, version %s\n\n",VERSION);
    printf("Settings for this run:\n");
    printf("  D      Sequence, Morph, Rest., Gene Freqs?  %s\n",
           (data == seqs       ) ? "Molecular sequences"      :
           (data == morphology ) ? "Discrete Morphology"      :
           (data == restsites)   ? "Restriction Sites"        :
           (data == genefreqs)   ? "Gene Frequencies" : "");
    if (data == restsites)
      printf("  E                       Number of enzymes?  %s\n",
             enzymes ? "Present in input file" :
                       "Not present in input file");
    if (data == genefreqs)
      printf("  A       All alleles present at each locus?  %s\n",
             all ? "Yes" : "No, one absent at each locus");
    if (data == morphology)
      printf("  F                 Use factors information?  %s\n",
           factors ? "Yes" : "No");

    printf("  J  Bootstrap, Jackknife, Permute, Rewrite?  %s\n",
           regular && jackknife ? "Delete-half jackknife" :
           (!regular) && jackknife ? "Delete-fraction jackknife" :
           permute   ? "Permute species for each character" :
           ild ? "Permute character order" :
           lockhart ? "Permute within species" :
           regular && bootstrap ? "Bootstrap" :
           (!regular) && bootstrap ? "Partial bootstrap" :
           "Rewrite data");
    if (bootstrap || jackknife) {
      printf("  ");
      putchar('%');
      printf("    Regular or altered sampling fraction?  ");
      if (regular)
        printf("regular\n");
      else {
        if (fabs(fracsample*100 - (int)(fracsample*100)) > 0.01) {
          if (fracsample < 1)
            printf("%2.1lf", 100.0*fracsample);
          else
            printf("%3.1lf", 100.0*fracsample);
        } else { if (fracsample < 1)
            printf("%2.0lf", 100.0*fracsample);
          else
            printf("%3.0lf", 100.0*fracsample);
        }
        putchar('%');
        printf(" sampled\n");
      }
    }
    if ((data == seqs)
        && !(jackknife || permute || bootstrap || ild || lockhart)) {
      printf("  P     PHYLIP, NEXUS, or XML output format?  %s\n",
           nexus ? "NEXUS" : xml ? "XML" : "PHYLIP");
      if (xml || ((data == seqs) && nexus)) {
        printf("  S             Type of molecular sequences?  " );
        switch (seq) {
          case (dna) : printf("DNA\n"); break;
          case (rna) : printf("RNA\n"); break;
          case (protein) : printf("Protein\n"); break;
        }
      }
    }
    if ((data == morphology) && !(jackknife || permute || ild
                                  || lockhart || bootstrap))
      printf("  P           PHYLIP or NEXUS output format?  %s\n",
                                         nexus ? "NEXUS" : "PHYLIP");
    if (bootstrap) {
      if (blocksize > 1)
      printf("  B      Block size for block-bootstrapping?  %ld\n", blocksize);
      else
        printf("  B      Block size for block-bootstrapping?  %ld (regular bootstrap)\n", blocksize);
    }
    if (bootstrap || jackknife || permute || ild || lockhart)
      printf("  R                     How many replicates?  %ld\n", reps);
    if (jackknife || bootstrap || permute) {
      printf("  W              Read weights of characters?  %s\n",
          (weights ? "Yes" : "No"));
      if(data == morphology){
        printf("  X                       Read mixture file?  %s\n",
               (mixture ? "Yes" : "No"));
        printf("  N                     Read ancestors file?  %s\n",
               (ancvar ? "Yes" : "No"));
      }
      if (data == seqs)
        printf("  C                Read categories of sites?  %s\n",
          (categories ? "Yes" : "No"));
      if ((!permute)) {
        printf("  S     Write out data sets or just weights?  %s\n",
          (justwts ? "Just weights" : "Data sets"));
      }
    }
    if (data == seqs || data == restsites)
      printf("  I             Input sequences interleaved?  %s\n",
             interleaved ? "Yes" : "No, sequential");
    printf("  0      Terminal type (IBM PC, ANSI, none)?  %s\n",
           ibmpc ? "IBM PC" : ansi  ? "ANSI" : "(none)");
    printf("  1       Print out the data at start of run  %s\n",
           printdata ? "Yes" : "No");
    if (printdata)
      printf("  .     Use dot-differencing to display them  %s\n",
           dotdiff ? "Yes" : "No");
    printf("  2     Print indications of progress of run  %s\n",
           progress ? "Yes" : "No");
    printf("\n  Y to accept these or type the letter for one to change\n");
#ifdef WIN32
    phyFillScreenColor();
#endif
    scanf("%c%*[^\n]", &ch);
    getchar();
    uppercase(&ch);
    if (ch == 'Y')
           break;
    if (
        (bootstrap && (strchr("ABCDEFSJPRWXNI%1.20",ch) != NULL)) ||
        (jackknife && (strchr("ACDEFSJPRWXNI%1.20",ch) != NULL)) ||
        (permute && (strchr("ACDEFSJPRWXNI%1.20",ch) != NULL)) ||
        ((ild || lockhart) && (strchr("ACDEFSJPRXNI%1.20",ch) != NULL)) ||
        ((!(bootstrap || jackknife || permute || ild || lockhart)) &&
         ((!xml) && (strchr("ADEFJPI1.20",ch) != NULL))) ||
        (((data == morphology) || (data == seqs))
         && (nexus || xml) && (strchr("ADEFJPSI1.20",ch) != NULL))
        ) {
      switch (ch) {
        
      case 'D':
        if (data == genefreqs)
          data = seqs;
        else
          data = (datatype)((long)data + 1);
        break;
        
      case 'A':
        all = !all;
        break;
        
      case 'E':
        enzymes = !enzymes;
        break;
        
      case 'J':
        if (permute) {
          permute = false;
          ild = true;
        } else if (ild) {
            ild = false;
            lockhart = true;
          } else if (lockhart)
              lockhart = false;
            else if (jackknife) {
                jackknife = false;
                permute = true;
            } else if (bootstrap) {
                  bootstrap = false;
                  jackknife = true;
              } else
                  bootstrap = true;
        break;
        
      case '%':
        regular = !regular;
        if (!regular) {
          loopcount2 = 0;
          do {
            printf("Samples as percentage of");
            if ((data == seqs) || (data == restsites))
              printf(" sites?\n");
            if (data == morphology)
            printf(" characters?\n");
            if (data == genefreqs)
              printf(" loci?\n");
            scanf("%lf%*[^\n]", &fracsample);
            getchar();
            done1 = (fracsample > 0.0);
            if (!done1) {
              printf("BAD NUMBER: must be positive\n");
            }
            fracsample = fracsample/100.0;
            countup(&loopcount2, 10);
          } while (done1 != true);
        }
        break;

      case 'P':
        if (data == seqs) {
          if (!xml && !nexus)
            nexus = true;
          else {
            if (nexus) {
              nexus = false;
              xml = true;
              }
            else xml = false;
            }
          }
        if (data == morphology) {
          nexus = !nexus;
          xml = false;
          }
        break;
        
      case 'S':
        if(jackknife || permute || bootstrap || ild || lockhart){
          justwts = !justwts;          
        } else {
          switch (seq) {
          case (dna): seq = rna; break;
          case (rna): seq = protein; break;
          case (protein): seq = dna; break;
          }
        }
        break;
        
      case 'B':
        loopcount2 = 0;
        do {
          printf("Block size?\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
          scanf("%ld%*[^\n]", &blocksize);
          getchar();
          done1 = (blocksize > 0);
          if (!done1) {
            printf("BAD NUMBER: must be positive\n");
          }
          countup(&loopcount2, 10);
        } while (done1 != true);
        break;

      case 'R':
        reps0 = reps;
        loopcount2 = 0;
        do {
          printf("Number of replicates?\n");
#ifdef WIN32
          phyFillScreenColor();
#endif
          scanf("%ld%*[^\n]", &reps);
          getchar();
          done1 = (reps > 0);
          if (!done1) {
            printf("BAD NUMBER: must be positive\n");
            reps = reps0;
          }
          countup(&loopcount2, 10);
        } while (done1 != true);
        break;

      case 'W':
        weights = !weights;
        break;

      case 'X':
        mixture = !mixture;
        break;

      case 'N':
        ancvar = !ancvar;
        break;

      case 'C':
        categories = !categories;
        break;

      case 'F':
        factors = !factors;
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
        
      case '.':
        dotdiff = !dotdiff;
        break;
        
      case '2':
        progress = !progress;
        break;
      }
    } else
      printf("Not a possible option!\n");
    countup(&loopcount, 100);
  }
  if (bootstrap || jackknife) {
    if (jackknife && regular)
      fracsample = 0.5;
    if (bootstrap && regular)
      fracsample = 1.0;
  }
  if (bootstrap || jackknife || permute || ild || lockhart)
    initseed(&inseed, &inseed0, seed);
  xml = xml && (data == seqs);
  categories = categories && (data == seqs);
  mixture = mixture && (data == morphology);
  ancvar = ancvar && (data == morphology);
}  /* getoptions */


void seqboot_inputnumbers()
{
  /* read numbers of species and of sites */
  long i;

  fscanf(infile, "%ld%ld", &spp, &sites);
  loci = sites;
  maxalleles = 1;
  if (data == restsites && enzymes)
    fscanf(infile, "%ld", &nenzymes);
  if (data == genefreqs) {
    alleles = (long *)Malloc(sites*sizeof(long));
    scan_eoln(infile);
    sites = 0;
    for (i = 0; i < (loci); i++) {
      if (eoln(infile)) 
        scan_eoln(infile);
      fscanf(infile, "%ld", &alleles[i]);
      if (alleles[i] > maxalleles)
         maxalleles = alleles[i];
      if (all)
         sites += alleles[i];
      else
         sites += alleles[i] - 1;
    }
    if (!all)
      maxalleles--;
    scan_eoln(infile);
  }
}  /* seqboot_inputnumbers */


void seqboot_inputfactors()
{
  long i, j;
  Char ch, prevch;

  prevch = ' ';
  j = 0;
  for (i = 0; i < (sites); i++) {
    do {
      if (eoln(factfile)) 
        scan_eoln(factfile);
      ch = gettc(factfile);
    } while (ch == ' ');
    if (ch != prevch)
      j++;
    prevch = ch;
    factorr[i] = j;
  }
  scan_eoln(factfile);
}  /* seqboot_inputfactors */


void inputoptions()
{
  /* input the information on the options */
  long weightsum, maxfactsize, i, j, k, l, m;

  if (data == genefreqs) {
    k = 0;
    l = 0;
    for (i = 0; i < (loci); i++) {
      if (all)
        m = alleles[i];
      else
        m = alleles[i] - 1;
      k++;
      for (j = 1; j <= m; j++) {
        l++;
        factorr[l - 1] = k;
      }
    }
  } else {
    for (i = 1; i <= (sites); i++)
      factorr[i - 1] = i;
  }
  if(factors){
    seqboot_inputfactors();
  }
  for (i = 0; i < (sites); i++)
    oldweight[i] = 1;
  if (weights)
    inputweights2(0, sites, &weightsum, oldweight, &weights, "seqboot");
  if (factors && printdata) {
    for(i = 0; i < sites; i++)
      factor[i] = (char)('0' + (factorr[i]%10));
    printfactors(outfile, sites, factor, " (least significant digit)");
  }
  if (weights && printdata)
    printweights(outfile, 0, sites, oldweight, "Sites");
  for (i = 0; i < (loci); i++)
    how_many[i] = 0;
  for (i = 0; i < (loci); i++)
    where[i] = 0;
  for (i = 1; i <= (sites); i++) {
    how_many[factorr[i - 1] - 1]++;
    if (where[factorr[i - 1] - 1] == 0)
      where[factorr[i - 1] - 1] = i;
  }
  groups = factorr[sites - 1];
  newgroups = 0;
  newsites = 0;
  maxfactsize = 0;
  for(i = 0 ; i < loci ; i++){
    if(how_many[i] > maxfactsize){
      maxfactsize = how_many[i];
    }
  }
  maxnewsites = groups * maxfactsize;
  allocnew();
  for (i = 0; i < (groups); i++) {
    if (oldweight[where[i] - 1] > 0) {
      newgroups++;
      newsites += how_many[i];
      newwhere[newgroups - 1] = where[i];
      newhowmany[newgroups - 1] = how_many[i];
    }
  }
}  /* inputoptions */


void seqboot_inputdata()
{
  /* input the names and sequences for each species */
  long i, j, k, l, m, n, basesread, basesnew=0;
  double x;
  Char charstate;
  boolean allread, done;

  if (data == genefreqs) {
    nodef = (double **)Malloc(spp*sizeof(double *));
    for (i = 0; i < (spp); i++)
      nodef[i] = (double *)Malloc(sites*sizeof(double));
  } else {
    nodep = (Char **)Malloc(spp*sizeof(Char *));
    for (i = 0; i < (spp); i++)
      nodep[i] = (Char *)Malloc(sites*sizeof(Char));
  }
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;
  if (printdata) {
    fprintf(outfile, "\nBootstrapping algorithm, version %s\n\n\n",VERSION);
    if (bootstrap)  {
      if (blocksize > 1) {
        if (regular)      
      fprintf(outfile, "Block-bootstrap with block size %ld\n\n", blocksize);
        else
          fprintf(outfile, "Partial (%2.0f%%) block-bootstrap with block size %ld\n\n",
                  100*fracsample, blocksize);
      } else {
        if (regular)
          fprintf(outfile, "Bootstrap\n\n");
        else 
          fprintf(outfile, "Partial (%2.0f%%) bootstrap\n\n", 100*fracsample);
      }
    } else {
      if (jackknife) {
        if (regular)
          fprintf(outfile, "Delete-half Jackknife\n\n");
        else
    fprintf(outfile, "Delete-%2.0f%% Jackknife\n\n", 100*(1.0-fracsample));
      } else {
        if (permute) {
          fprintf(outfile, "Species order permuted separately for each");
          if (data == genefreqs)
            fprintf(outfile, " locus\n\n");
          if (data == seqs)
            fprintf(outfile, " site\n\n");
          if (data == morphology)
            fprintf(outfile, " character\n\n");
          if (data == restsites)
            fprintf(outfile, " site\n\n");
        }
        else {
          if (ild) {
            if (data == genefreqs)
              fprintf(outfile, "Locus");
            if (data == seqs)
              fprintf(outfile, "Site");
            if (data == morphology)
              fprintf(outfile, "Character");
            if (data == restsites)
              fprintf(outfile, "Site");
            fprintf(outfile, " order permuted\n\n");
          } else {
            if (lockhart)
              if (data == genefreqs)
                fprintf(outfile, "Locus");
              if (data == seqs)
                fprintf(outfile, "Site");
              if (data == morphology)
                fprintf(outfile, "Character");
              if (data == restsites)
                fprintf(outfile, "Site");
         fprintf(outfile, " order permuted separately for each species\n\n");
          }
        }
      }
    }
    if (data == genefreqs)
      fprintf(outfile, "%3ld species, %3ld  loci\n\n", spp, loci);
    else {
      fprintf(outfile, "%3ld species, ", spp);
      if (data == seqs)
        fprintf(outfile, "%3ld  sites\n\n", sites);
        else if (data == morphology)
          fprintf(outfile, "%3ld  characters\n\n", sites);
          else if (data == restsites)
            fprintf(outfile, "%3ld  sites\n\n", sites);
    }
    fprintf(outfile, "Name");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "Data\n");
    fprintf(outfile, "----");
    for (i = 1; i <= j; i++)
      putc(' ', outfile);
    fprintf(outfile, "----\n\n");
  }
  interleaved = (interleaved && ((data == seqs) || (data == restsites)));
  if (data == genefreqs) {
    for (i = 1; i <= (spp); i++) {
      initname(i - 1);
      j = 1;
      while (j <= sites && !eoff(infile)) {
        if (eoln(infile)) 
          scan_eoln(infile);
        fscanf(infile, "%lf", &x);
        if ((unsigned)x > 1.0) {
          printf("GENE FREQ OUTSIDE [0,1] in species %ld\n", i);
          exxit(-1);
        } else {
          nodef[i - 1][j - 1] = x;
          j++;
        }
      }
      scan_eoln(infile);
    }
    return;
  }
  basesread = 0;
  allread = false;
  while (!allread) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile)) 
      scan_eoln(infile);
    i = 1;
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i-1);
      j = interleaved ? basesread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) ||eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' ||
              (data == seqs && charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          j++;
          if (charstate == '.')
            charstate = nodep[0][j-1];
          nodep[i-1][j-1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < sites) 
          scan_eoln(infile);
        else if (j == sites)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;
      scan_eoln(infile);
      if ((interleaved && j != basesnew) || ((!interleaved) && j != sites)){
        printf("\n\nERROR: sequences out of alignment at site %ld", j+1);
        printf(" of species %ld\n\n", i);
        exxit(-1);}
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > spp);
  }
  if (!printdata)
    return;
  if (data == genefreqs)
    m = (sites - 1) / 8 + 1;
  else
    m = (sites - 1) / 60 + 1;
  for (i = 1; i <= m; i++) {
    for (j = 0; j < spp; j++) {
      for (k = 0; k < nmlngth; k++)
        putc(nayme[j][k], outfile);
      fprintf(outfile, "   ");
      if (data == genefreqs)
        l = i * 8;
      else
        l = i * 60;
      if (l > sites)
        l = sites;
      if (data == genefreqs)
        n = (i - 1) * 8;
      else
        n = (i - 1) * 60;
      for (k = n; k < l; k++) {
        if (data == genefreqs)
          fprintf(outfile, "%8.5f", nodef[j][k]);
        else {
          if (j + 1 > 1 && nodep[j][k] == nodep[0][k])
            charstate = '.';
          else
            charstate = nodep[j][k];
          putc(charstate, outfile);
          if ((k + 1) % 10 == 0 && (k + 1) % 60 != 0)
            putc(' ', outfile);
        
        }
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  putc('\n', outfile);
}  /* seqboot_inputdata */


void allocrest()
{ /* allocate memory for bookkeeping arrays */

  oldweight = (steptr)Malloc(sites*sizeof(long));
  weight = (steptr)Malloc(sites*sizeof(long));
  if (categories)
    category = (steptr)Malloc(sites*sizeof(long));
  if (mixture)
    mixdata = (steptr)Malloc(sites*sizeof(long));
  if (ancvar)
    ancdata = (steptr)Malloc(sites*sizeof(long));
  where = (steptr)Malloc(loci*sizeof(long));
  how_many = (steptr)Malloc(loci*sizeof(long));
  factor = (Char *)Malloc(sites*sizeof(Char));
  factorr = (steptr)Malloc(sites*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
}  /* allocrest */

void allocnew(void)
{ /* allocate memory for arrays that depend on the lenght of the 
     output sequence*/
  long i;

  newwhere = (steptr)Malloc(loci*sizeof(long));
  newhowmany = (steptr)Malloc(loci*sizeof(long));
  newerwhere = (steptr)Malloc(loci*sizeof(long));
  newerhowmany = (steptr)Malloc(loci*sizeof(long));
  newerfactor = (steptr)Malloc(maxnewsites*maxalleles*sizeof(long));
  charorder = (steptr *)Malloc(spp*sizeof(steptr));
  for (i = 0; i < spp; i++)
    charorder[i] = (steptr)Malloc(maxnewsites*sizeof(long));
}

void doinput(int argc, Char *argv[])
{ /* reads the input data */
  getoptions();
  seqboot_inputnumbers();
  allocrest();
  if (weights)
    openfile(&weightfile,WEIGHTFILE,"input weight file",
               "r",argv[0],weightfilename);
  if (mixture){
    openfile(&mixfile,MIXFILE,"mixture file", "r",argv[0],mixfilename);
    openfile(&outmixfile,"outmixture","output mixtures file","w",argv[0],
               outmixfilename);
    seqboot_inputaux(mixdata, mixfile);
  }
  if (ancvar){
    openfile(&ancfile,ANCFILE,"ancestor file", "r",argv[0],ancfilename);
    openfile(&outancfile,"outancestors","output ancestors file","w",argv[0],
               outancfilename);
    seqboot_inputaux(ancdata, ancfile);
  }
  if (categories) {
    openfile(&catfile,CATFILE,"input category file","r",argv[0],catfilename);
    openfile(&outcatfile,"outcategories","output category file","w",argv[0],
               outcatfilename);
    inputcategs(0, sites, category, 9, "SeqBoot");
  }
  if (factors){
    openfile(&factfile,FACTFILE,"factors file","r",argv[0],factfilename);
    openfile(&outfactfile,"outfactors","output factors file","w",argv[0],
               outfactfilename);
  }
  if (justwts && !permute)
    openfile(&outweightfile,"outweights","output weight file",
               "w",argv[0],outweightfilename);
  else {
    openfile(&outfile,OUTFILE,"output data file","w",argv[0],outfilename);
  }
  inputoptions();
  seqboot_inputdata();
}  /* doinput */


void bootweights()
{ /* sets up weights by resampling data */
  long i, j, k, blocks;
  double p, q, r;

  ws = newgroups;
  for (i = 0; i < (ws); i++)
    weight[i] = 0;
  if (jackknife) {
    if (fabs(newgroups*fracsample - (long)(newgroups*fracsample+0.5))
       > 0.00001) {
      if (randum(seed)
          < (newgroups*fracsample - (long)(newgroups*fracsample))
            /((long)(newgroups*fracsample+1.0)-(long)(newgroups*fracsample)))
        q = (long)(newgroups*fracsample)+1;
      else
        q = (long)(newgroups*fracsample);
    } else
      q = (long)(newgroups*fracsample+0.5);
    r = newgroups;
    p = q / r;
    ws = 0;
    for (i = 0; i < (newgroups); i++) {
      if (randum(seed) < p) {
        weight[i]++;
        ws++;
        q--;
      }
      r--;
      if (i + 1 < newgroups)
        p = q / r;
    }
  } else if (permute) {
    for (i = 0; i < (newgroups); i++)
      weight[i] = 1;
  } else if (bootstrap) {
    blocks = fracsample * newgroups / blocksize;
    for (i = 1; i <= (blocks); i++) {
      j = (long)(newgroups * randum(seed)) + 1;
      for (k = 0; k < blocksize; k++) {
        weight[j - 1]++;
        j++;
        if (j > newgroups)
          j = 1;
      }
    }
  } else             /* case of rewriting data */
    for (i = 0; i < (newgroups); i++)
      weight[i] = 1;
  for (i = 0; i < (newgroups); i++)
    newerwhere[i] = 0;
  for (i = 0; i < (newgroups); i++)
    newerhowmany[i] = 0;
  newergroups = 0;
  newersites = 0;
  for (i = 0; i < (newgroups); i++) {
    for (j = 1; j <= (weight[i]); j++) {
      newergroups++;
      for (k = 1; k <= (newhowmany[i]); k++) {
        newersites++;
        newerfactor[newersites - 1] = newergroups;
      }
      newerwhere[newergroups - 1] = newwhere[i];
      newerhowmany[newergroups - 1] = newhowmany[i];
    }
  }
}  /* bootweights */


void sppermute(long n)
{ /* permute the species order as given in array sppord */
  long i, j, k;

  for (i = 1; i <= (spp - 1); i++) {
    k = (long)((i+1) * randum(seed));
    j = sppord[n - 1][i];
    sppord[n - 1][i] = sppord[n - 1][k];
    sppord[n - 1][k] = j;
  }
}  /* sppermute */


void charpermute(long m, long n)
{ /* permute the n+1 characters of species m+1 */
  long i, j, k;

  for (i = 1; i <= (n-1); i++) {
    k = (long)((i+1) * randum(seed));
    j = charorder[m][i];
    charorder[m][i] = charorder[m][k];
    charorder[m][k] = j;
  }
} /* charpermute */


void writedata()
{
  /* write out one set of bootstrapped sequences */
  long i, j, k, l, m, n, n2;
  double x;
  Char charstate;

  sppord = (long **)Malloc(newergroups*sizeof(long *));
  for (i = 0; i < (newergroups); i++)
    sppord[i] = (long *)Malloc(spp*sizeof(long));
  for (j = 1; j <= spp; j++)
    sppord[0][j - 1] = j;
  for (i = 1; i < newergroups; i++) {
    for (j = 1; j <= (spp); j++)
      sppord[i][j - 1] = sppord[i - 1][j - 1];
  }
  if (!justwts || permute) {
    if (data == restsites && enzymes)
      fprintf(outfile, "%5ld %5ld% 4ld\n", spp, newergroups, nenzymes);
    else if (data == genefreqs)
      fprintf(outfile, "%5ld %5ld\n", spp, newergroups);
    else {
      if ((data == seqs)
          && !(bootstrap || jackknife || permute || ild || lockhart) && xml)
        fprintf(outfile, "<alignment>\n");
      else 
        if (!(bootstrap || jackknife || permute || ild || lockhart) && nexus) {
          fprintf(outfile, "#NEXUS\n");
          fprintf(outfile, "BEGIN DATA;\n");
          fprintf(outfile, "  DIMENSIONS NTAX=%ld NCHAR=%ld;\n",
                    spp, newersites);
          fprintf(outfile, "  FORMAT");
          if (interleaved)
            fprintf(outfile, " interleave=yes");
          else
            fprintf(outfile, " interleave=no");
          fprintf(outfile, " DATATYPE=");
          if (data == seqs) {
            switch (seq) { 
              case (dna): fprintf(outfile, "DNA missing=N gap=-"); break;
              case (rna): fprintf(outfile, "RNA missing=N gap=-"); break;
              case (protein):
                fprintf(outfile, "protein missing=? gap=-");
                break;
              }
          }
          if (data == morphology)
            fprintf(outfile, "STANDARD");
          fprintf(outfile, ";\n  MATRIX\n");
          }
        else fprintf(outfile, "%5ld %5ld\n", spp, newersites);
    }
    if (data == genefreqs) {
      for (i = 0; i < (newergroups); i++)
        fprintf(outfile, " %3ld", alleles[factorr[newerwhere[i] - 1] - 1]);
      putc('\n', outfile);
    }
  }
  l = 1;
  if ((!(bootstrap || jackknife || permute || ild || lockhart || nexus))
       && ((data == seqs) || (data == restsites))) {
    interleaved = !interleaved;
    if (!(bootstrap || jackknife || permute || ild || lockhart) && xml)
      interleaved = false;
  }
  if (interleaved)
    m = 60;
  else
    m = newergroups;
  do {
    if (m > newergroups)
      m = newergroups;
    for (j = 0; j < spp; j++) {
      n = 0;
      if ((l == 1) || (interleaved && nexus)) {
        if (!(bootstrap || jackknife || permute || ild || lockhart) && xml) {
          fprintf(outfile, "   <sequence");
          switch (seq) {
            case (dna): fprintf(outfile, " type=\"dna\""); break;
            case (rna): fprintf(outfile, " type=\"rna\""); break;
            case (protein): fprintf(outfile, " type=\"protein\""); break;
          }
          fprintf(outfile, ">\n");
          fprintf(outfile, "      <name>");
        }
        n2 = nmlngth;
        if (!(bootstrap || jackknife || permute || ild || lockhart)
            && (xml || nexus)) {
          while (nayme[j][n2-1] == ' ')
            n2--;
        }
        if (nexus)
          fprintf(outfile, "  ");
        for (k = 0; k < n2; k++)
          if (nexus && (nayme[j][k] == ' ') && (k < n2))
            putc('_', outfile);
          else
            putc(nayme[j][k], outfile);
        if (!(bootstrap || jackknife || permute || ild || lockhart) && xml)
          fprintf(outfile, "</name>\n      <data>");
      } else {
        if (!(bootstrap || jackknife || permute || ild || lockhart) && xml) {
          fprintf(outfile, "      ");
        }
      }
      if (!xml) {
        for (k = 0; k < nmlngth-n2; k++)
          fprintf(outfile, " "); 
        fprintf(outfile, " "); 
      }
      for (k = l - 1; k < m; k++) {
        if (permute && j + 1 == 1)
          sppermute(newerfactor[n]);    /* we can assume chars not permuted */
        for (n2 = -1; n2 <= (newerhowmany[k] - 2); n2++) {
          n++;
          if (data == genefreqs) {
            if (n > 1 && (n & 7) == 1)
              fprintf(outfile, "\n              ");
            x = nodef[sppord[newerfactor[charorder[j][n - 1]] - 1][j] - 1]
                    [newerwhere[charorder[j][k]] + n2];
            fprintf(outfile, "%8.5f", x);
          } else {
            if (!(bootstrap || jackknife || permute || ild || lockhart) && xml
                && (n > 1) && (n % 60 == 1))
              fprintf(outfile, "\n            ");
            else if (!nexus && !interleaved && (n > 1) && (n % 60 == 1))
                fprintf(outfile, "\n           ");
            charstate = nodep[sppord[newerfactor[charorder[j][n - 1]] - 1]
                             [j] - 1][newerwhere[charorder[j][k]] + n2];
            putc(charstate, outfile);
            if (n % 10 == 0 && n % 60 != 0)
              putc(' ', outfile);
          }
        }
      }
      if (!(bootstrap || jackknife || permute || ild || lockhart ) && xml) {
        fprintf(outfile, "</data>\n   </sequence>\n");
      }
      putc('\n', outfile);
    }
    if (interleaved) {
      if ((m <= newersites) && (newersites > 60))
        putc('\n', outfile);
      l += 60;
      m += 60;
    }
  } while (interleaved && l <= newersites);
  if ((data == seqs) &&
      (!(bootstrap || jackknife || permute || ild || lockhart) && xml))
    fprintf(outfile, "</alignment>\n");
  if (!(bootstrap || jackknife || permute || ild || lockhart) && nexus)
    fprintf(outfile, "  ;\nEND;\n");
  for (i = 0; i < (newergroups); i++)
    free(sppord[i]);
  free(sppord);
}  /* writedata */


void writeweights()
{ /* write out one set of post-bootstrapping weights */
  long j, k, l, m, n, o;

  j = 0;
  l = 1;
  if (interleaved)
    m = 60;
  else
    m = sites;
  do {
    if(m > sites)
      m = sites;
    n = 0;
    for (k = l - 1; k < m; k++) {
      for(o = 0 ; o < how_many[k] ; o++){
        if(oldweight[k]==0){
          fprintf(outweightfile, "0");
          j++;
        }
        else{
          if (weight[k-j] < 10)
            fprintf(outweightfile, "%c", (char)('0'+weight[k-j]));
          else
            fprintf(outweightfile, "%c", (char)('A'+weight[k-j]-10));
          n++;
          if (!interleaved && n > 1 && n % 60 == 1) {
            fprintf(outweightfile, "\n");
            if (n % 10 == 0 && n % 60 != 0)
              putc(' ', outweightfile);
          }
        }
      }
    }
    putc('\n', outweightfile);
    if (interleaved) {
      l += 60;
      m += 60;
    }
  } while (interleaved && l <= sites);
}  /* writeweights */


void writecategories()
{
  /* write out categories for the bootstrapped sequences */
  long k, l, m, n, n2;
  Char charstate;
  if(justwts){
    if (interleaved)
      m = 60;
    else
      m = sites;
    l=1;
    do {
      if(m > sites)
        m = sites;
      n=0;
      for(k=l-1 ; k < m ; k++){
        n++;
        if (!interleaved && n > 1 && n % 60 == 1)
          fprintf(outcatfile, "\n ");
        charstate =  '0' + category[k];
        putc(charstate, outcatfile);
      }
      if (interleaved) {
        l += 60;
        m += 60;
      }
    }while(interleaved && l <= sites);
    fprintf(outcatfile, "\n");
    return;
  }

  l = 1;
  if (interleaved)
    m = 60;
  else
    m = newergroups;
  do {
    if (m > newergroups)
      m = newergroups;
    n = 0;
    for (k = l - 1; k < m; k++) {
      for (n2 = -1; n2 <= (newerhowmany[k] - 2); n2++) {
        n++;
        if (!interleaved && n > 1 && n % 60 == 1)
          fprintf(outcatfile, "\n ");
        charstate = '0' + category[newerwhere[k] + n2];
        putc(charstate, outcatfile);
        if (n % 10 == 0 && n % 60 != 0)
          putc(' ', outcatfile);
      }
    }
    if (interleaved) {
      l += 60;
      m += 60;
    }
  } while (interleaved && l <= newersites);
  fprintf(outcatfile, "\n");
}  /* writecategories */


void writeauxdata(steptr auxdata, FILE *outauxfile)
{
  /* write out auxiliary option data (mixtures, ancestors, ect) to
     appropriate file.  Samples parralel to data, or just gives one
     output entry if justwts is true */
  long k, l, m, n, n2;
  Char charstate;

  /* if we just output weights (justwts), and this is first set
     just output the data unsampled */
  if(justwts){
    if(firstrep){
      if (interleaved)
        m = 60;
      else
        m = sites;
      l=1;
      do {
        if(m > sites)
          m = sites;
        n = 0;
        for(k=l-1 ; k < m ; k++){
        n++;
        if (!interleaved && n > 1 && n % 60 == 1)
          fprintf(outauxfile, "\n ");
          charstate = auxdata[k];
          putc(charstate, outauxfile);
        }
        if (interleaved) {
          l += 60;
          m += 60;
        }
      }while(interleaved && l <= sites);
      fprintf(outauxfile, "\n");
    }
    return;
  }

  l = 1;
  if (interleaved)
    m = 60;
  else
    m = newergroups;
  do {
    if (m > newergroups)
      m = newergroups;
    n = 0;
    for (k = l - 1; k < m; k++) {
      for (n2 = -1; n2 <= (newerhowmany[k] - 2); n2++) {
        n++;
        if (!interleaved && n > 1 && n % 60 == 1)
          fprintf(outauxfile, "\n ");
        charstate = auxdata[newerwhere[k] + n2];
        putc(charstate, outauxfile);
        if (n % 10 == 0 && n % 60 != 0)
          putc(' ', outauxfile);
      }
    }
    if (interleaved) {
      l += 60;
      m += 60;
    }
  } while (interleaved && l <= newersites);
  fprintf(outauxfile, "\n");
}  /* writeauxdata */

void writefactors(void)
{
  long k, l, m, n, prevfact, writesites;
  char symbol;
  steptr wfactor;

  if(!justwts || firstrep){
    if(justwts){
      writesites = sites;
      wfactor = factorr;
    } else {
      writesites = newersites;
      wfactor = newerfactor;
    }
    prevfact = wfactor[0];
    symbol = '+';
    if (interleaved)
      m = 60;
    else
      m = writesites;
    l=1;
    do {
      if(m > writesites)
        m = writesites;
      n = 0;
      for(k=l-1 ; k < m ; k++){
        n++;
        if (!interleaved && n > 1 && n % 60 == 1)
          fprintf(outfactfile, "\n ");
        if(prevfact != wfactor[k]){
          symbol = (symbol == '+') ? '-' : '+';
          prevfact = wfactor[k];
        }
        putc(symbol, outfactfile);
        if (n % 10 == 0 && n % 60 != 0)
          putc(' ', outfactfile);
      }
      if (interleaved) {
        l += 60;
        m += 60;
      }
    }while(interleaved && l <= writesites);
    fprintf(outfactfile, "\n");
  }
} /* writefactors */


void bootwrite()
{ /* does bootstrapping and writes out data sets */
  long i, j, rr, repdiv10;

  if (!(bootstrap || jackknife || permute || ild || lockhart))
    reps = 1;
  repdiv10 = reps / 10;
  if (repdiv10 < 1)
    repdiv10 = 1;
  if (progress)
    putchar('\n');
  for (rr = 1; rr <= (reps); rr++) {
    for (i = 0; i < spp; i++)
      for (j = 0; j < maxnewsites; j++)
        charorder[i][j] = j;
    if(rr==1)
      firstrep = true;
    else
      firstrep = false;
    if (ild) {
      charpermute(0, maxnewsites);
      for (i = 1; i < spp; i++)
        for (j = 0; j < maxnewsites; j++)
          charorder[i][j] = charorder[0][j];
    }
    if (lockhart)
      for (i = 0; i < spp; i++)
        charpermute(i, maxnewsites);
    bootweights();
    if (!justwts || permute || ild || lockhart)
      writedata();
    if (justwts && !(permute || ild || lockhart))
      writeweights();
    if (categories)
      writecategories();
    if (factors)
      writefactors();
    if (mixture)
      writeauxdata(mixdata, outmixfile);
    if (ancvar)
      writeauxdata(ancdata, outancfile);
    if (progress && (bootstrap || jackknife || permute || ild || lockhart)
          && ((reps < 10) || rr % repdiv10 == 0)) {
      printf("completed replicate number %4ld\n", rr);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
  }
  if (progress) {
    if (justwts)
      printf("\nOutput weights written to file \"%s\"\n\n", outweightfilename);
    else
      printf("\nOutput written to file \"%s\"\n\n", outfilename);
  }
}  /* bootwrite */


void seqboot_inputaux(steptr dataptr, FILE* auxfile)
{ /* input auxiliary option data (mixtures, ancestors, ect) for 
     new style input, assumes that data is correctly formated
     in input files*/
  long i, j, k;
  Char ch;

  j = 0;
  k = 1;
  for (i = 0; i < (sites); i++) {
    do {
      if (eoln(auxfile))
        scan_eoln(auxfile);
      ch = gettc(auxfile);
      if (ch == '\n')
        ch = ' ';
    } while (ch == ' ');
    dataptr[i] = ch;
  }
  scan_eoln(auxfile);
}  /* seqboot_inputaux */


int main(int argc, Char *argv[])
{  /* Read in sequences or frequencies and bootstrap or jackknife them */
#ifdef MAC
  argc = 1;                /* macsetup("SeqBoot","");                */
  argv[0] = "SeqBoot";
#endif
  init(argc,argv);
  openfile(&infile, INFILE, "input file", "r", argv[0], infilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  doinput(argc, argv);
  bootwrite();
  FClose(infile);
  if (weights)
    FClose(weightfile);
  if (categories) {
    FClose(catfile);
    FClose(outcatfile);
  }
  if(mixture)
    FClose(outmixfile);
  if(ancvar)
    FClose(outancfile);
  if (justwts && !permute) {
    FClose(outweightfile);
  }
  else
    FClose(outfile);
#ifdef MAC
  fixmacfile(outfilename);
  if (justwts && !permute)
    fixmacfile(outweightfilename);
  if (categories)
    fixmacfile(outcatfilename);
  if (mixture)
    fixmacfile(outmixfilename);
#endif
  printf("Done.\n\n");
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}
