/***********************************************************
  RIsearch v 1.2   --   RNA-RNA interaction search
  Copyright 2012 Anne Wenzel <wenzel@rth.dk> (RIsearch v.1.0 and v.1.1)
  Copyright 2021 Giulia I Corsi <giulia@rth.dk> (Extension of RIsearch v.1.1 in RIsearch v.1.2)

  This file is part of RIsearch.

  RIsearch is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  RIsearch is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with RIsearch, see file COPYING.
  If not, see <http://www.gnu.org/licenses/>.

***********************************************************/

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdarg.h>
#include <ctype.h>
#include "fasta.h"
#include <unistd.h>
#include "set_modified.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define NEGINF INT_MIN/2
#define GAP 5			/* position of '-' in alphabet, not as define if read from matrix... */

typedef struct
{
  int qbeg, qend, tbeg, tend;
  int max;
  char *ali_seq1, *ali_seq2, *ali_ia;
} IA;

int RIs_linSpace (unsigned char *qseq, unsigned char *tseq, int m, int n,
		  short dsm[6][6][6][6], int extensionpenalty, int threshold,
		  char *qname, char *tname, char *matname);
int reverse (char *s, int len);
int getMat (char *matname, short *bA_nu);
int getWeights (char *arrayname, float *weights);
void printMat (int **mat, int rows, int cols, unsigned char *seq1,
	       unsigned char *seq2);
int seq2ix (int len, char *seq, unsigned char *retIx, char *name, char *type);
unsigned char nt2index (char nt);
char index2nt (unsigned char ix);
void usage (char *progname);
void getArgs (int argc, char *argv[]);
int **allocIntMatrix (int rows, int cols);
float **allocFloatMatrix (int rows, int cols);
void freeIntMatrix (int **m, int rows);
void freeFloatMatrix (float **m, int rows);
void RIs (unsigned char *qseq, unsigned char *tseq, int m, int n,
	  short dsm[6][6][6][6], IA * hit);
void RIs_force_start_end_init (unsigned char *qseqIx, unsigned char *tseqIx,
			       int len_seq1, int len_seq2,
			       short dsm[6][6][6][6], char *matname);
void RIs_force_start_end_weighted (unsigned char *qseq, unsigned char *tseq,
				   int m, int n, float *weights,
				   short dsm[6][6][6][6], IA * hit,
				   char *matname);
int min2 (int a, int b);
int max3 (int a, int b, int c);
int max4 (int a, int b, int c, int d);
int max5 (int a, int b, int c, int d, int e);
float max3f (float a, float b, float c);
float max4f (float a, float b, float c, float d);
float max5f (float a, float b, float c, float d, float e);
float find_max_value_f (float **M_col, float **Ix_col, float **Iy_col, int *k,
			int *i, int j, int n, short dsm[6][6][6][6],
			unsigned char *qseq, unsigned char *tseq);
char *allocate_char_array (int length);
void set_alignment_symbols (char query_nt, char target_nt,
			    char *query_alignment, char *target_alignment);

/* values to be overwritten by command line parameters */
char *mat_name = "t04", *seq1file_name, *seq2file_name, *seq1_cli, *seq2_cli;
int tblen = 40;			/* trace-back length, that many nucleotides before 'maxHit' */
int extPen = 0;			/* extension penalty; used to compute dsm */
int force_start_val = -1;	/* values used to unitialize the first column of the M matrix. If sufficiently high, can force the interaction to start at position 0 of the DNA. */
int vicinity = 0;		/* to omit neighboring hits (subalignments) */
char printShort = 0;		/* switch p to print 1 line per IA, only pos&E, not IA itself */
/*TODO possibly several print styles, Vienna-like (one line, but still IA) */
int minScore = INT_MAX;		/* score cutoff; if not set, only print best */
double maxEnergy = INT_MAX;	/*energy cutoff, not even print 'best' if it's not lower than that! */
int doSubopt = 0;		/* flag : minScore in use */
int filterE = 0;		/* flag : filterE in use */
int weighted_positions = 0;	/*true if each position of interaction has a certain weight. */
char *pos_weights = "CRISPR_20nt_3p_5p";	/*name of the array being used to assign weights to positions */

int
main (int argc, char *argv[])
{
  unsigned long len_seq1, len_seq2;
  char *one, *two;
  unsigned char *qseqIx, *tseqIx;
  FASTAFILE *ffpQ, *ffpT;	/*for query/target respectively */
  char *nameQ, *nameT;
  short dsm[6][6][6][6];
  int check;

  getArgs (argc, argv);

  getMat (mat_name, &dsm[0][0][0][0]);

  if (seq2file_name)
    {				/* target given as file - or STDIN */
      ffpT = OpenFASTA (seq2file_name);
      if (ffpT == NULL)
	{
	  fprintf (stderr, "Target file %s is not readable\n", seq2file_name);
	  return -1;
	}
      while (ReadFASTA (ffpT, &two, &nameT, &len_seq2))
	{
/*can be done already when reading in first place */
	  tseqIx = malloc ((len_seq2) * sizeof *tseqIx);
	  check = seq2ix (len_seq2, two, tseqIx, nameT, "target");
	  if (check > 0)
	    len_seq2 -= check;	/*removed gap characters */
	  if (check < 0)
	    continue;		/*non-alpha char in input */
	  free (two);		/*free'ing space for full seq, as we have it as ix version */

	  if (seq1file_name)
	    {			/* query given as file */

	      ffpQ = OpenFASTA (seq1file_name);
	      if (ffpQ == NULL)
		{
		  fprintf (stderr, "Query file %s is not readable\n",
			   seq1file_name);
		  CloseFASTA (ffpT);
		  free (tseqIx);
		  free (nameT);
		  return -1;
		}
	      while (ReadFASTA (ffpQ, &one, &nameQ, &len_seq1))
		{

		  qseqIx = malloc ((len_seq1) * sizeof *qseqIx);
		  check = seq2ix (len_seq1, one, qseqIx, nameQ, "query");
		  if (check > 0)
		    len_seq1 -= check;	/*removed gap characters */
		  if (check < 0)
		    continue;	/*non-alpha char in input */
		  free (one);	/*free'ing space for full seq, as we have it as ix version */
		  if (printShort < 2)
		    printf
		      ("\n\nquery %s (%lu nts) vs. target %s (%lu nts)\n\n",
		       nameQ, len_seq1, nameT, len_seq2);
		  if (weighted_positions || (force_start_val >= 0))
		    {
		      if (force_start_val < 0)
			{
			  fprintf (stderr,
				   "Parameter -f must be set when using weights (-w).\n");
			  exit (1);
			}
		      if (!weighted_positions)
			{
			  fprintf (stderr,
				   "Parameter -w must be set when using force start (-f). Use array of weights \"noweigths\" to avoid this error.\n");
			  exit (1);
			}
		      if (extPen || tblen != 40 || doSubopt || filterE
			  || printShort || vicinity)
			{
			  fprintf (stderr,
				   "Options -d -s -n -l -e -p are not available in combination with options -f -w \n");
			  exit (1);
			}
		      RIs_force_start_end_init (qseqIx, tseqIx, len_seq1,
						len_seq2, dsm, mat_name);
		    }
		  else
		    {
		      RIs_linSpace (qseqIx, tseqIx, len_seq1, len_seq2, dsm,
				    extPen, minScore, nameQ, nameT, mat_name);
		    }
		  free (qseqIx);
		  free (nameQ);
		}
	      CloseFASTA (ffpQ);

	    }
	  else if (seq1_cli)
	    {			/* query given as command line parameter */
	      len_seq1 = strlen (seq1_cli);
	      qseqIx = malloc ((len_seq1) * sizeof *qseqIx);
	      check =
		seq2ix (len_seq1, seq1_cli, qseqIx, "from command line",
			"query");
	      if (check > 0)
		len_seq1 -= check;	/*removed gap characters */
	      if (check < 0)
		return -1;	/*non-alpha char in input -- break would loop through all query seqs, no use */
	      if (printShort < 2)
		printf
		  ("\n\nquery from_cli (%lu nts) vs. target %s (%lu nts)\n\n",
		   len_seq1, nameT, len_seq2);
	      if (weighted_positions || (force_start_val >= 0))
		{
		  if (force_start_val < 0)
		    {
		      fprintf (stderr,
			       "Parameter -f must be set when using weights (-w).\n");
		      exit (1);
		    }
		  if (!weighted_positions)
		    {
		      fprintf (stderr,
			       "Parameter -w must be set when using force start (-f). Use array of weights \"noweigths\" to avoid this error.\n");
		      exit (1);
		    }
		  if (extPen || tblen != 40 || doSubopt || filterE
		      || printShort || vicinity)
		    {
		      fprintf (stderr,
			       "Options -d -s -n -l -e -p are not available in combination with options -f -w \n");
		      exit (1);
		    }
		  RIs_force_start_end_init (qseqIx, tseqIx, len_seq1,
					    len_seq2, dsm, mat_name);
		}
	      else
		{
		  RIs_linSpace (qseqIx, tseqIx, len_seq1, len_seq2, dsm,
				extPen, minScore, "from_cli", nameT,
				mat_name);
		}
	      free (qseqIx);

	    }
	  else
	    {
	      fprintf (stderr, "No query seq given!");
	      /* is caught in getArg already -- alternative run seq against itself!? */
	    }

	  free (tseqIx);
	  free (nameT);
	}
      CloseFASTA (ffpT);

    }
  else if (seq2_cli)
    {				/*target given as command line parameter */

      len_seq2 = strlen (seq2_cli);
      tseqIx = malloc ((len_seq2) * sizeof *tseqIx);
      check =
	seq2ix (len_seq2, seq2_cli, tseqIx, "from command line", "target");
      if (check > 0)
	len_seq2 -= check;	/*removed gap characters */
      if (check < 0)
	return -1;		/*non-alpha char in input */

      if (seq1file_name)
	{			/* query given as file */

	  ffpQ = OpenFASTA (seq1file_name);
	  if (ffpQ == NULL)
	    {
	      fprintf (stderr, "Query file %s is not readable\n",
		       seq1file_name);
	      free (tseqIx);
	      return -1;
	    }
	  while (ReadFASTA (ffpQ, &one, &nameQ, &len_seq1))
	    {

	      qseqIx = malloc ((len_seq1) * sizeof *qseqIx);
	      check = seq2ix (len_seq1, one, qseqIx, nameQ, "query");
	      if (check > 0)
		len_seq1 -= check;	/*removed gap characters */
	      if (check < 0)
		continue;	/*non-alpha char in input */
	      free (one);

	      if (printShort < 2)
		printf
		  ("\n\nquery %s (%lu nts) vs. target from_cli (%lu nts)\n\n",
		   nameQ, len_seq1, len_seq2);
	      if (weighted_positions || (force_start_val >= 0))
		{
		  if (force_start_val < 0)
		    {
		      fprintf (stderr,
			       "Parameter -f must be set when using weights (-w).\n");
		      exit (1);
		    }
		  if (!weighted_positions)
		    {
		      fprintf (stderr,
			       "Parameter -w must be set when using force start (-f). Use array of weights \"noweigths\" to avoid this error.\n");
		      exit (1);
		    }
		  if (extPen || tblen != 40 || doSubopt || filterE
		      || printShort || vicinity)
		    {
		      fprintf (stderr,
			       "Options -d -s -n -l -e -p are not available in combination with options -f -w \n");
		      exit (1);
		    }
		  RIs_force_start_end_init (qseqIx, tseqIx, len_seq1,
					    len_seq2, dsm, mat_name);
		}
	      else
		{
		  RIs_linSpace (qseqIx, tseqIx, len_seq1, len_seq2, dsm,
				extPen, minScore, nameQ, "from_cli",
				mat_name);
		}
	      free (qseqIx);
	      free (nameQ);
	    }
	  CloseFASTA (ffpQ);

	}
      else if (seq1_cli)
	{			/* query given as command line parameter */

	  len_seq1 = strlen (seq1_cli);
	  qseqIx = malloc ((len_seq1) * sizeof *qseqIx);
	  check =
	    seq2ix (len_seq1, seq1_cli, qseqIx, "from command line", "query");
	  if (check > 0)
	    len_seq1 -= check;	/*removed gap characters */
	  if (check < 0)
	    return -1;		/* non-alpha char in input -- break would loop through queries, no use */
	  if (printShort < 2)
	    printf
	      ("\n\nquery from_cli (%lu nts) vs. target from_cli (%lu nts)\n\n",
	       len_seq1, len_seq2);
	  if (weighted_positions || (force_start_val >= 0))
	    {
	      if (force_start_val < 0)
		{
		  fprintf (stderr,
			   "Parameter -f must be set when using weights (-w).\n");
		  exit (1);
		}
	      if (!weighted_positions)
		{
		  fprintf (stderr,
			   "Parameter -w must be set when using force start (-f). Use array of weights \"noweigths\" to avoid this error.\n");
		  exit (1);
		}
	      if (extPen || tblen != 40 || doSubopt || filterE || printShort
		  || vicinity)
		{
		  fprintf (stderr,
			   "Options -d -s -n -l -e -p are not available in combination with options -f -w \n");
		  exit (1);
		}
	      RIs_force_start_end_init (qseqIx, tseqIx, len_seq1, len_seq2,
					dsm, mat_name);
	    }
	  else
	    {
	      RIs_linSpace (qseqIx, tseqIx, len_seq1, len_seq2, dsm, extPen,
			    minScore, "from_cli", "from_cli", mat_name);
	    }
	  free (qseqIx);
	}
      else
	{
	  fprintf (stderr, "No query seq given!");
	  /* is caught in getArg already -- alternative run seq against itself!? */
	}

      free (tseqIx);

    }
  else
    {
      fprintf (stderr, "No target seq given!");
      /* is caught in getArg already -- alternative run seq against itself!? */
    }

  return 0;
}

/* can save first traversal as length is known before! */
int
reverse (char *str_beg, int j)
{
  char tmp;
  char *str_end = &str_beg[j];

  while (str_end > str_beg)
    {
      tmp = *str_beg;
      *str_beg++ = *str_end;
      *str_end-- = tmp;
    }

  return 0;
}

void
printMat (int **mat, int rows, int cols, unsigned char *seq1,
	  unsigned char *seq2)
{
  int i, j;
  printf ("\t-");
  for (j = 0; j < cols - 1; j++)
    {
      printf ("\t%c", index2nt (*(seq2 + j)));
    }
  for (i = 0; i < rows; i++)
    {
      printf ("\n%c", (i == 0 ? '-' : index2nt (*(seq1 + i - 1))));
      for (j = 0; j < cols; j++)
	{
	  printf ("\t%d", (mat[i][j] == NEGINF ? -8 : mat[i][j]));	/* -8 as dummy for -inf */
	}
    }
  printf ("\n\n");
}

void
printfloatMat (float **mat, int rows, int cols, unsigned char *seq1,
	       unsigned char *seq2)
{
  int i, j;
  printf ("\t-");
  for (j = 0; j < cols - 1; j++)
    {
      printf ("\t%c", index2nt (*(seq2 + j)));
    }
  for (i = 0; i < rows; i++)
    {
      printf ("\n%c", (i == 0 ? '-' : index2nt (*(seq1 + i - 1))));
      for (j = 0; j < cols; j++)
	{
	  printf ("\t%f", (mat[i][j] == NEGINF ? -8 : mat[i][j]));	/* -8 as dummy for -inf */
	}
    }
  printf ("\n\n");
}

int
seq2ix (int len, char *seq, unsigned char *retIx, char *name, char *type)
{
  int i;
  int gapcnt = 0;
  for (i = 0; i < len; i++)
    {
      switch (seq[i])
	{
	case 'A':
	case 'a':
	  *(retIx + i - gapcnt) = 0;
	  break;
	case 'C':
	case 'c':
	  *(retIx + i - gapcnt) = 1;
	  break;
	case 'G':
	case 'g':
	  *(retIx + i - gapcnt) = 2;
	  break;
	case 'T':
	case 't':
	case 'U':
	case 'u':
	  *(retIx + i - gapcnt) = 3;
	  break;
	case 'N':
	case 'n':
	  *(retIx + i - gapcnt) = 4;
	  break;
	case '-':
	case '.':		/*discard gaps from input  --  also add "case ' ' :"??? */
	  gapcnt++;
	  break;
	default:
	  if (isalpha (seq[i]))
	    {
	      fprintf (stderr,
		       "Nonstandard nucleotide code '%c' in %s sequence '%s'. Replaced with 'N'\n",
		       seq[i], type, name);
	      *(retIx + i - gapcnt) = 4;
	      break;
	    }
	  else
	    {			/*skip sequence!? */
	      fprintf (stderr,
		       "Unexpected character '%c' in %s sequence '%s'. Skipping sequence.\n",
		       seq[i], type, name);
	      return -1;
	    }
	}
    }
  return gapcnt;
}

unsigned char
nt2index (char nt)
{
  switch (nt)
    {
    case 'A':
    case 'a':
      return 0;
    case 'C':
    case 'c':
      return 1;
    case 'G':
    case 'g':
      return 2;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      return 3;
    case 'N':
    case 'n':
      return 4;
    case '-':
      return 5;			/*any case!? */
    default:
      fprintf (stderr, "Nonstandard nucleotide code: %c\n", nt);
      return 4;
    }
}

char
index2nt (unsigned char ix)
{
  switch (ix)
    {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'U';
    case 4:
      return 'N';
    case 5:
      return '-';		/*any case!? */
    default:
      fprintf (stderr, "\nUnknown symbol >>%d<< found.\n", ix);
      exit (1);			/*just skip it!? not fail? */
    }
}

void
usage (char *progname)
{
  fprintf (stderr,
	   "====== RIsearch1 ver 1.2 ======\n= RNA-RNA interaction search =\n");
  fprintf (stderr,
	   "=   Contact: wenzel@rth.dk   =\n==============================\n\n");
  fprintf (stderr, "Usage: \t%s [ARGUMENTS]\n", progname);
  fprintf (stderr, "\n   [INPUT]\n");
  fprintf (stderr, "\t-q <file> Fasta file containing query sequence(s)\n");
  fprintf (stderr,
	   "\t-t <file> Fasta file containing target sequence(s) or specify '-t -' to pass the fasta formatted input to STDIN\n");
  fprintf (stderr, "\t-Q <str>  Query sequence only, direct as string\n");
  fprintf (stderr, "\t-T <str>  Target sequence as string\n");
  fprintf (stderr,
	   "   If q is given, Q will be ignored; same for t over T.\n");
  fprintf (stderr, "\n   [OPTIONS]\n");
  fprintf (stderr,
	   "\t-d <int>  per-nucleotide extension penalty given in dacal/mol (recommended: 30, default: 0)\n");
  fprintf (stderr,
	   "\t-s <int>  threshold for suboptimal duplexes (minimum score, NOT energy)\n");
  fprintf (stderr,
	   "\t-n <int>  neighborhood (only backtrack from best position within this range - to omit many overlapping results), default is 0, backtrack all\n");
  fprintf (stderr,
	   "\t-f <int> force the interaction to start and end, respectively, at the 3'and 5' end of target and end at the query's 3' end. Takes input int >0.\n");
  fprintf (stderr,
	   "\t            Use a high value, e.g. 200*max(length(query),length(target))\n");
  fprintf (stderr,
	   "\t            It is currently not possible to compute weighted interactions without a fixed start and end.\n");
  fprintf (stderr, "\t            This option requires -w.\n");
  fprintf (stderr, "\t-l <int>  max trace back length (default: 40)\n");
  fprintf (stderr,
	   "\t-m <str>  matrix to use, t99 or t04(def) or su95 or su95_noGU or sl04_noGU or modified e1 e2\n");
  fprintf (stderr,
	   "\t-w <str>  weights vector to use, CRISPR_20nt_5p_3p or noweights. \n");
  fprintf (stderr,
	   "\t            Weights length mush be >= than the length of the query -1.\n");
  fprintf (stderr,
	   "\t            This option requires -f. Please set -f appropriately, do not try to use a low -f value to avoid the force start. \n");
  fprintf (stderr,
	   "\t            Launched with this option, RIsearch does not perform step 1 (memory optimization) of the algorithm (see Wenzel et al. 2012).\n");
  fprintf (stderr, "\t-e <num>  energy threshold, checked after backtrack\n");
  fprintf (stderr,
	   "\t            An interaction is only printed if the predicted hybridization energy is lower than (or equal to) this threshold.\n");
  fprintf (stderr,
	   "\t            Also the 'best hit' per query/target pair might be filtered out, appears only once if at all (and not necessarily as first).\n");
  fprintf (stderr,
	   "\t-p          switch for short output, for backwards compatibility, same as -p1\n");
  fprintf (stderr, "\t-p[1-3]     different shorter output modes:\n");
  fprintf (stderr,
	   "\t\t-p1     one line per hit, incl. interaction string, still header for each pair (query / target)\n");
  fprintf (stderr,
	   "\t\t-p2     one line per hit, tab seperated 'Qname Qbeg Qend Tname Tbeg Tend score energy'; no header; 'best' hit only once, not first\n");
  fprintf (stderr,
	   "\t\t-p3     one line per pair with number of hits that would have been printed, tab seperated 'Qname Tname hit-count'\n");
  fprintf (stderr,
	   "\t\t        without p (and with p1), the 'best' interaction (per pair) is always shown first, and repeated in the list of results\n");
  fprintf (stderr, "\n\n");
  exit (1);
}

void
getArgs (int argc, char *argv[])
{
  char c;
  while ((c = getopt (argc, argv, "q:t:Q:T:d:X:m:s:e:n:w:l:f:p::")) != -1)
    switch (c)
      {
      case 'q':
	seq1file_name = optarg;
	break;
      case 't':
	seq2file_name = optarg;
	break;
      case 'Q':
	seq1_cli = optarg;
	break;
      case 'T':
	seq2_cli = optarg;
	break;
      case 'd':
	extPen = atoi (optarg);
	break;
      case 'X':
	break;			/*silent var */
      case 'm':
	mat_name = optarg;
	break;
      case 's':
	minScore = atoi (optarg);
	doSubopt = 1;
	break;
      case 'e':
	maxEnergy = atof (optarg);
	filterE = 1;
	break;
      case 'n':
	vicinity = atoi (optarg);
	break;
      case 'w':
	weighted_positions = 1;
	pos_weights = optarg;
	break;
      case 'l':
	tblen = atoi (optarg);
	break;
      case 'f':
	force_start_val = atoi (optarg);
	break;
      case 'p':
	if (optarg)
	  printShort = atoi (optarg);
	else
	  printShort = 1;
	break;
      case '?':
	usage (argv[0]);
	break;
      default:
	abort ();
      }

  if (!((seq1file_name || seq1_cli) && (seq2file_name || seq2_cli)))
    {
      fprintf (stderr,
	       "\nYou need to provide a query (see -Q or -q option) and a target (-T/-t)\n\n");
      usage (argv[0]);
    }
  if (seq1file_name && (!strcmp (seq1file_name, "-")))
    {
      fprintf (stderr,
	       "\nQuery can currently not be read from STDIN, only target can!\n\n");
      usage (argv[0]);
    }
}

int
getMat (char *matname, short *bA_nu)
{
  short *bA_bas, *bA_ext;
  int i;
  extern short dsm_extend[6][6][6][6];
  bA_ext = &dsm_extend[0][0][0][0];
  double e1,e2;  
  if(sscanf(matname, "modified %lf %lf", &e1, &e2) == 2){
  	set_modified();
  	extern short modified_dsm[6][6][6][6];
	extern short modified_dsm_extend[6][6][6][6];
  	bA_ext = &modified_dsm_extend[0][0][0][0];
  	bA_bas = &modified_dsm[0][0][0][0];
  }
  else if (!strcmp (matname, "t04"))
    {
      extern short dsm_t04[6][6][6][6];
      bA_bas = &dsm_t04[0][0][0][0];
    }
  else if (!strcmp (matname, "t99"))
    {
      extern short dsm_t99[6][6][6][6];
      bA_bas = &dsm_t99[0][0][0][0];
    }
  else if (!strcmp (matname, "su95"))
    {
      extern short dsm_su95_rev_wGU_pos[6][6][6][6];
      bA_bas = &dsm_su95_rev_wGU_pos[0][0][0][0];
    }
  else if (!strcmp (matname, "sl04_noGU"))
    {
      extern short dsm_slh04_woGU_pos[6][6][6][6];
      bA_bas = &dsm_slh04_woGU_pos[0][0][0][0];
    }
  else if (!strcmp (matname, "su95_noGU"))
    {
      extern short dsm_su95_rev_woGU_pos[6][6][6][6];
      bA_bas = &dsm_su95_rev_woGU_pos[0][0][0][0];
    }
  else
    {
      fprintf (stderr,
	       "Undefined matrix, -m needs to be set to either t99 or t04 for RNA-RNA interaction, su95 or su95_noGU for RNA-DNA interaction or sl04_noGU for DNA interaction or modified e1 e2 \n");
      exit (1);
    }


  /* create dsm from   dsm_base - d * dsm_extend   */
  for (i = 0; i < 1296; i++)
    {
      *(bA_nu + i) = *(bA_bas + i) - extPen * *(bA_ext + i);	/* bA_nu[i] =  ... also works */
    }
  return 0;
}

int
max3 (int a, int b, int c)
{
  int max = a;
  if (b > max)
    {
      max = b;
    }
  if (c > max)
    {
      max = c;
    }
  return max;
}

int
min2 (int a, int b)
{
  return b > a ? a : b;
}

int
max4 (int a, int b, int c, int d)
{
  int max = a;
  if (b > max)
    {
      max = b;
    }
  if (c > max)
    {
      max = c;
    }
  if (d > max)
    {
      max = d;
    }
  return max;
}

int
max5 (int a, int b, int c, int d, int e)
{
  int max = a;
  if (b > max)
    {
      max = b;
    }
  if (c > max)
    {
      max = c;
    }
  if (d > max)
    {
      max = d;
    }
  if (e > max)
    {
      max = e;
    }
  return max;
}

float
max3f (float a, float b, float c)
{
  float max = a;
  if (b > max)
    {
      max = b;
    }
  if (c > max)
    {
      max = c;
    }
  return max;
}

float
max4f (float a, float b, float c, float d)
{
  float max = a;
  if (b > max)
    {
      max = b;
    }
  if (c > max)
    {
      max = c;
    }
  if (d > max)
    {
      max = d;
    }
  return max;
}

float
max5f (float a, float b, float c, float d, float e)
{
  float max = a;
  if (b > max)
    {
      max = b;
    }
  if (c > max)
    {
      max = c;
    }
  if (d > max)
    {
      max = d;
    }
  if (e > max)
    {
      max = e;
    }
  return max;
}

float
find_max_value_f (float **M, float **Ix, float **Iy, int *k, int *i, int j,
		  int n, short dsm[6][6][6][6], unsigned char *qseq,
		  unsigned char *tseq)
{
  float max = 0.0;
  *i = 0;			/*row in which maxval is found */
  *k = 0;			/*matrix in which the maximum has been found; 0 1 2 if M Ix Iy */
  if (M[n][j] + dsm[qseq[n - 1]][GAP][tseq[j - 1]][GAP] > max)
    {
      max = M[n][j] + dsm[qseq[n - 1]][GAP][tseq[j - 1]][GAP];
      *i = n;
    }
  /* allow to end with a gap */
  if (Ix[n][j] + dsm[qseq[n - 1]][GAP][GAP][GAP] > max)
    {
      max = Ix[n][j] + dsm[qseq[n - 1]][GAP][GAP][GAP];
      *k = 1;
      *i = n;
    }
  if (Iy[n][j] + dsm[GAP][GAP][tseq[n - 1]][GAP] > max)
    {
      max = Iy[n][j] + dsm[GAP][GAP][tseq[n - 1]][GAP];
      *k = 2;
      *i = n;
    }
  return max;
}

void
set_alignment_symbols (char query_nt, char target_nt, char *query_alignment,
		       char *target_alignment)
{
/*this function takes two nucleotides as inputs and set characters in the alignment string to:
 | : if the two nt can base pair,
 M : if there is a mismatch,
 W : if there is a wobble base pair*/
  if ('A' == query_nt || 'A' == target_nt)
    {
      if ('U' == query_nt || 'U' == target_nt)
	{
	  *query_alignment = '|';
	  *target_alignment = '|';
	}
      else
	{
	  *query_alignment = 'M';
	  *target_alignment = 'M';
	}
    }
  else if ('G' == query_nt || 'G' == target_nt)
    {
      if ('U' == query_nt || 'U' == target_nt)
	{
	  *query_alignment = 'W';
	  *target_alignment = 'W';
	}
      else if ('C' == query_nt || 'C' == target_nt)
	{
	  *query_alignment = '|';
	  *target_alignment = '|';
	}
      else
	{
	  *query_alignment = 'M';
	  *target_alignment = 'M';
	}
    }
  else
    {
      *query_alignment = 'M';
      *target_alignment = 'M';
    }

}

int **
allocIntMatrix (int rows, int cols)
{
  /* function to allocate an Integer Matrix of size rows x cols */
  int **m;
  int i;

  m = malloc (rows * sizeof (int *));

  if (!m)
    {
      printf ("Cannot allocate integer matrix with %d rows\n", rows);
      exit (1);
    }

  for (i = 0; i < rows; i++)
    {
      m[i] = malloc (cols * sizeof (int));
      if (!m[i])
	{
	  printf
	    ("Cannot allocate column %d of matrix with %d rows and %d cols\n",
	     i, rows, cols);
	  exit (1);
	}
    }
  return m;
}

float **
allocFloatMatrix (int rows, int cols)
{
  /* function to allocate an Integer Matrix of size rows x cols */
  float **m;
  int i;

  m = malloc (rows * sizeof (float *));

  if (!m)
    {
      printf ("Cannot allocate integer matrix with %d rows\n", rows);
      exit (1);
    }

  for (i = 0; i < rows; i++)
    {
      m[i] = malloc (cols * sizeof (float));
      if (!m[i])
	{
	  printf
	    ("Cannot allocate column %d of matrix with %d rows and %d cols\n",
	     i, rows, cols);
	  exit (1);
	}
    }
  return m;
}

char *
allocate_char_array (int length)
{
  char *string;
  int i;

  string = malloc (length + 1 * sizeof (char));
  for (i = 0; i < length; i++)
    {
      string[i] = 'X';
    }
  string[length] = '\0';
  return string;
}

void
freeIntMatrix (int **m, int rows)
{
  while (rows--)
    free ((char *) (m[rows]));
  free ((char *) (m));
}

void
freeFloatMatrix (float **m, int rows)
{
  while (rows--)
    free ((char *) (m[rows]));
  free ((char *) (m));
}

void
RIs (unsigned char *qseq,	/* query sequence - numeric representation */
     unsigned char *tseq,	/* target sequence - reversed */
     int m,			/* query seq length */
     int n,			/* target seq length */
     short dsm[6][6][6][6],	/* scoring matrix */
     IA * hit			/* pointer to struct, fill results */
  )
{
  int **M, **Ix, **Iy;		/* matrices for alignment scores ending in different states */
  int maxi, maxj;		/* k is 0,1,2 for M,Ix,Iy !? - max will never be found in gapped anyway!? */
  int maxval, maxk, tmp;	/* maxk not needed!? k itself could even be char... */
  int i, j, k;
  int l = 0;			/* alilen so far -used in backtrack */
  int mVal, xVal, yVal, nVal;	/* values coming from M, Ix, Iy, or starting a NEW alignment */

  M = allocIntMatrix (m + 1, n + 1);	/* (Mis)Match */
  Ix = allocIntMatrix (m + 1, n + 1);	/* Insertion in x(=query), so x paired to gap (in y) */
  Iy = allocIntMatrix (m + 1, n + 1);	/* Insertion(=bulge) in y(=target) */
  maxi = maxj = maxk = 0;

  M[0][0] = Ix[0][0] = Iy[0][0] = 0;

  /*init first row (j=0) -- this is COL - change throughout!? */

  /* explicitly handling of the boundary condition since i-2 is not defined for i = 1
     NOT needed anymore, is now just 0 anyway
     Ix[1][0] = 0; //MAX(0, dsm[GAP][qseq[0]][GAP][GAP]); 
     Iy[1][0] = M[1][0] = NEGINF; // not possible to occure 
   */
  for (i = 1; i <= m; i++)
    {
      Iy[i][0] = M[i][0] = NEGINF;	/* not possible before beginning of target seq */
      Ix[i][0] = 0;		/*MAX(0, dsm[qseq[i-2]][qseq[i-1]][GAP][GAP]); *//* require to always have a match first! */
    }

  /*init first col (i=0) --- this is ROW - change throughout!?
     Iy[0][1] = 0; // MAX(0, dsm[GAP][GAP][GAP][tseq[0]]);
     Ix[0][1] = M[0][1] = NEGINF; 
   */
  for (j = 1; j <= n; j++)
    {
      Ix[0][j] = M[0][j] = NEGINF;
      Iy[0][j] = 0;		/*MAX(0, dsm[GAP][GAP][tseq[j-2]][tseq[j-1]]); */
    }

  /* The initialization of i=1 column and j=1 row have to be handled explicitly
     since at this point we do not have two residues to use.

     Handle the (1,1) cell explicitly since the boundary recursion includes (i-2) or (j-2) cases.
   */

  M[1][1] = dsm[GAP][qseq[0]][GAP][tseq[0]];	/*MAX(0,dsm[GAP][qseq[0]][GAP][tseq[0]]); */
  maxval = M[1][1] + dsm[qseq[0]][GAP][tseq[0]][GAP];	/* MAX(0,dsm[qseq[0]][GAP][tseq[0]][GAP]); */
  /* (1,1) cell can not be in Ix or Iy state. */
  Ix[1][1] = Iy[1][1] = NEGINF;

  /* init j=1 row */
  for (i = 2; i <= m; i++)
    {
      /* value for M matrix, case we have a pair here (k=0) */
      /* first letter of target sequence, so it can ONLY come from gapped */
      /* if previous cell has not been 0, than add current pair */
      /* if the case is removed on top, it can be here as well, save checks...

         xVal = Ix[i-1][0] != 0 ? Ix[i-1][0] + dsm[qseq[i-2]][qseq[i-1]][GAP][tseq[0]] : -1; // reflects (i-1,i;-,1)  coming from Ix (seq1 paired to gap)
         otherwise - could start NEW ali here OR 0 if nothing else scores positive
         nVal = dsm[GAP][qseq[i-1]][GAP][tseq[0]];
         M[i][1] = max3(xVal,nVal,0); */
      M[i][1] = dsm[GAP][qseq[i - 1]][GAP][tseq[0]];	/* MAX(0, ); */

      tmp = M[i][1] + dsm[qseq[i - 1]][GAP][tseq[0]][GAP];
      if (tmp > maxval)
	{
	  maxval = tmp;
	  maxi = i;
	  maxj = 1;
	  maxk = 0;
	}

      /* value for Ix matrix, case query sequence paired to gap (k=1) */
      /* prev. match, now gap - add (Xi-1, Xi; Y1, -) */
      mVal =
	M[i - 1][1] !=
	0 ? M[i - 1][1] + dsm[qseq[i - 2]][qseq[i - 1]][tseq[0]][GAP] : -1;
      /* extending existing gap  - add (Xi-1, Xi; -, -) */
      xVal =
	Ix[i - 1][1] !=
	0 ? Ix[i - 1][1] + dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] : -1;
      /* OR start new alignment that starts in gap, reflected by (-, Xi; -, -)  */
      /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */
      /*Ix[i][1] = max4(mVal,xVal,nVal,0); */
      Ix[i][1] = max3 (mVal, xVal, 0);	/*removed nVal option */
/* do not allow alignment ending in Xi!
    tmp = Ix[i][1] + dsm[qseq[i-1]][GAP][GAP][GAP];   // (Xi, -; -, -)
    if (tmp > maxval) {
      maxval = tmp;
      maxi = i; maxj = 1; maxk = 1;
    }
*/
      /* value for Iy matrix, case target sequence paired to gap (k=2) */
      /* not possible in this row */
      Iy[i][1] = NEGINF;

    }

  /* init i=1 column */
  for (j = 2; j <= n; j++)
    {
      /* value for M matrix, case we have a pair here (k=0) */
      /* coming from gap in qseq (Iy-matrix), adding (-, X1; Yj-1, Yj) */
      /* not possible, set to 0 before for all of them! */
      /* yVal = Iy[0][j-1] != 0 ? Iy[0][j-1] + dsm[GAP][qseq[0]][tseq[j-2]][tseq[j-1]] : -1; */
      /* OR starting a new alignment with (-, X1; -, Yj) */
      /* nVal = dsm[GAP][qseq[0]][GAP][tseq[j-1]]; */
      /* coming from gap in tseq NOT possible as we're looking at first position of query!? */
      /* M[1][j] = max3(yVal,nVal,0); */
      M[1][j] = dsm[GAP][qseq[0]][GAP][tseq[j - 1]];	/* MAX(0, ); */

      tmp = M[1][j] + dsm[qseq[0]][GAP][tseq[j - 1]][GAP];
      if (tmp > maxval)
	{
	  maxval = tmp;
	  maxi = 1;
	  maxj = j;
	  maxk = 0;
	}

      /* value for Ix matrix, case query sequence paired to gap (k=1) */
      /* not possible in this column */
      Ix[1][j] = NEGINF;

      /* value for Iy matrix, case target sequence paired to gap (k=2) */
      /*coming from match, opening gap - add (X1, -; Yj-1, Yj) */
      mVal =
	M[1][j - 1] !=
	0 ? M[1][j - 1] + dsm[qseq[0]][GAP][tseq[j - 2]][tseq[j - 1]] : -1;
      /* extending an existing gap - add (-, -; Yj-1, Yj) */
      yVal =
	Iy[1][j - 1] !=
	0 ? Iy[1][j - 1] + dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] : -1;
      /* or start new ali, begin with gap -- forbidden ! */
      /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */
      Iy[1][j] = max3 (mVal, yVal, 0);	/*removed option nVal */
/* do not allow alignments ending in gap 
    tmp = Iy[1][j] + dsm[GAP][GAP][tseq[j-1]][GAP];
    if (tmp > maxval) {
      maxval = tmp;
      maxi = 1; maxj = j; maxk = 2;
    }
*/
    }
  /* initialization of first rows and columns completed
     recursion to complete alignment with two residues follows
   */
  for (i = 2; i <= m; i++)
    {				/*alt bed: *p1  */
      for (j = 2; j <= n; j++)
	{

	  /* value for M matrix, case we have a pair here (k=0) */
	  /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
	  mVal =
	    M[i - 1][j - 1] !=
	    0 ? M[i - 1][j - 1] +
	    dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] : -1;
	  /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
	  xVal =
	    Ix[i - 1][j - 1] !=
	    0 ? Ix[i - 1][j - 1] +
	    dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq[j - 1]] : -1;
	  /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
	  yVal =
	    Iy[i - 1][j - 1] !=
	    0 ? Iy[i - 1][j - 1] +
	    dsm[GAP][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] : -1;
	  /* starting a new alignment with this pair: (-, Xi; -, Yj) */
	  nVal = dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]];

	  M[i][j] = max4 (mVal, xVal, yVal, nVal);
	  tmp = M[i][j] + dsm[qseq[i - 1]][GAP][tseq[j - 1]][GAP];
	  if (tmp > maxval)
	    {
	      maxval = tmp;
	      maxi = i;
	      maxj = j;
	      maxk = 0;
	    }

	  /* value for Ix matrix, case query paired to gap (k=1) */
	  /*coming from match, add (Xi-1, Xi; Yj, -) */
	  mVal =
	    M[i - 1][j] !=
	    0 ? M[i - 1][j] +
	    dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 1]][GAP] : -1;
	  /*extend existing gap, add (Xi-1, Xi; -, -) */
	  xVal =
	    Ix[i - 1][j] !=
	    0 ? Ix[i - 1][j] + dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] : -1;
	  /* start new alignment - starts with gap -> ILLEGAL! */
	  /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */

	  Ix[i][j] = max3 (mVal, xVal, 0);	/*removed option nVal */
/* do not allow alignments ending in gap
      tmp = Ix[i][j] + dsm[qseq[i-1]][GAP][GAP][GAP];
      if (tmp > maxval) {
        maxval = tmp;
        maxi = i; maxj = j; maxk = 1;
      }
*/
	  /* value for Iy matrix, case target paired to gap (k=2) */
	  /*coming from match, add (Xi, -; Yj-1, Yj) */
	  mVal =
	    M[i][j - 1] !=
	    0 ? M[i][j - 1] +
	    dsm[qseq[i - 1]][GAP][tseq[j - 2]][tseq[j - 1]] : -1;
	  /*extend existing gap, add (-, -; Yj-1, Yj) */
	  yVal =
	    Iy[i][j - 1] !=
	    0 ? Iy[i][j - 1] + dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] : -1;
	  /* start new alignment - starts with gap LEGAL!? */
	  /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */

	  Iy[i][j] = max3 (mVal, yVal, 0);	/*removed option nVal */
/* do not allow alignments ending in gap
      tmp = Iy[i][j] + dsm[GAP][GAP][tseq[j-1]][GAP];
      if (tmp > maxval) {
        maxval = tmp;
        maxi = i; maxj = j; maxk = 2;
      }
*/
	}
    }
#ifdef VERBOSE
  printf ("found maxval %d on pos %d/%d in mat %d\n", maxval, maxi, maxj, maxk);	/* pos are 1-based */
  printf ("print F[0] matrix:\n");
  printMat (M, m + 1, n + 1, qseq, tseq);
  printf ("print F[1] matrix:\n");
  printMat (Ix, m + 1, n + 1, qseq, tseq);
  printf ("print F[2] matrix:\n");
  printMat (Iy, m + 1, n + 1, qseq, tseq);
#endif

/*backtrack*/
  tmp = (int) (1.5 * tblen);

  i = maxi;
  j = maxj;
  k = maxk;			/* 0-1-2 M Ix Iy - should always be 0 to begin with */
  if (k != 0)
    {
      fprintf (stderr,
	       "\nErr: Found highest value in one of the gap matrices (k=%d)!?\n\n",
	       k);
    }

  tmp -= 2;
  while ((i > 0 && j > 0) && (M[i][j] > 0 || Ix[i][j] > 0 || Iy[i][j] > 0))
    {
      if (l > tmp)
	{
	  printf
	    ("Interaction longer than max, so the following is only the end of the full alignment:\n");
	  /*alt: stop here / reallocate? prevent creation of longer alignments in first place? */
	  break;
	}
#ifdef DEBUG
#if VERBOSE>1
      fprintf (stderr, "in backtrack step with i=%d, j=%d, k=%d\n", i, j, k);
#endif
#endif

/*l++ in end of while instead of every sub? */

      /* find which cell gave score */
      if (k == 0)
	{			/* highest score in M matrix, having a match! */
#ifdef DEBUG
#if VERBOSE>1
	  /*printf("check next if M[i][j] (%d) ==  M[i-1][j-1] (%d) + dsm[qseq[i-2]][qseq[i-1]][tseq[j-2]][tseq[j-1]] (%d) \n", M[i][j], M[i-1][j-1], dsm[qseq[i-2]][qseq[i-1]][tseq[j-2]][tseq[j-1]] );  //this line will segfault when at i=1 OR j=1 */
	  printf ("check next if M[i][j] (%d) \n", M[i][j]);
	  printf (" == dsm[GAP][qseq[i-1]][GAP][tseq[j-1]] (%d)\n",
		  dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]]);
#endif
#endif
	  if (M[i][j] == dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]])
	    {
/*      started new alignment */
#ifdef DEBUG
#if VERBOSE>1
	      printf ("means we started alignment at M[%d][%d] \n", i, j);
#endif
#endif
	      k = 3;
	      /* need to print here already OR later check k and break */
	      hit->ali_seq1[l] = index2nt (qseq[--i]);
	      hit->ali_seq2[l] = index2nt (tseq[--j]);
	      if (qseq[i] + tseq[j] == 3)
		{		/* AU or CG pair */
		  hit->ali_ia[l] = '|';
		}
	      else if (qseq[i] + tseq[j] == 5)
		{		/*GU wobble pair */
		  hit->ali_ia[l] = '.';
		}
	      else
		{
		  hit->ali_ia[l] = ' ';
		}
	      l++;
	      break;
	    }
	  else if (M[i][j] ==
		   M[i - 1][j - 1] +
		   dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]])
	    {
	      /* access is save here as matrix always build from subset only, so we need to come to an end with check before!? */
/*      previous was also match */
#ifdef DEBUG
#if VERBOSE>1
	      printf ("Not the case, so check if it's \n");
	      printf (" ==  M[i-1][j-1] (%d)\n", M[i - 1][j - 1]);
	      printf
		(" + dsm[qseq[i-2]][qseq[i-1]][tseq[j-2]][tseq[j-1]] (%d) \n",
		 dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]]);
#endif
#endif
	      k = 0;
	    }
	  else if (M[i][j] ==
		   Ix[i - 1][j - 1] +
		   dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq[j - 1]])
	    {
/*      coming from gap in seq2/target */
	      k = 1;
	    }
	  else if (M[i][j] ==
		   Iy[i - 1][j - 1] +
		   dsm[GAP][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]])
	    {
/*      coming from gap in seq1/query */
	      k = 2;
	    }
	  else
	    {
	      printf ("unexpected value in k=0.\n");
	    }
	  hit->ali_seq1[l] = index2nt (qseq[--i]);
	  hit->ali_seq2[l] = index2nt (tseq[--j]);
/*      printf("test to calc %c (%d or %d) + %c (%d or %d) = %d or %d\n", hit.ali_seq1[l], hit.ali_seq1[l],qseq[i], hit.ali_seq2[l], hit.ali_seq2[l],tseq[j], hit.ali_seq1[l]+hit.ali_seq2[l], qseq[i]+tseq[j] ); */
	  if (qseq[i] + tseq[j] == 3)
	    {
	      hit->ali_ia[l] = '|';
	    }
	  else if (qseq[i] + tseq[j] == 5)
	    {
	      hit->ali_ia[l] = '.';
	    }
	  else
	    {
	      hit->ali_ia[l] = ' ';
	    }
	  l++;
#ifdef DEBUG
#if VERBOSE>1
	  printf ("print pos %d / %d = %c/ %c\n", i, j, hit->ali_seq1[l - 1], hit->ali_seq2[l - 1]);
#endif
#endif

	}
      else if (k == 1)
	{			/* seq1(query) paired to a gap (in target) */
	  if (Ix[i][j] ==
	      M[i - 1][j] + dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 1]][GAP])
	    {
	      k = 0;		/* open a new gap coming from match */
	    }
	  else if (Ix[i][j] ==
		   Ix[i - 1][j] + dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP])
	    {
	      k = 1;		/* extend existing gap */
	    }
	  else if (Ix[i][j] == dsm[GAP][qseq[i - 1]][GAP][GAP])
	    {
	      k = 3;		/* start new alignment with gap; not possible, prevented by scoring... */
	      fprintf (stderr,
		       "\nErr: This alignment starts in a gap - not even an option!?\n");
	      hit->ali_seq1[l] = index2nt (qseq[--i]);
	      hit->ali_ia[l] = ' ';
	      hit->ali_seq2[l++] = '-';
	      break;
	    }
	  else
	    {
	      printf ("unexpected case in k=1 : %d\n", Ix[i][j]);
	    }
	  hit->ali_seq1[l] = index2nt (qseq[--i]);
	  hit->ali_ia[l] = ' ';
	  hit->ali_seq2[l++] = '-';

	}
      else if (k == 2)
	{			/* seq2(target) paired to a gap (in query) */
	  if (Iy[i][j] ==
	      M[i][j - 1] + dsm[qseq[i - 1]][GAP][tseq[j - 2]][tseq[j - 1]])
	    {
	      k = 0;		/* open a new gap coming from match */
	    }
	  else if (Iy[i][j] ==
		   Iy[i][j - 1] + dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]])
	    {
	      k = 2;		/* extend existing gap */
	    }
	  else if (Iy[i][j] == dsm[GAP][GAP][GAP][tseq[j - 1]])
	    {
	      k = 3;		/* start new alignment - with gap LEGAL!? */
	      fprintf (stderr,
		       "\nErr: This alignment starts in a gap - not even an option!?\n");
	      hit->ali_seq1[l] = '-';
	      hit->ali_ia[l] = ' ';
	      hit->ali_seq2[l++] = index2nt (tseq[--j]);
	      break;
	    }
	  else
	    {
	      printf ("unexpected case in k=2 : %d\n", Iy[i][j]);
	    }
	  hit->ali_seq1[l] = '-';
	  hit->ali_ia[l] = ' ';
	  hit->ali_seq2[l++] = index2nt (tseq[--j]);

	}
      else
	{
	  fprintf (stderr, "\nThis should really NEVER happen!\n");
	}
    }
  hit->ali_seq1[l] = '\0';
  hit->ali_ia[l] = '\0';
  hit->ali_seq2[l] = '\0';

/* reverse sequences in the end*/
/* printf("my length = %d ; altern. test = %d\n", l, strlen(ali_seq1)); */

  reverse (hit->ali_seq1, l - 1);
/*  std::reverse(str, &str[l]); */
  reverse (hit->ali_ia, l - 1);
  reverse (hit->ali_seq2, l - 1);

  hit->qbeg = i + 1;
  hit->qend = maxi;
  hit->tbeg = n + 1 - maxj;
  hit->tend = n - j;
  hit->max = maxval;
/*
  printf("%d - %d\n", i+1, maxi); //alignment in seq1 from to
  printf("%s\n%s\n", ali_seq1, ali_seq2);
  printf("%lu - %lu (3' <-- 5')\n", n-j, n+1-maxj);  //n+1-(j+1)

  printf("Score2fakeE (not considering extpen): %.2f\n", (maxval-559.0)/(-100.0));
  printf("no of nucls in ia: %d + %d = %d\n", maxi-i , maxj-j, maxi-i + maxj-j);
*/
  freeIntMatrix (M, m + 1);
  freeIntMatrix (Ix, m + 1);
  freeIntMatrix (Iy, m + 1);

}

int
RIs_linSpace (unsigned char *qseq,	/* query sequence - numeric representation */
	      unsigned char *tseq,	/* target sequence */
	      int m,		/* query seq length */
	      int n,		/* target seq length */
	      short dsm[6][6][6][6],	/* scoring matrix -- TODO variable length!? */
	      int extensionpenalty,	/* as used in dsm, to calc Score2fakE -- now also a global */
	      int threshold,	/* give out hits higher than that */
	      char *qname,	/* query name */
	      char *tname,	/* target name */
	      char *matname	/* name of the scoring matrix */
  )
{
/* create a hit-struct instead and print from main or other sub?*/
  int **M, **Ix, **Iy;		/* matrices for alignment scores ending in different states */
  int maxi, maxj;
  int maxval, maxk, testmax;	/* maxk not needed!? - max will never be found in gapped anyway! */
  int i, j;

  unsigned char currentRow, lastRow;	/* alternating 0;1 */
  int rowMax_score, *hits_score;
  int rowMax_pos, *hits_pos;
  int tmpQbeg, tmpTend;

  unsigned char *tmpQseq, *tmpTseq;
  short tmpQlen, tmpTlen;
  int nt_count;			/*number of nt in ia to recalc Score2fakE - only tmp no need to store... */
  int tmp, locMax, tmp_min_j;

  IA *maxHit;

  int hitcount = 0;
  double energy;

#if VERBOSE>1
  int tmpi;
#endif

  /* should test if malloc successful! */
  hits_score = malloc (n * sizeof (*hits_score));
  hits_pos = malloc (n * sizeof (*hits_pos));

  tmpQseq = malloc (tblen * sizeof (*tmpQseq));
  tmpTseq = malloc (tblen * sizeof (*tmpTseq));

  maxHit = malloc (sizeof (*maxHit));
  testmax = (int) (1.5 * tblen);
  maxHit->ali_seq1 = malloc (testmax * sizeof (char));
  maxHit->ali_seq2 = malloc (testmax * sizeof (char));
  maxHit->ali_ia = malloc (testmax * sizeof (char));

  M = allocIntMatrix (2, m + 1);	/* (Mis)Match */
  Ix = allocIntMatrix (2, m + 1);	/*Insertion in x(=query), so x paired to gap (in y) */
  Iy = allocIntMatrix (2, m + 1);	/*Insertion(=bulge) in y(=target) */
  maxi = maxj = maxk = 0;
  currentRow = lastRow = 1;

#if VERBOSE>1
  printf ("\t-");
  for (tmpi = 0; tmpi < m; tmpi++)
    printf ("\t%c", index2nt (*(qseq + tmpi)));
#endif

  M[0][0] = Ix[0][0] = Iy[0][0] = 0;

  /*init first row (j=0) -- refers to "-" before first nt of target */

  /* explicitly handling of the boundary condition since i-2 is not defined for i = 1
     No need anymore!
     Ix[0][1] = 0; // MAX(0, dsm[GAP][qseq[0]][GAP][GAP]);      - do not bother with that
     //S(-,*;-,-) should ALWAYS be negative! OR dangling end which can not be incorp. in scoring scheme (not w/o looking further back)
     Iy[0][1] = M[0][1] = NEGINF; // not possible to occure -- TODO: in fact 0 works just as well
   */
  for (i = 1; i <= m; i++)
    {
      Iy[0][i] = M[0][i] = NEGINF;	/* not possible before beginning of target seq */
      Ix[0][i] = 0;		/*MAX(0, dsm[qseq[i-2]][qseq[i-1]][GAP][GAP]); */
      /* do not bother, require to always have a match first!? */
    }

#if VERBOSE>1			/* j = 0 */
  printf ("\n-\t0");
  for (tmpi = 1; tmpi <= m; tmpi++)
    printf ("\t-NI");
/* has been - get back to this if decide not to use NEGINF at all...
    printf("\n-");
    for (tmpi=0; tmpi<=m; tmpi++)
      printf("\t%d", (M[0][tmpi] == NEGINF ? -8 : M[0][tmpi]) ); // -8 as dummy for -inf
  }
 */
#endif

  /* init j=1 row */
  /*init first col (i=0) */
  Iy[1][0] = 0;			/* MAX(0, dsm[GAP][GAP][GAP][tseq[n-1]]); *//*n-1 is last nt in target; used to be 0 after reversion */
  Ix[1][0] = M[1][0] = NEGINF;

  M[1][1] = dsm[GAP][qseq[0]][GAP][tseq[n - 1]];
  /*    maxval = M[1][1] + MAX(0,dsm[qseq[0]][GAP][tseq[n-1]][GAP]); ---- why should we allow a special treatment here!? */
  rowMax_score = M[1][1] + dsm[qseq[0]][GAP][tseq[n - 1]][GAP];
  rowMax_pos = 1;

  /* (1,1) cell can not be in Ix or Iy state. */
  Ix[1][1] = Iy[1][1] = NEGINF;

  for (i = 2; i <= m; i++)
    {

      /* value for M matrix, case we have a pair here (k=0) */
      M[1][i] = dsm[GAP][qseq[i - 1]][GAP][tseq[n - 1]];
/*had been 3 possibilities before, removed case
 Ix[0][i-1] != 0 ? Ix[0][i-1] + dsm[qseq[i-2]][qseq[i-1]][GAP][tseq[n-1]] : -1,
 as above all Ix[0][i] are set to 0
*/

      /* max so far? */
      if ((testmax =
	   M[1][i] + dsm[qseq[i - 1]][GAP][tseq[n - 1]][GAP]) > rowMax_score)
	{
	  rowMax_score = testmax;
	  rowMax_pos = i;
	}

      /* value for Ix matrix, case query sequence paired to gap (k=1) */
      /* removed one option, namely: start new alignment that starts in gap, reflected by (-, Xi; -, -)  */
      Ix[1][i] = max3 (0,
		       /* prev. match, now gap (no match possible before!?)  - add (Xi-1, Xi; Y1, -) */
		       M[1][i - 1] !=
		       0 ? M[1][i - 1] +
		       dsm[qseq[i - 2]][qseq[i - 1]][tseq[n - 1]][GAP] : -1,
		       /* extending existing gap  - add (Xi-1, Xi; -, -) */
		       Ix[1][i - 1] !=
		       0 ? Ix[1][i - 1] +
		       dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] : -1
		       /* OR start new alignment that starts in gap, reflected by (-, Xi; -, -) 
		          dsm[GAP][qseq[i-1]][GAP][GAP] */
	);
/* do NOT allow max other than match state!? -- however (Xi, -; -, -) is positive!?
    testmax = Ix[1][i] + dsm[qseq[i-1]][GAP][GAP][GAP];
    if (testmax > maxval) {
      maxval = testmax;
      maxi = i; maxj = 1; maxk = 1;
    }
*/
      /* value for Iy matrix, case target sequence paired to gap (k=2) */
      /* not possible in this row */
      Iy[1][i] = NEGINF;

    }

  /* initialization of first rows and columns completed
     recursion to complete alignment with two residues follows
   */

#if VERBOSE>1
  printf ("\n%c\t-NI", (index2nt (*(tseq + n - 1))));
  for (tmpi = 1; tmpi <= m; tmpi++)
    printf ("\t%d", M[1][tmpi]);
#endif

  /*no need to store rowMax_score in the first place, can access hits_score[] just as fast!? (keep pointer to current) */
  hits_score[0] = rowMax_score;
  hits_pos[0] = rowMax_pos;
  maxval = rowMax_score;
  maxi = rowMax_pos;
  maxj = 1;			/*can be set in init */

  for (j = 2; j <= n; j++)
    {				/* new row */
      lastRow = currentRow;
      currentRow = j % 2;
      /*changes for pruning:
         switch rows and col / mat[i][j] becomes mat[j][i]
         matIndex  j-1       becomes  lastRow
         matIndex   j        becomes  currentRow
         last nt  tseq[j-2]  becomes  tseq[n-j+1]
         this nt  tseq[j-1]  becomes  tseq[n-j]
       */

      /* TODO better solved with additional running pointer to *(n-j) ?? */

      /* init i=0 column */
      /* no need do do this all over - will not change - no diff if NEGINF or 0 - so already set by init of row j=0 and j=1
         Ix[currentRow][0] = M[currentRow][0] = NEGINF;
         Iy[currentRow][0] = 0; *//*MAX(0, dsm[GAP][GAP][tseq[n-j+1]][tseq[n-j]]); */

      /* init i=1 column */

      /* value for M matrix, case we have a pair here (k=0) */
      M[currentRow][1] = MAX (0, dsm[GAP][qseq[0]][GAP][tseq[n - j]]);
      /* starting a new alignment with (-, X1; -, Yj)  OR  0 if not possible */
      /* coming from gap in tseq NOT possible as we're looking at first position of query! */
      /*had been 3 possibilities before, removed case 'coming from gap in qseq (Iy-matrix), adding (-, X1; Yj-1, Yj)'
         Iy[lastRow][0] != 0 ? Iy[lastRow][0] + dsm[GAP][qseq[0]][tseq[n-j+1]][tseq[n-j]] : -1,
         as above all Iy[j][0] are set to 0
       */

      rowMax_score = M[currentRow][1] + dsm[qseq[0]][GAP][tseq[n - j]][GAP];	/*stricter only if not 0 before!? */
      rowMax_pos = 1;
      /* this would be a single pair w/o stacking, could save this check as well!? */

      /* value for Ix matrix, case query sequence paired to gap (k=1) */
      /* not possible in this column */
      Ix[currentRow][1] = NEGINF;

      /* value for Iy matrix, case target sequence paired to gap (k=2) */
      /* removed one option, namely: start new ali starting w/ gap, (-, -; Yj, -)  // dsm[GAP][GAP][GAP][tseq[n-j]] */
      Iy[currentRow][1] = MAX (
				/*coming from match, opening gap - add (X1, -; Yj-1, Yj) */
				M[lastRow][1] +
				dsm[qseq[0]][GAP][tseq[n - j + 1]][tseq
								   [n - j]],
				/* extending an existing gap - add (-, -; Yj-1, Yj) */
				Iy[lastRow][1] +
				dsm[GAP][GAP][tseq[n - j + 1]][tseq[n - j]]);
/*  do NOT allow max other than match state!? -- however (-, -; *, -) is positive!?
    testmax = Iy[currentRow][1] + dsm[GAP][GAP][tseq[n-j]][GAP];
    if (testmax > maxval) {
      maxval = testmax;
      maxi = 1; maxj = j; maxk = 2;
    }
*/
      /* finished init of i=1 col */

      for (i = 2; i <= m; i++)
	{			/* cols *//* alt bed: *p1  */

	  /* value for M matrix, case we have a pair here (k=0) */
	  M[currentRow][i] = max4 (
				    /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
				    M[lastRow][i - 1] !=
				    0 ? M[lastRow][i - 1] +
				    dsm[qseq[i - 2]][qseq[i - 1]][tseq
								  [n - j +
								   1]][tseq[n
									    -
									    j]]
				    : -1,
				    /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
				    Ix[lastRow][i - 1] +
				    dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq
								       [n -
									j]],
				    /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
				    Iy[lastRow][i - 1] +
				    dsm[GAP][qseq[i - 1]][tseq[n - j + 1]]
				    [tseq[n - j]],
				    /* starting a new alignment with this pair: (-, Xi; -, Yj) */
				    dsm[GAP][qseq[i - 1]][GAP][tseq[n - j]]);
	  if ((testmax =
	       M[currentRow][i] + dsm[qseq[i - 1]][GAP][tseq[n - j]][GAP]) >
	      rowMax_score)
	    {
	      rowMax_score = testmax;
	      rowMax_pos = i;
	    }
	  /* value for Ix matrix, case query paired to gap (k=1) */
	  /* removed one option, namely: start new alignment that starts in gap, reflected by (-, Xi; -, -) */
	  Ix[currentRow][i] = MAX (
				    /*coming from match, add (Xi-1, Xi; Yj, -) */
				    M[currentRow][i - 1] +
				    dsm[qseq[i - 2]][qseq[i - 1]][tseq[n - j]]
				    [GAP],
				    /*extend existing gap, add (Xi-1, Xi; -, -) */
				    Ix[currentRow][i - 1] +
				    dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP]
				    /* start new alignment - starts with gap LEGAL!? 
				       dsm[GAP][qseq[i-1]][GAP][GAP] */
	    );
/* only M[][] can be max, point of backtrack...
      testmax = Ix[currentRow][i] + dsm[qseq[i-1]][GAP][GAP][GAP];
      if (testmax > maxval) {
        maxval = testmax;
        maxi = i; maxj = j; maxk = 1;
      }
*/
	  /* value for Iy matrix, case target paired to gap (k=2) */
	  /* removed one option, namely: start new ali starting w/ gap, (-, -; Yj, -)  // dsm[GAP][GAP][GAP][tseq[n-j]] */
	  Iy[currentRow][i] = MAX (
				    /*coming from match, add (Xi, -; Yj-1, Yj) */
				    M[lastRow][i] +
				    dsm[qseq[i - 1]][GAP][tseq[n - j + 1]]
				    [tseq[n - j]],
				    /*extend existing gap, add (-, -; Yj-1, Yj) */
				    Iy[lastRow][i] +
				    dsm[GAP][GAP][tseq[n - j + 1]][tseq
								   [n - j]]
				    /* start new alignment - starts with gap LEGAL!? 
				       dsm[GAP][GAP][GAP][tseq[n-j]] */
	    );
/*  only M[][] can be max, point of backtrack...
      testmax = Iy[currentRow][i] + dsm[GAP][GAP][tseq[n-j]][GAP];
      if (testmax > maxval) {
        maxval = testmax;
        maxi = i; maxj = j; maxk = 2;
      }
*/
	}

      hits_score[j - 1] = rowMax_score;
      hits_pos[j - 1] = rowMax_pos;

      if (rowMax_score > maxval)
	{
	  maxval = rowMax_score;
	  maxi = rowMax_pos;
	  maxj = j;
	}
#if VERBOSE>1
      printf ("\n%c\t-NI", (index2nt (*(tseq + n - j))));
      for (tmpi = 1; tmpi <= m; tmpi++)
	printf ("\t%d", M[currentRow][tmpi]);
#endif

    }				/*next row j */

#if VERBOSE>1
  printf ("\n");
#endif

/*
  if checking for subopts && set energy threshold, do not guarantee to print best hit first, but only once!
  if checking for subopts and p[2-4], do not check best first.
*/
  if (!(doSubopt && (filterE || printShort > 1)))
    {

#ifdef VERBOSE
      printf ("found maxval %d on pos %d/%d in mat %d\n", maxval, maxi, maxj, maxk);	/* pos are 1-based */
      printf ("equals end pos in query: %d - start pos in target: %d\n", maxi,
	      n + 1 - maxj);
#endif

/*do backtrack for this one only! by recomputing whole matrix for this subsection */
/* max going back tblen */
      tmpQbeg = maxi > tblen - 1 ? maxi - (tblen - 1) : 1;
      tmpTend = maxj > tblen - 1 ? maxj - (tblen - 1) : 1;
      tmpQlen = maxi - tmpQbeg + 1;
      tmpTlen = maxj - tmpTend + 1;

#ifdef VERBOSE
      printf
	("going to realign query %d - %d (len: %hd) vs. target: %d - %d (len: %hd)\n",
	 tmpQbeg, maxi, tmpQlen, n + 1 - tmpTend, n + 1 - maxj, tmpTlen);
#endif

      for (i = 0; i < tmpQlen; i++)
	{
	  tmpQseq[i] = qseq[tmpQbeg - 1 + i];
#ifdef VERBOSE
	  printf ("%c", index2nt (tmpQseq[i]));
#endif

	}
#ifdef VERBOSE
      printf ("\n");
#endif
      for (i = 0; i < tmpTlen; i++)
	{
	  tmpTseq[i] = tseq[n - tmpTend - i];
#ifdef VERBOSE
	  printf ("%c", index2nt (tmpTseq[i]));
#endif
	}
#ifdef VERBOSE
      printf ("\n");
#endif

/*strncpy will not work, as 0 is needed, no terminate etc. -
* no need to reset string as we need to give length anyway... otherwise like follows*/
/*  memset(input_str, '\0', sizeof( input_str )); */

/*With -249 in the Sugimoto case the results we obtain are the same that we get
 from the scripts used for the off-target paper (without weights). 249 corresponds to the value -G-C in the su95 matrix,
 removing it means we do not want to add the energy of adding a firs GC on top of nothing.
 No initialization contribution is subtracted (which instead was the case for the -559 case of the Turner matrices that is
 -150 for the -G-C cell in the matrix, plus -409 as initiation contribution). In the case of the Santa Lucia DNA-DNA
 matrix we applied the same reasoning used for the Sugimoto case, therefore we do
  -363 that is the -G-C case in the matrix. In case another matrix is used options need to be added.*/

      RIs (tmpQseq, tmpTseq, tmpQlen, tmpTlen, dsm, maxHit);

#ifdef VERBOSE
      printf ("%d - %d\n", maxHit->qbeg, maxHit->qend);	/*alignment in subseq1 from to */
      printf ("%s\n%s\n%s\n", maxHit->ali_seq1, maxHit->ali_ia,
	      maxHit->ali_seq2);
      printf ("%d - %d (3' <-- 5')\n", maxHit->tend, maxHit->tbeg);	/* n+1-(j+1) */
#endif
#ifdef VERBOSE
      printf ("full final ali:\n");
#endif
      nt_count =
	maxHit->qend - maxHit->qbeg + 1 + maxHit->tend - maxHit->tbeg + 1;
      //   
      double e1,e2;  
      if(sscanf(matname, "modified %lf %lf", &e1, &e2) == 2){
        energy = (maxHit->max + extensionpenalty * nt_count - e1) / (-e2);          
      }
      else if (!(strcmp (matname, "t99")) || !(strcmp (matname, "t04")))
	{
	  energy =
	    (maxHit->max + extensionpenalty * nt_count - 559.0) / (-100.0);
	}
      else if (!(strcmp (matname, "su95"))
	       || !(strcmp (matname, "su95_noGU")))
	{
	  energy =
	    (maxHit->max + extensionpenalty * nt_count - 249.0) / (-100.0);
	}
      else
	{
	  energy =
	    (maxHit->max + extensionpenalty * nt_count - 363.0) / (-100.0);
	}
    /** TODO :: use maxHit->max OR maxval==hits_score[maxj-1] in output !?!? **/
      if (energy <= maxEnergy)
	{
	  if (printShort == 1)
	    {
	      printf ("%d\t%d\t%d\t%d\t%.2f\n", tmpQbeg + maxHit->qbeg - 1,
		      tmpQbeg + maxHit->qend - 1, n - maxj + maxHit->tbeg,
		      n - maxj + maxHit->tend, energy);
/* to be consistent with other output:
        printf("%d\t%d\t%d\t%d\t%.2f\t%s\n", tmpQbeg+maxHit->qbeg-1, tmpQbeg+maxHit->qend-1, n-maxj+maxHit->tend, n-maxj+maxHit->tbeg, energy, maxHit->ali_ia); 
*/
	    }
	  else if (printShort == 2)
	    {
	      printf ("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", qname,
		      tmpQbeg + maxHit->qbeg - 1, tmpQbeg + maxHit->qend - 1,
		      tname, n - maxj + maxHit->tbeg, n - maxj + maxHit->tend,
		      maxval, energy);
	    }
	  else if (printShort == 3)
	    {
	      hitcount += 1;
	    }
	  else
	    {
	      printf ("Free energy [kcal/mol]: %.2f (%d)\n", energy, maxval);
/*      printf("no of nucls in ia: %lu + %lu = %lu\n", maxHit->qend-maxHit->qbeg+1 , maxHit->tend-maxHit->tbeg+1 , maxHit->qend-maxHit->qbeg+1+maxHit->tend-maxHit->tbeg+1); */

	      printf ("%d - %d\n", tmpQbeg + maxHit->qbeg - 1, tmpQbeg + maxHit->qend - 1);	/*alignment in seq1 from to */
	      printf ("%s\n%s\n%s\n", maxHit->ali_seq1, maxHit->ali_ia,
		      maxHit->ali_seq2);
	      printf ("%d - %d (3' <-- 5')\n", n - maxj + maxHit->tend,
		      n - maxj + maxHit->tbeg);
	    }
	}
    }

/* break here if -s not set */
  if (!doSubopt)
    {
      if (printShort == 3)
	printf ("%s\t%s\t%d\n", qname, tname, hitcount);

      free (maxHit->ali_seq1);
      free (maxHit->ali_seq2);
      free (maxHit->ali_ia);
      free (maxHit);

      free (hits_score);
      free (hits_pos);

      free (tmpQseq);
      free (tmpTseq);

      freeIntMatrix (M, 2);
      freeIntMatrix (Ix, 2);
      freeIntMatrix (Iy, 2);

      return 0;
    }

/* handle suboptimals - BEGIN */

/*MOST naive implementation, puts out one ia for each position in the target,
 given that it is higher than the threshold */

  j = n;
  /* highest array index is j-1, lowest 0 */
  /* as opposed to before were it was 1-n; so have to adjust by 1 here! */
  /* OR change indices before and have array size n+1 */
  /* what about i - still from 1-m as before !? */
  while (j--)
    {

      if (hits_score[j] > threshold)
	{
#ifdef DEBUG
	  printf
	    ("j=%d with %d better than threshold %d - testing neighbors\n", j,
	     hits_score[j], threshold);
#endif
	  tmp = MIN (vicinity, j);
	  tmp_min_j = j - tmp++;
	  locMax = 0;
	  while (--tmp)
	    {
	      if (hits_score[j - tmp] > hits_score[j - locMax])
		{
#ifdef DEBUG
		  printf ("better result %d in distance %d\n",
			  hits_score[j - tmp], tmp);
#endif
		  locMax = tmp;
		}
	    }
	  j -= locMax;

	  /* maxj => j+1 ;; maxi => hits_pos[j] ;; maxval => hits_score[j] */

	  /* do backtrack for this hit, recompute whole matrix in subsection tblen */
	  tmpQbeg = hits_pos[j] > tblen - 1 ? hits_pos[j] - (tblen - 1) : 1;
	  tmpTend = j + 1 > tblen - 1 ? j + 1 - (tblen - 1) : 1;
	  tmpQlen = hits_pos[j] - tmpQbeg + 1;
	  tmpTlen = j + 1 - tmpTend + 1;

#ifdef VERBOSE
	  printf
	    ("\nfound another high hit %d at: end pos in query: %d , start pos in target: %d\n",
	     hits_score[j], hits_pos[j], n - j);
	  printf
	    ("going to realign query %d - %d (len: %hd) vs. target: %d - %d (len: %hd)\n",
	     tmpQbeg, hits_pos[j], tmpQlen, n + 1 - tmpTend, n - j, tmpTlen);
#endif

	  for (i = 0; i < tmpQlen; i++)
	    {
	      tmpQseq[i] = qseq[tmpQbeg - 1 + i];
	    }
	  for (i = 0; i < tmpTlen; i++)
	    {
	      tmpTseq[i] = tseq[n - tmpTend - i];
	    }

	  RIs (tmpQseq, tmpTseq, tmpQlen, tmpTlen, dsm, maxHit);

	  nt_count =
	    maxHit->qend - maxHit->qbeg + 1 + maxHit->tend - maxHit->tbeg + 1;
	  double e1,e2;  
    if(sscanf(matname, "modified %lf %lf", &e1, &e2) == 2){
        energy = (maxHit->max + extensionpenalty * nt_count - e1) / (-e2);          
      }  
	  else if (!(strcmp (matname, "t99")) || !(strcmp (matname, "t04")))
	    {
	      energy =
		(maxHit->max + extensionpenalty * nt_count -
		 559.0) / (-100.0);
	    }
	  else if (!(strcmp (matname, "su95"))
		   || !(strcmp (matname, "su95_noGU")))
	    {
	      energy =
		(maxHit->max + extensionpenalty * nt_count -
		 249.0) / (-100.0);
	    }
	  else
	    {
	      energy =
		(maxHit->max + extensionpenalty * nt_count -
		 363.0) / (-100.0);
	    }

/* TODO anything about this or ignore !?
      if (maxHit->max != hits_score[j]) {
        printf("did some realignment here - probably resulting in a duplicate...\n");
      }
*/
/** TODO :: use maxHit->max OR hits_score[j] in output !?!? **/
	  if (energy <= maxEnergy)
	    {
	      if (printShort == 1)
		{
		  printf ("%d\t%d\t%d\t%d\t%.2f\t%s\n",
			  tmpQbeg + maxHit->qbeg - 1,
			  tmpQbeg + maxHit->qend - 1,
			  n - (j + 1) + maxHit->tend,
			  n - (j + 1) + maxHit->tbeg, energy, maxHit->ali_ia);
		}
	      else if (printShort == 2)
		{
		  printf ("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2f\n", qname,
			  tmpQbeg + maxHit->qbeg - 1,
			  tmpQbeg + maxHit->qend - 1, tname,
			  n - (j + 1) + maxHit->tbeg,
			  n - (j + 1) + maxHit->tend, hits_score[j], energy);
		}
	      else if (printShort == 3)
		{
		  hitcount += 1;
		}
	      else
		{
		  printf ("Free energy [kcal/mol]: %.2f (%d)\n", energy,
			  hits_score[j]);
/*      printf("no of nucls in ia: %lu + %lu = %lu\n", maxHit->qend-maxHit->qbeg+1 , maxHit->tend-maxHit->tbeg+1 , maxHit->qend-maxHit->qbeg+1+maxHit->tend-maxHit->tbeg+1); */
		  printf ("%d - %d\n", tmpQbeg + maxHit->qbeg - 1, tmpQbeg + maxHit->qend - 1);	/*alignment in seq1 from to */
		  printf ("%s\n%s\n%s\n", maxHit->ali_seq1, maxHit->ali_ia,
			  maxHit->ali_seq2);
		  printf ("%d - %d (3' <-- 5')\n", n - (j + 1) + maxHit->tend,
			  n - (j + 1) + maxHit->tbeg);
		}
	    }
	  tmp = j + 1 - maxHit->tend;
	  j = tmp_min_j;	/*next check will start at first position after(before) the range that was tested here */
	  /* alternative: set to end of that backtraced alignment!? */
#ifdef DEBUG
	  printf
	    ("set j=%d, equals pos in seq: %d -- what about starting at end pos of ali: %d (j=%d) -- next attempt at j-1 resp. (pos in seq)++\n",
	     j, n - j, n - tmp, tmp);
#endif
/*      printf ("next attempt at %d not %d", j, j-1-locMax);*/
	}
    }

/* handle suboptimals - END */
  if (printShort == 3)
    printf ("%s\t%s\t%d\n", qname, tname, hitcount);

  free (maxHit->ali_seq1);
  free (maxHit->ali_seq2);
  free (maxHit->ali_ia);
  free (maxHit);

  free (hits_score);
  free (hits_pos);

  free (tmpQseq);
  free (tmpTseq);

  freeIntMatrix (M, 2);
  freeIntMatrix (Ix, 2);
  freeIntMatrix (Iy, 2);

  return 0;
}

void
RIs_force_start_end_init (unsigned char *qseqIx,	/* query sequence - numeric representation */
			  unsigned char *tseqIx,	/* target sequence  */
			  int len_seq1,	/* query seq length */
			  int len_seq2,	/* target seq length */
			  short dsm[6][6][6][6],	/* scoring matrix */
			  char *matname	/* name of the scoring matrix */
  )
{
  int i;
  unsigned char *tmp;
  int testmax;
  float *noweight;

  IA *maxHit;
  tmp = malloc ((len_seq2) * sizeof (tseqIx));
  maxHit = malloc (sizeof (IA));
  testmax = (int) (1.5 * len_seq2);
  maxHit->ali_seq1 = malloc (testmax * sizeof (char));
  maxHit->ali_seq2 = malloc (testmax * sizeof (char));
  maxHit->ali_ia = malloc (testmax * sizeof (char));

  /*reverting seq2 (target) */
  for (i = 0; i < len_seq2; i++)
    {
      tmp[i] = tseqIx[len_seq2 - 1 - i];
    }
  if (strcmp (pos_weights, "noweights") && force_start_val == 0)
    {
      fprintf (stderr,
	       "ERR: weighted interactions must be forced to start at position 0 in the target. Please increase -f parameter, or use array \"noweights\" for option -w.\n");
      exit (1);
    }
  if (!strcmp (pos_weights, "CRISPR_20nt_5p_3p"))
    {
      extern float wC20_5p_3p[19];
      extern int size_wC20_5p_3p;
      if (size_wC20_5p_3p < len_seq1 - 1)
	{
	  fprintf (stderr,
		   "ERR: the array of weights is too short for the given query.\n Weights length mush be >= than the length of the query -1.\n");
	  exit (1);
	}
      if (size_wC20_5p_3p == len_seq1 - 1)
	{
	  RIs_force_start_end_weighted (qseqIx, tmp, len_seq1, len_seq2,
					&wC20_5p_3p[0], dsm, maxHit, matname);
	}
      else
	{
	  RIs_force_start_end_weighted (qseqIx, tmp, len_seq1, len_seq2,
					&wC20_5p_3p[size_wC20_5p_3p + 1 -
						    len_seq1], dsm, maxHit,
					matname);
	}
    }				/*
				   else if (!strcmp(pos_weights, "test")){
				   extern float test[4];
				   extern int size_test;
				   if (size_test<len_seq1-1){
				   fprintf(stderr,"ERR: the array of weights is too short for the given query.\n Weights length mush be >= than the length of the query -1.\n");
				   exit(1);
				   }
				   if (size_test==len_seq1-1){
				   RIs_force_start_end_weighted(qseqIx,tmp,len_seq1,len_seq2,&test[0],size_test,dsm,maxHit, matname);
				   }
				   else{
				   RIs_force_start_end_weighted(qseqIx,tmp,len_seq1,len_seq2,&test[size_test+1-len_seq1],len_seq1-1,dsm,maxHit, matname);
				   }
				   } */
  else if (!strcmp (pos_weights, "noweights"))
    {
      noweight = malloc ((len_seq1 - 1) * sizeof (float));
      for (i = 0; i < (len_seq1 - 1); ++i)
	{
	  noweight[i] = 1.0;
	}
      RIs_force_start_end_weighted (qseqIx, tmp, len_seq1, len_seq2,
				    &noweight[0], dsm, maxHit, matname);
    }
  else
    {
      fprintf (stderr,
	       "Undefined weights array. Existing weights verctors are CRISPR_20nt_5p_3p and noweights. To add a new weights vector, create an array in weights.c and an ad-hoc if clause. Use noweights to set all weights to 1.\n");
      exit (1);
    }
}

void
RIs_force_start_end_weighted (unsigned char *qseq,	/* query sequence - numeric representation */
			      unsigned char *tseq,	/* target sequence - reversed */
			      int m,	/* query seq length */
			      int n,	/* target seq length */
			      float *weights,	/*array of weights */
			      short dsm[6][6][6][6],	/* scoring matrix */
			      IA * hit,	/* pointer to struct, fill results */
			      char *matname)
{

  float **M, **Ix, **Iy;	/* matrices for alignment scores ending in different states */
  int i, j, k, colj;
  float w1, w2;
  int startpos_found;		/* used in backtrack, true if start pos has been reached */
  float mVal, xVal, yVal;	/* values coming from M, Ix, Iy. Not possible to start a new alignment at any position, due to the force start. */
  char *query_alignment, *target_alignment;	/*used to store query and target alignment symbols in backtracking */
  float max_score;		/*max score found in a certain column */
  float energy;			/*energy corresponding to a given max_score */

  M = allocFloatMatrix (m + 1, n + 1);	/* (Mis)Match */
  Ix = allocFloatMatrix (m + 1, n + 1);	/* Insertion in x(=query), so x paired to gap (in y) */
  Iy = allocFloatMatrix (m + 1, n + 1);	/* Insertion(=bulge) in y(=target) */

  Ix[0][0] = Iy[0][0] = 0;
  M[0][0] = force_start_val;
  /*init first row (j=0) -- this is COL - change throughout!? */

  /*filling column zero */
  for (i = 1; i <= m; i++)
    {
      Iy[i][0] = NEGINF;	/* not possible before beginning of target seq */
      Ix[i][0] = 0;		/*MAX(0, dsm[qseq[i-2]][qseq[i-1]][GAP][GAP]); *//* require to always have a match first! */
      M[i][0] = force_start_val;	/*To force interaction to start at column 1 (-N in the DNA) */
    }

  /*filling row zero */
  for (j = 1; j <= n; j++)
    {
      Ix[0][j] = M[0][j] = NEGINF;
      Iy[0][j] = 0;		/*MAX(0, dsm[GAP][GAP][tseq[j-2]][tseq[j-1]]); */
    }

  /* The initialization of i=1 column and j=1 row have to be handled explicitly
     since at this point we do not have two residues to use.
     Handle the (1,1) cell explicitly since the boundary recursion includes (i-2) or (j-2) cases. */
  M[1][1] = force_start_val + (float) dsm[GAP][qseq[0]][GAP][tseq[0]];	/*no weight used, first base pair */
  Ix[1][1] = Iy[1][1] = NEGINF;

  /* Filling first column */
  for (i = 2; i <= m; i++)
    {
      /* value for M matrix, case we have a pair here (k=0) */
      /* first letter of target sequence, so it can ONLY come from gapped (upstream nucleotides in query have not been used) */
      M[i][1] = force_start_val + (float) dsm[GAP][qseq[i - 1]][GAP][tseq[0]];	/* MAX(0, ); no weight used, first base pair */
      /* value for Ix matrix, case query sequence paired to gap (k=1) */
      mVal =			/* prev. match, now gap - add (Xi-1, Xi; Y1, -) */
	M[i - 1][1] !=		/*adding a gap on top of a bp, weights used */
	force_start_val ? M[i - 1][1] +
	dsm[qseq[i - 2]][qseq[i - 1]][tseq[0]][GAP] * weights[i - 2] : -1;
      xVal =			/* extending existing gap  - add (Xi-1, Xi; -, -) */
	Ix[i - 1][1] != 0 ?	/*adding a gap on top of a gap, weights used */
	Ix[i - 1][1] + dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] * weights[i - 2] : -1;
      /* OR start new alignment that starts in gap, reflected by (-, Xi; -, -). Forbidden! */
      /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */
      /*Ix[i][1] = max4(mVal,xVal,nVal,0); */
      Ix[i][1] = max3f (mVal, xVal, 0);	/*removed nVal option */
      /* value for Iy matrix, case target sequence paired to gap (k=2) not possible in this row */
      Iy[i][1] = NEGINF;
    }

  /* Filling first row */
  for (j = 2; j <= n; j++)
    {
      /* value for M matrix, case we have a pair here (k=0 */
      /* coming from gap in qseq (Iy-matrix); First letter in query sequence, adding (-, X1; Yj-1, Yj) */
      M[1][j] = dsm[GAP][qseq[0]][GAP][tseq[j - 1]];	/* MAX(0, ); no weight used, first base pair */
      /* value for Iy matrix, case target sequence paired to gap (k=2) */
      /*coming from match, opening gap - add (X1, -; Yj-1, Yj) */
      mVal = M[1][j - 1] !=	/*adding a gap on top of a bp, weights used */
	0 ? M[1][j - 1] +
	(float) dsm[qseq[0]][GAP][tseq[j - 2]][tseq[j - 1]] * weights[0] : -1;
      /* extending an existing gap - add (-, -; Yj-1, Yj) */
      yVal = Iy[1][j - 1] != 0 ?	/*adding a gap on top of a gap, weights used */
	Iy[1][j - 1] + dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] *
	weights[0] : -1;
      /* or start new ali, begin with gap -- forbidden ! */
      /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */
      Iy[1][j] = max3f (mVal, yVal, 0);	/*removed nVal option */
      /* value for Ix matrix, case query sequence paired to gap (k=1) */
      /* not possible in this column */
      Ix[1][j] = NEGINF;
    }
  /* initialization of first rows and columns completed
     recursion to complete alignment with two residues follows */
  for (i = 2; i <= m; i++)
    {
      for (j = 2; j <= n; j++)
	{
	  /* value for M matrix, case we have a pair here (k=0) */
	  /* coming from a match, add (Xi-1, Xi; Yi-1, Yi) */
	  mVal =
	    M[i - 1][j - 1] !=
	    0 ? M[i - 1][j - 1] +
	    ((float) dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] *
	     weights[i - 2]) : -1;
	  /* coming from gap in target, add (Xi-1, Xi; -, Yi) */
	  xVal =
	    Ix[i - 1][j - 1] !=
	    0 ? Ix[i - 1][j - 1] +
	    ((float) dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq[j - 1]] *
	     weights[i - 2]) : -1;

	  /* coming from gap in query, add (-, Xi; Yi-1, Yi) */
	  yVal =
	    Iy[i - 1][j - 1] !=
	    0 ? Iy[i - 1][j - 1] +
	    ((float) dsm[GAP][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] *
	     weights[i - 2]) : -1;

	  /* starting a new alignment with this pair: (-, Xi; -, Yj) */
	  M[i][j] = max4f (mVal, xVal, yVal, 0);
	  /* value for Ix matrix, case query paired to gap (k=1) */
	  /*coming from match, add (Xi-1, Xi; Yj, -) */
	  mVal =		/*add a bulge on top of a bp */
	    M[i - 1][j] !=
	    0 ? M[i - 1][j] +
	    (float) dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 1]][GAP] *
	    weights[i - 2] : -1;
	  /*extend existing gap, add (Xi-1, Xi; -, -) */
	  xVal =
	    Ix[i - 1][j] !=
	    0 ? Ix[i - 1][j] +
	    (float) dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] * weights[i -
								      2] : -1;
	  /* start new alignment - starts with gap -> ILLEGAL! */
	  /* nVal = dsm[GAP][qseq[i-1]][GAP][GAP]; */

	  Ix[i][j] = max3f (mVal, xVal, 0);	/*removed option nVal */
	  /* value for Iy matrix, case target paired to gap (k=2) */
	  /*coming from match, add (Xi, -; Yj-1, Yj) */
	  w2 = i - 2 < m - 1 - 1 ? weights[i - 1] : weights[m - 1 - 1];	/*number of weights=length query -1) */
	  mVal = M[i][j - 1] != 0 ? M[i][j - 1] + (float) dsm[qseq[i - 1]][GAP][tseq[j - 2]][tseq[j - 1]] * ((weights[i - 2] + w2) / 2) : -1;	/*use average of weights for bulges in target */
	  /*extend existing gap, add (-, -; Yj-1, Yj) */
	  yVal = Iy[i][j - 1] != 0 ? Iy[i][j - 1] + (float) dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] * ((weights[i - 2] + w2) / 2) : -1;	/* No weight used for continuing bulges */
	  /* start new alignment - starts with gap LEGAL!? */
	  /* nVal = dsm[GAP][GAP][GAP][tseq[j-1]]; */

	  Iy[i][j] = max3f (mVal, yVal, 0);	/*removed option nVal */
#ifdef DEBUG
	  printf ("M %f; %i; %i\n", weights[i - 2], i, j);
	  printf ("Y %f; %i; %i\n", weights[i - 2], i, j);
	  printf ("X %f; %i; %i\n", weights[i - 2], i, j);
#endif
	}
    }
#ifdef DEBUG
  printf ("M matrix:\n");
  printfloatMat (M, m + 1, n + 1, qseq, tseq);
  printf ("Bq matrix:\n");
  printfloatMat (Ix, m + 1, n + 1, qseq, tseq);
  printf ("Bt matrix:\n");
  printfloatMat (Iy, m + 1, n + 1, qseq, tseq);
#endif
  printf ("***Structures and Energies***\n");
  /*backtrack */
  for (colj = n; colj <= n; colj++)
    {				/*n is the number of columns, m is the number of rows.  NOTE: here we put colj = n because we only want the interaction that ends at position n (ends at the 5' of the target sequence, which is in 3'-5' direction) */
      /*printf("%d\n",m);
         printf("%d\n",n); */
      query_alignment = allocate_char_array (m);	/* no need to allocate +1 for \0 since the matrix already contain 1 position more for - */
      target_alignment = allocate_char_array (n);	/* no need to allocate +1 for \0 since the matrix already contain 1 position more for - */
      /*printf("col : %d\n",colj); */
      max_score = find_max_value_f (M, Ix, Iy, &k, &i, colj, m, dsm, qseq, tseq);	/*given arrays representing columns, we want to find the row position in which the max score is reached. K is 0-1-2 M Ix Iy - should always be 0 to begin with. NOTE: in this new version we search only the maximums in row m (last row) because we want to force the interaction to end at 3' of query */
      /*note that i and k are modified in place in the function, since we are passing the pointers */
      j = colj;
      startpos_found = 0;
#ifdef DEBUG
      printf ("max score %f\n", max_score);
#endif
      while (i > 0 && j > 0 && !startpos_found)
	{
	  /*printf("%s\n",query_alignment);
	     printf("%s\n",target_alignment); */
	  /* find which cell gave score */
	  if (k == 0)
	    {			/* highest score in M matrix, having a match! */
	      if (M[i][j] ==
		  (dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]]) + force_start_val
		  && j == 1)
		{		/*      started new alignment */
		  /*printf("means we started alignment at row [%d] col [%d] of matrix %d \n", i, j, k); */
		  startpos_found = 1;
		}
	      else if (M[i][j] == (dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]])
		       && j != 1)
		{		/*      started new alignment */
		  /*printf("means we started alignment at row [%d] col [%d] of matrix %d \n", i, j, k); */
		  fprintf (stderr,
			   "Force start option did not work: try to increase the number given to -f.\nTry with a number > %d\n",
			   max3 (2000, m * 500, n * 500));
		  exit (1);
		}
	      else if (M[i][j] ==
		       M[i - 1][j - 1] +
		       ((float)
			dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 2]][tseq
								   [j -
								    1]] *
			weights[i - 2]) && i != 1 && j != 1)
		{
		  /*printf("previous was also match\n"); *//*      previous was also match */
		}
	      else if (M[i][j] ==
		       Ix[i - 1][j - 1] +
		       dsm[qseq[i - 2]][qseq[i - 1]][GAP][tseq[j - 1]] *
		       weights[i - 2] && i != 1)
		{
		  /*printf("coming from gap in seq2/target\n"); *//*      coming from gap in seq2/target */
		  k = 1;
		}
	      else if (M[i][j] ==
		       Iy[i - 1][j - 1] +
		       dsm[GAP][qseq[i - 1]][tseq[j - 2]][tseq[j - 1]] *
		       weights[i - 2] && j != 1)
		{
		  /*printf("coming from gap in seq1/query\n"); *//*      coming from gap in seq1/query */
		  k = 2;
		}
	      else if (M[i][j] ==
		       (dsm[GAP][qseq[i - 1]][GAP][tseq[j - 1]]) +
		       force_start_val && j != 1)
		{		/*      started new alignment */
		  fprintf (stderr,
			   "Force start option did not work: try to increase the number given to -f.\nTry with a number > %d\n",
			   max3 (2000, m * 500, n * 500));
		  exit (1);
		}
	      else
		{
		  fprintf (stderr, "ERR: unexpected value in k=0.\n");
		  exit (1);
		}
	      set_alignment_symbols (index2nt (qseq[i - 1]),
				     index2nt (tseq[j - 1]),
				     &query_alignment[i - 1],
				     &target_alignment[j - 1]);
	      i = i - 1;
	      j = j - 1;
	    }
	  else if (k == 1)
	    {			/* seq1(query) paired to a gap (in target). Bulge on top of bp */
	      if (Ix[i][j] ==
		  M[i - 1][j] +
		  (float) dsm[qseq[i - 2]][qseq[i - 1]][tseq[j - 1]][GAP] *
		  weights[i - 2])
		{
		  k = 0;	/* open a new gap coming from match */
		  /*printf("Gap in Ix coming from match %f\n", Ix[i][j]); */
		}
	      else
		{
		  if (Ix[i][j] ==
		      Ix[i - 1][j] +
		      (float) dsm[qseq[i - 2]][qseq[i - 1]][GAP][GAP] *
		      weights[i - 2])
		    {		/*gap on top of gap */
		      if (i == 1)
			{
			  k = 3;
			  startpos_found = 1;
			}
		      else
			{
			  k = 1;	/* extend existing gap */
			}
		      /*printf("Gap in Ix coming from gap %f\n", Ix[i][j]); */
		    }
		  else
		    {
		      fprintf (stderr,
			       "ERR: unexpected case in k=1 : %f in row %d and col %d\n",
			       Ix[i][j], i, j);
		      exit (1);
		    }
		}
	      query_alignment[i - 1] = 'B';
	      i = i - 1;
	    }
	  else if (k == 2)
	    {			/* seq2(target) paired to a gap (in query). bulge on top of bp. */
	      if (i == 1)
		{
		  w1 = weights[0];
		  w2 = weights[0];
		}
	      else
		{
		  w1 = weights[i - 2];
		  w2 = i - 2 < m - 1 - 1 ? weights[i - 1] : weights[m - 1 - 1];	/*number of weights=length query -1) */
		}
	      if (Iy[i][j] ==
		  M[i][j - 1] +
		  (float) dsm[qseq[i - 1]][GAP][tseq[j - 2]][tseq[j - 1]] *
		  ((w1 + w2) / 2))
		{
		  k = 0;	/* open a new gap coming from match */
		  /*printf("Gap in Iy coming from match\n"); */
		}
	      else
		{
		  if (Iy[i][j] ==
		      Iy[i][j - 1] +
		      (float) dsm[GAP][GAP][tseq[j - 2]][tseq[j - 1]] *
		      ((w1 + w2) / 2))
		    {
		      if (j == 1)
			{
			  k = 3;
			  startpos_found = 1;
			}
		      else
			{
			  k = 2;	/* extend existing gap */
			}
		    }
		  else
		    {
		      fprintf (stderr, "unexpected case in k=2 : %f\n",
			       Iy[i][j]);
		      exit (1);
		    }
		}
	      target_alignment[j - 1] = 'B';
	      j = j - 1;
	    }
	  else
	    {
	      fprintf (stderr, "\nThis should really NEVER happen!\n");
	      exit (1);
	    }
	}
      if (!startpos_found)
	{
	  if (j > 0)
	    {
	      fprintf (stderr,
		       "Force start option did not work: try to increase the number given to -f.\nTry with a number > %d.\n",
		       (int) max3 (2000, m * 200, n * 200));
	      exit (1);
	    }
	}
      printf ("%s\t%s\n", query_alignment, target_alignment);
	  double e1,e2;  
    	if(sscanf(matname, "modified %lf %lf", &e1, &e2) == 2){
        energy = (max_score - force_start_val - e1) / (e2);         
      }   
   else if (!(strcmp (matname, "t99")) || !(strcmp (matname, "t04")))
	{
	  energy = (max_score - force_start_val - 559.0) / (-100.0);
	}
      else if (!(strcmp (matname, "su95"))
	       || !(strcmp (matname, "su95_noGU")))
	{
	  energy = (max_score - force_start_val - 249.0) / (-100.0);
	}
      else
	{
	  energy = (max_score - force_start_val - 363.0) / (-100.0);
	}
      printf ("Free energy [kcal/mol] (No extension penalty): %.2f (%f)\n",
	      energy, max_score - force_start_val);
    }

  freeFloatMatrix (M, m + 1);
  freeFloatMatrix (Ix, m + 1);
  freeFloatMatrix (Iy, m + 1);
}
