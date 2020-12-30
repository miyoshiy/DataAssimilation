#include <stdio.h>
#include <ctype.h>
#include <string.h>
#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */
#include <stdlib.h>
#endif /* ANSI */

#define MAXLINE 256
#define MAXFILE 256

#define SEPCHAR '/'
#define SEPSTR "/"
#define SUFFIX ".f"
#define NAME "nrtodp.dat"
#ifndef LIBRARY_DIR
#define LIBRARY_DIR "/usr/lib"
#endif
static char usage[]=
"Usage: nrtodp path file [file...]\n\
  Performs a pro-forma conversion of Numerical Recipes FORTRAN files\n\
  to DOUBLE PRECISION.  Input file names are the 2nd and subsequent\n\
  arguments.  Output files are obtained by prepending \"path\" and the\n\
  letter \"d\".  Example: the command\n\
     nrtodp /nr/double/ lubksb.f\n\
  creates an output file \"/nr/double/dlubksb.f\n";

#define TAB '\011'

#ifdef NEED_STRSTR
char *strstr();
#endif /* NEED_STRSTR */

main(argc,argv)
int argc;
char *argv[];
{
	char lin[MAXLINE], lout[MAXLINE], outnam[MAXFILE], barnam[MAXFILE];
	char *keytext, *next, *endtext, ch, *loutc, *preal, *pname, *xstring;
	char *xname, oldtext[MAXLINE], newtext[MAXLINE];
	char *lini, *louti;
	int il,ip,ir,ie,inc,j,k,iarg,exquote;
	long int xlength;
	FILE *infile, *outfile, *xfile;
	
	if (argc < 3) { fputs(usage,stderr); exit(1); }
	/* open and digest the data file */

	if ((xfile = fopen(NAME,"r")) == NULL) {
		/* First try current directory, then library directory */
		if ((xname = (char *) malloc((unsigned) (strlen(LIBRARY_DIR) + strlen(NAME) + 2))) == NULL) {
			fprintf(stderr,"Can't get enough memory for date file name.\n");
			exit(1);
		}
		strcpy(xname,LIBRARY_DIR);
		strcat(xname,SEPSTR);
		strcat(xname,NAME);
		if ((xfile = fopen(xname,"r")) == NULL) {
			fprintf(stderr,"Can't open data file %s\n",xname);
			exit(1);
		}
	}
	fseek(xfile,0L,2);
	xlength = ftell(xfile);
	fseek(xfile,0L,0);
	if ((xstring = (char *) malloc((unsigned) (xlength+2L))) == NULL) {
		fprintf(stderr,"Can't get enough memory for data file.\n");
		exit(1);
	}
	if (fread(xstring,1,xlength,xfile) == 0) {
		fprintf(stderr,"Error on reading data file.\n");
		exit(1);
	}
	xstring[xlength] = '\0';
	close(xfile);
	
	/* now process the individual files*/
	for (iarg=2;iarg<argc;iarg++) {
		if ((infile = fopen(argv[iarg],"r")) == NULL) continue;
		/* parse out the bare routine name */
		strcpy(barnam,argv[iarg]);
		if ((preal = strrchr(barnam,'.')) != NULL) *preal = '\0';
		if ((preal = strrchr(barnam,SEPCHAR)) != NULL) pname = preal + 1;
		else pname = barnam;
		strcpy(outnam,argv[1]);
		strcat(outnam,"d");
		strcat(outnam,pname);
		strcat(outnam,SUFFIX);
		if ((outfile = fopen(outnam,"w")) == NULL) {
			fprintf(stderr,"Can't open output file %s\n",outnam);
			exit(1);
		}
		fprintf(stderr,"Pro-forma converting %s to %s\n",argv[iarg],outnam);
		/* get the first exception */
		strcpy(newtext,"&");
		strcat(newtext,pname);
		keytext=strstr(xstring,newtext);
		if (keytext) {
			endtext=strchr(keytext+1,'&');
			keytext=strchr(keytext,':')+1;  /* 1st char of old */
			next=strchr(keytext,':');       /* separator */
			*next='\0';
			strcpy(oldtext,keytext);
			*next=':';
			keytext=next+1;                  /* 1st char of new */
			next=strchr(keytext,':');       /* terminator */
			*next='\0';
			strcpy(newtext,keytext);
			*next=':';
			keytext=strchr(next+1,':')+1;      /* 1st char of next old */
		}
		for (;;) {
			if (fgets(lin,MAXLINE,infile) == NULL) break;
			while (keytext && (next=strstr(lin,oldtext)) != NULL) {
				ch = *next;
				*next = '\0';
				strcpy(lout,lin);
				strcat(lout,newtext);
				strcat(lout,next+strlen(oldtext));
				strcpy(lin,lout);
				if (keytext > endtext) keytext = NULL;
				else {
					next=strchr(keytext,':');       /* separator */
					*next='\0';
					strcpy(oldtext,keytext);
					*next=':';
					keytext=next+1;                  /* 1st char of new */
					next=strchr(keytext,':');       /* terminator */
					*next='\0';
					strcpy(newtext,keytext);
					*next=':';
					keytext=strchr(next+1,':')+1;      /* 1st char of next old */
				}
			}
			if ((preal = strstr(lin,"REAL ")) != NULL) {
				ie = 1;
				for (loutc = preal-1; loutc >= lin; loutc--) {
					if (*loutc != ' ' && *loutc != TAB) ie=0;
				}
				if (ie) {
					*preal = '\0';
					strcpy(lout,lin);
					strcat(lout,"DOUBLE PRECISION ");
					strcat(lout,preal+5);
					strcpy(lin,lout);
				}
			}
			if ((preal = strstr(lin,"COMPLEX ")) != NULL) {
				ie = 1;
				for (loutc = preal-1; loutc >= lin; loutc--) {
					if (*loutc != ' ' && *loutc != TAB) ie=0;
				}
				if (ie) {
					*preal = '\0';
					strcpy(lout,lin);
					strcat(lout,"COMPLEX*16 ");
					strcat(lout,preal+8);
					strcpy(lin,lout);
				}
			}
#define REPLACE_SIMPLE(old,new,length_old) \
			while ((preal = strstr(lin,old)) != NULL) \
				strncpy(preal,new,length_old)

			REPLACE_SIMPLE("sngl(","dble(",5);
			REPLACE_SIMPLE("float("," dble(",6);
			REPLACE_SIMPLE("aimag(","dimag(",6);

/* This is for the case where the original is a subset of the replacement */
/* Note this will not catch a mixed-precision line of the form dx() + x() */
#define REPLACE_SUB(old,new,length_old) \
			lini = lin; louti = lout; \
			while ((preal = strstr(lini,old)) != NULL \
			       && *(preal-1) != 'd' ) { \
				*preal = '\0'; \
				inc = preal-lini+length_old+1; \
				strcpy(louti,lini); \
				strcat(louti,new); \
				strcat(louti,preal+length_old); \
				strcpy(lini,louti); \
				lini += inc; \
				louti += inc; \
			}
			REPLACE_SUB("real(","dreal(",5);
			REPLACE_SUB("cmplx(","dcmplx(",6);

			/* this is a weird special case in ran4: */
			if ((preal = strstr(lin,"DOUBLE PRECISION ftemp")) != NULL)
			strcpy(preal,"REAL ftemp\n");
			/* end of special case */
			loutc = lout;
			exquote = 1;
			for (ip=0; lin[ip] != '\0'; ip++) {
				*loutc++ = lin[ip];
				if (lin[ip] == '\'') exquote = 1 - exquote;
				if (exquote && lin[ip] == '.') {
					il = ir = ip;
					while (il>0 && isdigit(lin[il-1])) --il;
					while (isdigit(lin[ir+1])) ++ir;
					if (il == ir) continue;
					if (il > 0 && isalpha(lin[il-1])) continue;
					ie = 0;
					if (toupper(lin[ir+1]) == 'E') {
						ie = ir+1;
						j = ie;
						if (lin[j+1] == '+' || lin[j+1] == '-') j++;
						if ( ! isdigit(lin[j+1]) ) ie = 0;
					}
					if (ie == 0 && isalpha(lin[ir+1])) continue;
					if (ie) {
						lin[ie] = 'd';
					} else {
						while (ip < ir) *loutc++ = lin[++ip];
						*loutc++ = 'd';
						*loutc++ = '0';
					}
				}
			}
			*loutc = '\0';
			k=strlen(lout);
			if (k <= 73  || lout[0] == 'C') {
				/* 73 because of the newline at the end */
				fprintf(outfile,"%s",lout);
			} else {
				j = 72;
				while ((isalnum(lout[j]) || lout[j]=='.') && j>0) j--;
				ch = lout[j];
				lout[j] = '\0';
				fprintf(outfile,"%s\n     *%c%s",lout,ch,&lout[j+1]);
			}
		}
		fclose(infile);
		fclose(outfile);
	}
	exit(0);
}
