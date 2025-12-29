static char adSid[]="$Id: adStack.c 3723 2011-02-24 13:34:42Z llh $";
/*
Replac the faulty PUSH/POPs by assignments to/from local arrays,
each of them being declared localized to its thread, make it
thread-safe for OPENMP.

If the OPENMP THREAD NUMBER IS GREATER THAN 256 on a single node,
the STACK_POOL_SEZE needs to be increased.

Please beware of this portion that just amended the necessary codes
which used by WRFPLUS to make it thread-safe, NOT THE WHOLE CODES.
Check it carefully before using for your own case.

2011-10-10 jliu@ucar.edu

*/

#ifndef CRAY
# ifdef NOUNDERSCORE
#      define PUSHINTEGER4ARRAY pushinteger4array
#      define POPINTEGER4ARRAY popinteger4array
#      define PUSHREAL8ARRAY pushreal8array
#      define POPREAL8ARRAY popreal8array
# else
#      define PUSHINTEGER4ARRAY pushinteger4array_
#      define POPINTEGER4ARRAY popinteger4array_
#      define PUSHREAL8ARRAY pushreal8array_
#      define POPREAL8ARRAY popreal8array_
# endif
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
#ifdef _OPENMP
*/
#include <omp.h>
#define STACK_POOL_SIZE 256
/*
#else
#define omp_get_thread_num() 0
#define STACK_POOL_SIZE 1
#endif
*/

#define ONE_BLOCK_SIZE 16384
#ifndef STACK_SIZE_TRACING
#define STACK_SIZE_TRACING 1
#endif
/* The main stack is a double-chain of DoubleChainedBlock objects.
 * Each DoubleChainedBlock holds an array[ONE_BLOCK_SIZE] of char. */
typedef struct _doubleChainedBlock{
  struct _doubleChainedBlock *prev ;
  char                       *contents ;
  struct _doubleChainedBlock *next ;
} DoubleChainedBlock ;

/* Globals that define the current position in the stack: */
static DoubleChainedBlock *curStack[STACK_POOL_SIZE] = {NULL} ;
static char               *curStackTop[STACK_POOL_SIZE]    = {NULL} ;
/* Globals that define the current LOOKing position in the stack: */
static DoubleChainedBlock *lookStack[STACK_POOL_SIZE] = {NULL} ;
static char               *lookStackTop[STACK_POOL_SIZE] = {NULL} ;

static long int mmctraffic = 0 ;
static long int mmctrafficM = 0 ;
#ifdef STACK_SIZE_TRACING
long int bigStackSize = 0;
#endif

/* PUSHes "nbChars" consecutive chars from a location starting at address "x".
 * Resets the LOOKing position if it was active.
 * Checks that there is enough space left to hold "nbChars" chars.
 * Otherwise, allocates the necessary space. */
void pushNarray(char *x, unsigned int nbChars) {
  int thread_id = omp_get_thread_num();
  unsigned int nbmax = (curStack[thread_id])?ONE_BLOCK_SIZE-(curStackTop[thread_id]-(curStack[thread_id]->contents)):0 ;
#ifdef STACK_SIZE_TRACING
  bigStackSize += nbChars;
#endif

  mmctraffic += nbChars ;
  while (mmctraffic >= 1000000) {
     mmctraffic -= 1000000 ;
     mmctrafficM++ ;
  }

  lookStack[thread_id] = NULL ;
  if (nbChars <= nbmax) {
    memcpy(curStackTop[thread_id],x,nbChars) ;
    curStackTop[thread_id]+=nbChars ;
  } else {
    char *inx = x+(nbChars-nbmax) ;
    if (nbmax>0) memcpy(curStackTop[thread_id],inx,nbmax) ;
    while (inx>x) {
      if ((curStack[thread_id] == NULL) || (curStack[thread_id]->next == NULL)) {
        /* Create new block: */
	DoubleChainedBlock *newStack ;
	char *contents = (char*)malloc(ONE_BLOCK_SIZE*sizeof(char)) ;
	newStack = (DoubleChainedBlock*)malloc(sizeof(DoubleChainedBlock)) ;
	if ((contents == NULL) || (newStack == NULL)) {
	  DoubleChainedBlock *stack = curStack[thread_id] ;
	  int nbBlocks = (stack?-1:0) ;
	  while(stack) {
	      stack = stack->prev ;
	      nbBlocks++ ;
	  }
	  printf("Out of memory (allocated %i blocks of %i bytes)\n",
		 nbBlocks, ONE_BLOCK_SIZE) ;
          exit(0);
	}
	if (curStack[thread_id] != NULL) curStack[thread_id]->next = newStack ;
	newStack->prev = curStack[thread_id] ;
	newStack->next = NULL ;
	newStack->contents = contents ;
	curStack[thread_id] = newStack ;
        /* new block created! */
      } else
	curStack[thread_id] = curStack[thread_id]->next ;
      inx -= ONE_BLOCK_SIZE ;
      if(inx>x)
	memcpy(curStack[thread_id]->contents,inx,ONE_BLOCK_SIZE) ;
      else {
	unsigned int nbhead = (inx-x)+ONE_BLOCK_SIZE ;
	curStackTop[thread_id] = curStack[thread_id]->contents ;
	memcpy(curStackTop[thread_id],x,nbhead) ;
	curStackTop[thread_id] += nbhead ;
      }
    }
  }
}

/* POPs "nbChars" consecutive chars to a location starting at address "x".
 * Resets the LOOKing position if it was active.
 * Checks that there is enough data to fill "nbChars" chars.
 * Otherwise, pops as many blocks as necessary. */
void popNarray(char *x, unsigned int nbChars) {
  int thread_id = omp_get_thread_num();
  unsigned int nbmax = curStackTop[thread_id]-(curStack[thread_id]->contents) ;
#ifdef STACK_SIZE_TRACING
  bigStackSize -= nbChars;
#endif
  lookStack[thread_id] = NULL ;
  if (nbChars <= nbmax) {
    curStackTop[thread_id]-=nbChars ;
    memcpy(x,curStackTop[thread_id],nbChars);
  } else {
    char *tlx = x+nbChars ;
    if (nbmax>0) memcpy(x,curStack[thread_id]->contents,nbmax) ;
    x+=nbmax ;
    while (x<tlx) {
      curStack[thread_id] = curStack[thread_id]->prev ;
      if (curStack==NULL) printf("Popping from an empty stack!!!") ;
      if (x+ONE_BLOCK_SIZE<tlx) {
	memcpy(x,curStack[thread_id]->contents,ONE_BLOCK_SIZE) ;
	x += ONE_BLOCK_SIZE ;
      } else {
	unsigned int nbtail = tlx-x ;
	curStackTop[thread_id]=(curStack[thread_id]->contents)+ONE_BLOCK_SIZE-nbtail ;
	memcpy(x,curStackTop[thread_id],nbtail) ;
	x = tlx ;
      }
    }
  }
}

/* LOOKs "nbChars" consecutive chars to a location starting at address "x".
 * Activates the LOOKing position if it was reset.
 * LOOKing is just like POPping, except that the main pointer
 * remains in place, so that the value is not POPped.
 * Further PUSHs or POPs will start from the same place as if
 * no LOOK had been made. */
void lookNarray(char *x, unsigned int nbChars) {
  int thread_id = omp_get_thread_num();
  unsigned int nbmax ;
  if (lookStack[thread_id] == NULL) {
    lookStack[thread_id] = curStack[thread_id] ;
    lookStackTop[thread_id] = curStackTop[thread_id] ;
  }
  nbmax = lookStackTop[thread_id]-(lookStack[thread_id]->contents) ;
  if (nbChars <= nbmax) {
    lookStackTop[thread_id]-=nbChars ;
    memcpy(x,lookStackTop[thread_id],nbChars);
  } else {
    char *tlx = x+nbChars ;
    if (nbmax>0) memcpy(x,lookStack[thread_id]->contents,nbmax) ;
    x+=nbmax ;
    while (x<tlx) {
      lookStack[thread_id] = lookStack[thread_id]->prev ;
      if (lookStack[thread_id]==NULL) printf("Looking into an empty stack!!!") ;
      if (x+ONE_BLOCK_SIZE<tlx) {
	memcpy(x,lookStack[thread_id]->contents,ONE_BLOCK_SIZE) ;
	x += ONE_BLOCK_SIZE ;
      } else {
	unsigned int nbtail = tlx-x ;
	lookStackTop[thread_id]=(lookStack[thread_id]->contents)+ONE_BLOCK_SIZE-nbtail ;
	memcpy(x,lookStackTop[thread_id],nbtail) ;
	x = tlx ;
      }
    }
  }
}

void resetadlookstack_() {
  int thread_id = omp_get_thread_num();
  lookStack[thread_id]=NULL ;
}

/****** Exported PUSH/POP/LOOK functions for ARRAYS: ******/

void pushcharacterarray_(char *x, unsigned int *n) {
  pushNarray(x,*n) ;
}
void popcharacterarray_(char *x, unsigned int *n) {
  popNarray(x,*n) ;
}
void lookcharacterarray_(char *x, unsigned int *n) {
  lookNarray(x,*n) ;
}

void pushbooleanarray_(char *x, unsigned int *n) {
  pushNarray(x,(*n*4)) ;
}
void popbooleanarray_(char *x, unsigned int *n) {
  popNarray(x,(*n*4)) ;
}
void lookbooleanarray_(char *x, unsigned int *n) {
  lookNarray(x,(*n*4)) ;
}

void PUSHINTEGER4ARRAY(char *x, unsigned int *n) {
  pushNarray(x,(*n*4)) ;
}
void POPINTEGER4ARRAY(char *x, unsigned int *n) {
  popNarray(x,(*n*4)) ;
}
void lookinteger4array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushinteger8array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*8)) ;
}
void popinteger8array_(char *x, unsigned int *n) {
  popNarray(x,(*n*8)) ;
}
void lookinteger8array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushinteger16array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*16)) ;
}
void popinteger16array_(char *x, unsigned int *n) {
  popNarray(x,(*n*16)) ;
}
void lookinteger16array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushreal4array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*4)) ;
}
void popreal4array_(char *x, unsigned int *n) {
  popNarray(x,(*n*4)) ;
}
void lookreal4array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*4)) ;
}

void PUSHREAL8ARRAY(char *x, unsigned int *n) {
  pushNarray(x,(*n*8)) ;
}
void POPREAL8ARRAY(char *x, unsigned int *n) {
  popNarray(x,(*n*8)) ;
}
void lookreal8array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushreal16array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*16)) ;
}
void popreal16array_(char *x, unsigned int *n) {
  popNarray(x,(*n*16)) ;
}
void lookreal16array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushreal32array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*32)) ;
}
void popreal32array_(char *x, unsigned int *n) {
  popNarray(x,(*n*32)) ;
}
void lookreal32array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*32)) ;
}

void pushcomplex4array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*4)) ;
}
void popcomplex4array_(char *x, unsigned int *n) {
  popNarray(x,(*n*4)) ;
}
void lookcomplex4array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*4)) ;
}

void pushcomplex8array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*8)) ;
}
void popcomplex8array_(char *x, unsigned int *n) {
  popNarray(x,(*n*8)) ;
}
void lookcomplex8array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*8)) ;
}

void pushcomplex16array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*16)) ;
}
void popcomplex16array_(char *x, unsigned int *n) {
  popNarray(x,(*n*16)) ;
}
void lookcomplex16array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*16)) ;
}

void pushcomplex32array_(char *x, unsigned int *n) {
  pushNarray(x,(*n*32)) ;
}
void popcomplex32array_(char *x, unsigned int *n) {
  popNarray(x,(*n*32)) ;
}
void lookcomplex32array_(char *x, unsigned int *n) {
  lookNarray(x,(*n*32)) ;
}

/****** Exported PUSH/POP/LOOK functions for F95 POINTERS: ******/

/* IMPORTANT: Don't forget to add the following interface into each calling routines:

      INTERFACE
         SUBROUTINE PUSHPOINTER(pp)
           REAL, POINTER :: pp
         END SUBROUTINE PUSHPOINTER
         SUBROUTINE POPPOINTER(pp)
           REAL, POINTER :: pp
         END SUBROUTINE POPPOINTER
      END INTERFACE

*/

void pushpointer_(char *ppp) {
  pushNarray(ppp, 4) ;
}

void poppointer_(char *ppp) {
  popNarray(ppp, 4) ;
}


/************* Debug displays of the state of the stack: ***********/

void printbigbytes(long int nbblocks, long int blocksz, long int nbunits) {
  long int a3, b3, res3, res6, res9, res12 ;
  int a0, b0, res0 ;
  int printzeros = 0 ;
  a0 = (int)nbblocks%1000 ;
  a3 = nbblocks/1000 ;
  b0 = (int)blocksz%1000 ;
  b3 = blocksz/1000 ;
  res0 = ((int)(nbunits%1000)) + a0*b0 ;
  res3 = nbunits/1000 + a3*b0 + a0*b3 ;
  res6 = a3*b3 ;
  res3 += ((long int)(res0/1000)) ;
  res0 = res0%1000 ;
  res6 += res3/1000 ;
  res3 = res3%1000 ;
  res9 = res6/1000 ;
  res6 = res6%1000 ;
  res12 = res9/1000 ;
  res9 = res9%1000 ;
  if (res12>0) {
    printf("%li ", res12) ;
    printzeros = 1 ;
  }
  if ((res9/100)>0 || printzeros) {
    printf("%li",res9/100) ;
    printzeros = 1 ;
    res9 = res9%100 ;
  }
  if ((res9/10)>0 || printzeros) {
    printf("%li",res9/10) ;
    printzeros = 1 ;
    res9 = res9%10 ;
  }
  if (res9>0 || printzeros) {
    printf("%li ",res9) ;
    printzeros = 1 ;
  }
  if ((res6/100)>0 || printzeros) {
    printf("%li",res6/100) ;
    printzeros = 1 ;
    res6 = res6%100 ;
  }
  if ((res6/10)>0 || printzeros) {
    printf("%li",res6/10) ;
    printzeros = 1 ;
    res6 = res6%10 ;
  }
  if (res6>0 || printzeros) {
    printf("%li ",res6) ;
    printzeros = 1 ;
  }
  if ((res3/100)>0 || printzeros) {
    printf("%li",res3/100) ;
    printzeros = 1 ;
    res3 = res3%100 ;
  }
  if ((res3/10)>0 || printzeros) {
    printf("%li",res3/10) ;
    printzeros = 1 ;
    res3 = res3%10 ;
  }
  if (res3>0 || printzeros) {
    printf("%li ",res3) ;
    printzeros = 1 ;
  }
  if ((res0/100)>0 || printzeros) {
    printf("%i",res0/100) ;
    printzeros = 1 ;
    res0 = res0%100 ;
  }
  if ((res0/10)>0 || printzeros) {
    printf("%i",res0/10) ;
    printzeros = 1 ;
    res0 = res0%10 ;
  }
  printf("%i",res0) ;
}

void printctraffic_() {
  printf(" C Traffic: ") ;
  printbigbytes(mmctrafficM, 1000000, mmctraffic) ;
  printf(" bytes\n") ;
}

void printftrafficinc_(long int *mmfM, int *mmfsz, int *mmf) {
  printf(" F Traffic: ") ;
  printbigbytes(*mmfM, (long int)*mmfsz, (long int)*mmf) ;
  printf(" bytes\n") ;
}

void printtopplace_() {
    int thread_id = omp_get_thread_num();
    DoubleChainedBlock *stack = curStack[thread_id] ;
    int nbBlocks = (stack?-1:0) ;
    int remainder = 0;
    while(stack) {
	stack = stack->prev ;
	nbBlocks++ ;
    }
    if (curStack[thread_id] && curStackTop[thread_id]) remainder = curStackTop[thread_id]-(curStack[thread_id]->contents) ;
    printf(" Stack size: ") ;
    printbigbytes((long int)nbBlocks, ONE_BLOCK_SIZE, (long int)remainder) ;
    printf(" bytes\n") ;
}

void printtopplacenum_(int *n) {
    int thread_id = omp_get_thread_num();
    DoubleChainedBlock *stack = curStack[thread_id] ;
    int nbBlocks = (stack?-1:0) ;
    int remainder = 0;
    while(stack) {
	stack = stack->prev ;
	nbBlocks++ ;
    }
    if (curStack[thread_id] && curStackTop[thread_id]) remainder = curStackTop[thread_id]-(curStack[thread_id]->contents) ;
    printf(" Stack size at location %i : ", *n) ;
    printbigbytes((long int)nbBlocks, ONE_BLOCK_SIZE, (long int)remainder) ;
    printf(" bytes\n") ;
}

void printstackmax_() {
    int thread_id = omp_get_thread_num();
    DoubleChainedBlock *stack = curStack[thread_id] ;
    int nbBlocks = (stack?-2:0) ;
    int remainder = 0;
    long int totalsz ;
    while(stack) {
	stack = stack->prev ;
	nbBlocks++ ;
    }
    stack = curStack[thread_id] ;
    while(stack) {
	stack = stack->next ;
	nbBlocks++ ;
    }

    printf(" Max Stack size (%i blocks): ", nbBlocks) ;
    printbigbytes((long int)nbBlocks, ONE_BLOCK_SIZE, (long int)0) ;
    printf(" bytes\n") ;
}

void printlookingplace_() {
    int thread_id = omp_get_thread_num();
    if (lookStack[thread_id] == NULL)
	printtopplace_() ;
    else {
	DoubleChainedBlock *stack = lookStack[thread_id] ;
	int nbBlocks = (stack?-1:0) ;
	while(stack) {
	    stack = stack->prev ;
	    nbBlocks++ ;
	}
        printf(" Stack look at: ") ;
        printbigbytes((long int)nbBlocks, ONE_BLOCK_SIZE,
                      ((long int)(lookStackTop[thread_id]-(lookStack[thread_id]->contents)))) ;
        printf(" bytes\n") ;
    }
}

void showrecentcstack_() {
  int thread_id = omp_get_thread_num();
  if (curStack[thread_id] && curStackTop[thread_id]) {
    int totalNumChars = 30 ;
    DoubleChainedBlock *stack = curStack[thread_id] ;
    char *stackTop = curStackTop[thread_id] ;
    unsigned short int *st1 ;
    printf("TOP OF C STACK  : ") ;
    while (totalNumChars>0 && stackTop>(stack->contents)) {
      stackTop-- ;
      st1 = (unsigned short int *)stackTop ;
      printf("%02X,",*st1%256) ;
      totalNumChars-- ;
    }
    while (totalNumChars>0 && stack->prev) {
      printf(" || ") ;
      stack = stack->prev ;
      stackTop = (stack->contents)+ONE_BLOCK_SIZE ;
      while (totalNumChars>0 && stackTop>(stack->contents)) {
        stackTop-- ;
        st1 = (unsigned short int *)stackTop ;
        printf("%02X,",*st1%256) ;
        totalNumChars-- ;
      }
    }
    if (stack->prev || stackTop>(stack->contents))
      printf(" ...\n") ;
    else
      printf(" || BOTTOM\n") ;
  } else {
    printf("NOTHING IN C STACK.\n") ;
  }
}

void getnbblocksinstack_(int *nbblocks) {
  int thread_id = omp_get_thread_num();
  DoubleChainedBlock *stack = curStack[thread_id] ;
  *nbblocks = 0 ;
  while(stack) {
    stack = stack->prev ;
    (*nbblocks)++ ;
  }
}
