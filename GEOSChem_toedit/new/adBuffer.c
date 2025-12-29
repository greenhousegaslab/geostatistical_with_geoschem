/*
C implementation of C$Id: adBuffer.f 3723 2011-02-24 13:34:42Z llh $

Replac the faulty PUSH/POPs by assignments to/from local arrays,
each of them being declared localized to its thread, make it
thread-safe for OPENMP.

If the OPENMP THREAD NUMBER IS GREATER THAN 256 on a single node,
the STACK_POOL_SEZE needs to be increased.

This portion just ported the necessary codes which used by WRFPLUS
from adBuffer.f, NOT THE WHOLE CODES.
Check it carefully before using for your own case.

2011-10-10 jliu@ucar.edu
*/

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
#define BUFF_SIZE 512

#define bit_set(x,y)   ((x) |  ((unsigned long)0x00000001 << (y)))
#define bit_clear(x,y) ((x) & ~((unsigned long)0x00000001 << (y)))
#define bit_test(x,y)  ((x) >> (y) & (unsigned long)0x00000001)

typedef int logical;

static int c__1 = 1;
static int c__512 = 512;
static int c_n512 = -512;
static int c__4 = 4;
static int c_n2048 = -2048;
static int c__8 = 8;
static int c_n4096 = -4096;
static int c__16 = 16;
static int c_n8192 = -8192;
static int c_b101 = 1000000;
static int c__9 = 9;
static int c__3 = 3;

static int adbitbuf[STACK_POOL_SIZE] = {0};
static int adbitibuf[STACK_POOL_SIZE] = {0};

static int *adi4buf[STACK_POOL_SIZE] = {NULL};
static int adi4ibuf[STACK_POOL_SIZE] = {0};

static double *adr8buf[STACK_POOL_SIZE] = {NULL};
static int    adr8ibuf[STACK_POOL_SIZE] = {0};

static int i4arrnum = 0;
static int r8arrnum = 0;

static int i4num = 0;
static int r8num = 0;

void print_ulong_bin(const unsigned long * const var, int bits) {
        int i;

        #if defined(__LP64__) || defined(_LP64)
                if( (bits > 64) || (bits <= 0) )
        #else
                if( (bits > 32) || (bits <= 0) )
        #endif
                return;

        for(i = 0; i < bits; i++) {
                printf("%lu", (*var >> (bits - 1 - i)) & 0x01);
        }
        printf("\n");
}

/* ========================= BITS ============================*/
int pushbit_( logical *bit )
{
  int thread_id = omp_get_thread_num();

  //printf("thread %i push %i bit %i\n", thread_id, *bit, adbitibuf[thread_id]);

  if (*bit) {
    adbitbuf[thread_id] = bit_set(adbitbuf[thread_id],
                                  adbitibuf[thread_id]);
  } else {
    adbitbuf[thread_id] = bit_clear(adbitbuf[thread_id],
                                    adbitibuf[thread_id]);
  }
  if (adbitibuf[thread_id] == 31) {
    //printf("thread %02i : pushbit : pushinteger4 %03i : ",thread_id, i4num++);
    //print_ulong_bin(&adbitbuf[thread_id],32);
    pushinteger4_(&adbitbuf[thread_id]);
    adbitbuf[thread_id]=0;
    adbitibuf[thread_id]=0;
  } else {
    adbitibuf[thread_id]++;
  }
  return 0;
}

logical popbit_(void)
{
  logical ret_val;
  int thread_id = omp_get_thread_num();
  if ( adbitibuf[thread_id] == 0 ) {
    //printf("thread %02i : popbit : popinteger4 %03i : ",thread_id, --i4num);
    //print_ulong_bin(&adbitbuf[thread_id],32);
    popinteger4_(&adbitbuf[thread_id]);
    adbitibuf[thread_id]=31;
  } else {
    adbitibuf[thread_id]--;
  }
  ret_val = bit_test(adbitbuf[thread_id], adbitibuf[thread_id]);
  //printf("thread %i pop %i bit %i\n", thread_id, ret_val, adbitibuf[thread_id]);
  return ret_val;
}

/* ========================= CONTROL ========================= */

int pushcontrol1b_(int *cc)
{
  logical L;
  L = *cc != 0;
  pushbit_(&L);
  return 0;
}

int popcontrol1b_(int *cc)
{
  if (popbit_()) {
    *cc = 1;
  } else {
    *cc = 0;
  }
  return 0;
}

int pushcontrol2b_(int *cc)
{
  logical L;
  L = bit_test(*cc,0);
  pushbit_(&L);
  L = bit_test(*cc,1);
  pushbit_(&L);
  return 0;
}

int popcontrol2b_(int *cc)
{
  if (popbit_()) {
    *cc = 2;
  } else {
    *cc = 0;
  }
  if (popbit_()) {
    *cc = bit_set(*cc,0);
  }
  return 0;
}

int pushcontrol3b_(int *cc)
{
  logical L;
  L = bit_test(*cc,0);
  pushbit_(&L);
  L = bit_test(*cc,1);
  pushbit_(&L);
  L = bit_test(*cc,2);
  pushbit_(&L);
  return 0;
}

int popcontrol3b_(int *cc)
{
  if (popbit_()) {
    *cc = 4;
  } else {
    *cc = 0;
  }
  if (popbit_()) {
    *cc = bit_set(*cc,1);
  }
  if (popbit_()) {
    *cc = bit_set(*cc,0);
  }
  return 0;
}

int pushcontrol4b_(int *cc)
{
  logical L;
  L = bit_test(*cc,0);
  pushbit_(&L);
  L = bit_test(*cc,1);
  pushbit_(&L);
  L = bit_test(*cc,2);
  pushbit_(&L);
  L = bit_test(*cc,3);
  pushbit_(&L);
  return 0;
}

int popcontrol4b_(int *cc)
{
  if (popbit_()) {
    *cc = 8;
  } else {
    *cc = 0;
  }
  if (popbit_()) {
    *cc = bit_set(*cc,2);
  }
  if (popbit_()) {
    *cc = bit_set(*cc,1);
  }
  if (popbit_()) {
    *cc = bit_set(*cc,0);
  }
  return 0;
}

int pushcontrol5b_(int *cc)
{
  logical L;
  L = bit_test(*cc,0);
  pushbit_(&L);
  L = bit_test(*cc,1);
  pushbit_(&L);
  L = bit_test(*cc,2);
  pushbit_(&L);
  L = bit_test(*cc,3);
  pushbit_(&L);
  L = bit_test(*cc,4);
  pushbit_(&L);
  return 0;
}

int popcontrol5b_(int *cc)
{
  if (popbit_()) {
    *cc = 16;
  } else {
    *cc = 0;
  }
  if (popbit_()) {
    *cc = bit_set(*cc,3);
  }
  if (popbit_()) {
    *cc = bit_set(*cc,2);
  }
  if (popbit_()) {
    *cc = bit_set(*cc,1);
  }
  if (popbit_()) {
    *cc = bit_set(*cc,0);
  }
  return 0;
}

/* ======================= INTEGER*4 =========================: */

int pushinteger4_(int *x)
{
  int thread_id = omp_get_thread_num();

  if ( adi4buf[thread_id] == NULL ) {
    int *buf = (int*)calloc(BUFF_SIZE, sizeof(int));
    if ( buf == NULL ) {
      puts("Memory allocation failed.");
      exit(0);
    }
    adi4buf[thread_id] = buf;
  }
  adi4buf[thread_id][adi4ibuf[thread_id]] = *x;
  ++adi4ibuf[thread_id];
  if ( adi4ibuf[thread_id] == BUFF_SIZE) {
     //printf("thread %02i : pushinteger4 : pushinteger4array : %04i : ", thread_id, ++i4arrnum);
     //printf("[0]: %i [%i] : %i\n",adi4buf[thread_id][0], adi4ibuf[thread_id]-1, adi4buf[thread_id][adi4ibuf[thread_id]-1]);
     pushinteger4array_(adi4buf[thread_id], &c__512);
     adi4ibuf[thread_id] = 0  ;
  }
  //printf("thread %i pushinteger4 %i buff %i\n ", thread_id, *x, adi4ibuf[thread_id]);
  return 0;
}

int popinteger4_(int *x)
{
  int thread_id = omp_get_thread_num();

  if ( adi4ibuf[thread_id] == 0 ) {
    popinteger4array_(adi4buf[thread_id], &c__512);
    adi4ibuf[thread_id] = BUFF_SIZE;
    //printf("thread %02i : popinteger4  : popinteger4array  : %04i : ", thread_id, i4arrnum--);
    //printf("[0]: %i [%i] : %i\n",adi4buf[thread_id][0], adi4ibuf[thread_id]-1, adi4buf[thread_id][adi4ibuf[thread_id]-1]);
  }
  --adi4ibuf[thread_id];
  *x = adi4buf[thread_id][adi4ibuf[thread_id]];
  //printf("thread %i popinteger4 %i buff %i\n", thread_id, *x, adi4ibuf[thread_id]);
  return 0;
}

/* ======================= REAL*8 ========================= */
int pushreal8_(double *x)
{
  int thread_id = omp_get_thread_num();

  if (adr8buf[thread_id] == NULL ) {
    double  *buf = (double*)calloc(BUFF_SIZE, sizeof(double));
    if ( buf == NULL ) {
      puts("Memory allocation failed.");
      exit(0);
    }
    adr8buf[thread_id] = buf;
  }
  adr8buf[thread_id][adr8ibuf[thread_id]] = *x;
  ++adr8ibuf[thread_id];
  if ( adr8ibuf[thread_id] == BUFF_SIZE) {

     //printf("thread %02i : pushreal8    : pushreal8array    : %04i : ", thread_id, ++r8arrnum);
     //printf("[0]: %f [%i] : %f\n",adr8buf[thread_id][0], adr8ibuf[thread_id]-1, adr8buf[thread_id][adr8ibuf[thread_id]-1]);

     pushreal8array_(adr8buf[thread_id], &c__512);
     adr8ibuf[thread_id] = 0 ;
  }
  //printf("thread %i pushreal8 %f buff %i\n", thread_id, *x,adr8ibuf[thread_id] );
  return 0;
}

int popreal8_(double *x)
{
  int thread_id = omp_get_thread_num();

  if ( adr8ibuf[thread_id] == 0 ) {
    popreal8array_(adr8buf[thread_id], &c__512);
    adr8ibuf[thread_id] = BUFF_SIZE;

   //printf("thread %02i : popreal8     : popreal8array     : %04i : ", thread_id, r8arrnum--);
   //printf("[0]: %f [%i] : %f\n",adr8buf[thread_id][0], adr8ibuf[thread_id]-1, adr8buf[thread_id][adr8ibuf[thread_id]-1]);

  }
  --adr8ibuf[thread_id];
  *x = adr8buf[thread_id][adr8ibuf[thread_id]];
  //printf("thread %i popreal8 %f buff %i\n", thread_id, *x, adr8ibuf[thread_id]);
  return 0;
}

