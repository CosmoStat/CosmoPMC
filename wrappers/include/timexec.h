#ifndef __TIMEXEC__
#define __TIMEXEC__

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <poll.h>
#include <signal.h>

#include "errorlist.h"
#include "io.h"



typedef void child_loop(void*,error** err);
typedef void  free_specific(void**);
typedef int  parent_end(void*,error** err);

#define ATTACH_IN 1
#define ATTACH_OUT 2

typedef struct {
  int running; 
  int pin, pout, perr;
  int attach;
  long timout;
  pid_t pid;
  void* specific_data;
  child_loop *own_loop;
  free_specific *own_free;
  parent_end *parent_init;
  char* name;
} forked_state;
 
typedef struct {
  char* arg0;
  char** argv;
} exec_data;

#define DEAD -1000
#define TIMEOUT -1001
#define HANGUP -1002
#define WRITEERROR -1003
#define READERROR -1004

/* generic code */
void launchMe(forked_state *self,error** err);
void free_forked(forked_state **self);
forked_state* init_forked(child_loop* loop,void* specific_data,long timeout, free_specific* fri,parent_end *parent_init,char* name, int attach,error** err);
int alive(forked_state * self);
void killIt(forked_state *self);
int sendAndCheck(forked_state * self, char * buf, int lenbuf, long timout);
void pipe_cleanup(int fd);

/* specific spawn executable */
forked_state* init_exec(int argc, char *arg0, char** argv,long timout,parent_end* parent_init,char* name,error** err);
void exec_loop(void*,error** err);
void free_exec(void** data);

/* error numbers */
#define te_base      -10000
#define te_pipeerror -1 + te_base
#define te_forkerror -2 + te_base
#define te_allocate  -3 + te_base
#define te_execvp    -4 + te_base

#endif

