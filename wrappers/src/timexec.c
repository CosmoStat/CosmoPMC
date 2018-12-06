#include <errno.h>

#include "timexec.h"


/* generic code */


/* Fork the program, create and link the pipes and launch the child loop */ 
void launchMe(forked_state *self, error **err)
{
   int errV;
   int pip[2],pop[2];
   pid_t pid;

   /* Close old pipes */
   if (self->pin>=0)
     close(self->pin);
   if (self->pout>=0)
     close(self->pout);
 
   /* Create communication pipes */
   errV = pipe(pip);
   if (errV!=0) {
      /* Cannot create pip */
      *err = addError(te_pipeerror, "Cannot create input pipe", *err, __LINE__);
      return;
   }
   errV = pipe(pop);
   if (errV!=0) {
      /* Cannot create pop */
      *err = addError(te_pipeerror, "Cannot create output pipe", *err, __LINE__);
      return;
   }
  
   self->pin  = pip[1];
   self->pout = pop[0];
  
   errno = 0;
   pid   = fork();
  
   if (pid<0) {
      /* Error */
      perror("Call to fork() failed:");
      *err = addErrorVA(te_forkerror, "Call to fork() failed, errno=%d", *err, __LINE__, errno);
      return;
   }
   if (pid==0) {
      /* Child */
    
      /* Close parent's end of the pipe */
      close(pip[1]);
      close(pop[0]);

      /* Duplicate stdin and stdout */
      self->pin=pip[0];
      if (self->attach & ATTACH_IN) {
	 dup2(pip[0],0);
	 close(pip[0]);   
	 self->pin=0;
      }
      self->pout=pop[1];
      if (self->attach & ATTACH_OUT) {
	 dup2(pop[1],1);
	 close(pop[1]);
	 self->pout=1;

      }

#ifdef COMM_DEBUG
      fprintf(stderr, "Calling own_loop %p\n", self->own_loop);
#endif

      self->own_loop(self, err);
      quitOnError(*err, __LINE__, stderr);
      *err = addError(te_forkerror, "Should not reach this line. Loop function returned prematurely(?)", *err, __LINE__);
      return;
   }
  
   /* Parent */
  
   /* Close child end of the pipe */
   close(pip[0]);
   close(pop[1]);

   self->pid     = pid;
   self->running = 1;
  
   if (alive(self)) {
      self->running = 1;
      if (self->parent_init!=NULL) {
	 self->running = self->parent_init((void*)self,err);
	 forwardError(*err,__LINE__,);
      }
      return;
   }

   *err = addError(te_forkerror, "Should not reach this line", *err, __LINE__);
   return;
}


void free_forked(forked_state **self)
{
   killIt(*self);
   if ((*self)->own_free !=NULL) {
      (*self)->own_free(&((*self)->specific_data)); 
   } else {
      free((*self)->specific_data);
   }
   free((*self)->name);
   free(*self);
   *self=NULL;
}

forked_state* init_forked(child_loop *loop, void *specific_data, long timout, free_specific *fri,
			  parent_end *parent_init, char *name, int attach, error **err)
{
   forked_state * self;

#ifdef COMM_DEBUG
   fprintf(stderr, "init_forked start\n");
#endif
  
   /* Allocate memory */
   self = (forked_state *)malloc_err(sizeof(forked_state), err);
   forwardError(*err, __LINE__, NULL);

   self->name = malloc_err(sizeof(char)*(strlen(name)+1), err);
   forwardError(*err, __LINE__, NULL);
   strcpy(self->name,name);
  
   /* Set data */
   self->attach        = attach;
   self->specific_data = specific_data;
   self->own_loop      = loop;
   self->own_free      = fri;
   self->parent_init   = parent_init;
   self->timout        = timout;
   self->pin           = -1;
   self->pout          = -1;

   launchMe(self, err);
   forwardError(*err,__LINE__,NULL);

   return self;
}
  
int alive(forked_state * self) {
   pid_t wid;
   int sts;
  
   if (self->running!=1) {
      return 1==0;
   }
   wid=waitpid(self->pid,&sts,WNOHANG);
   if (wid == self->pid) { 
      /* not alive ! */
      self->running=0;
      return 1==0;
   }
   return 1==1;
} 


void killIt(forked_state *self) {
   //if (alive(self)) {
   kill(self->pid,SIGKILL);
   self->running=0;
   //}
   self->running=0;
}

int sendAndCheck(forked_state *self, char *buf, int lenbuf, long timout)
{
   int errV;
   struct pollfd fds;  
   long ttim;

   if (!alive(self)) {
      self->running=0;
#ifdef COMM_DEBUG
      fprintf(stderr, "Process is dead, returning with error\n");
#endif
      return DEAD;
   }

   fds.fd     = self->pout;
   fds.events = POLLIN;
   ttim       = timout;

   if ((timout<=0) || (self->timout<timout)) {
      ttim = self->timout;
   }
    
   if (lenbuf>0) {
      errV = write(self->pin, buf, lenbuf);

#ifdef COMM_DEBUG
      fprintf(stderr, "*** Wrote: lenght=%d (errV=%d)\n", lenbuf, errV);
#endif

      if (errV!=lenbuf) {
	 killIt(self);
	 return WRITEERROR;
      }
   }

   errno = 0;
   ttim = -1;
#ifdef COMM_DEBUG
   fprintf(stderr, "Calling 'poll', fds=(%d %d %d) ttim=%ld\n", fds.fd, fds.events, fds.revents, ttim); fflush(stderr);
#endif

   errV  = poll(&fds,1,ttim);

#ifdef COMM_DEBUG
   fprintf(stderr, "'poll' returned, fds=(%d %d %d), errV=%d (%s)\n", fds.fd, fds.events, fds.revents,
	   errV, errV==1? "succes" : "failure");
   fflush(stderr);
#endif


   if (errno) {
      perror("Error from poll");
      errno = 0;
   }

   if (errV==0) {
      /* Time-out ! */
      killIt(self);
      return TIMEOUT;
   } else if (errV<0) {
      killIt(self);
      return READERROR;
   } else {
      /* Got something or HANGUP, or ERROR */
      if (fds.revents & POLLHUP) {
#ifdef COMM_DEBUG
	 fprintf(stderr, "POLLHUP error\n");
#endif
	 killIt(self);
	 return HANGUP;
      }
      if (fds.revents & POLLERR) {
#ifdef COMM_DEBUG
	 fprintf(stderr, "POLLERR error\n");
#endif
	 killIt(self);
	 return READERROR;
      }

#ifdef COMM_DEBUG
      //fprintf(stderr, "Returning successfully\n");
#endif
      /* Success */
      return 0;
   }  
}

void pipe_cleanup(int fd) {
   struct pollfd fds;  
   char buf[2];
   fds.fd=fd;
   fds.events=POLLIN;
   while(poll(&fds,1,0)==1) {
      if (fds.revents & (POLLHUP | POLLERR)) {
	 break;
      }
      read(fd,buf,1);
   }
}

/* specific spawn executable */
void free_exec(void** data) {
   exec_data **self;
   int i;
  
   self=(exec_data**) data;
   free((*self)->arg0);
   i=0;
   while ((*self)->argv[i]!=NULL) {
      free((*self)->argv[i]);
      i++;
   }
   free((*self)->argv);
   free(*self);
   *self=NULL;
}

forked_state * init_exec(int argc, char *arg0, char** argv, long timout, parent_end *parent_init, char* name, error** err)
{
   int i;
   exec_data* self_data;
   forked_state *res;

#ifdef COMM_DEBUG
   fprintf(stderr, "init_exec start\n");
#endif
  
   /* Allocate memory */
   self_data = (exec_data *) malloc(sizeof(exec_data));
   if (self_data==NULL) {
      *err = addError(te_allocate,"Cannot allocate memory",*err,__LINE__);
      return NULL;
   }

   /* copy args */
   self_data->arg0 = malloc_err(sizeof(char)*(strlen(arg0)+1), err);
   forwardError(*err, __LINE__, NULL);

   strcpy(self_data->arg0,arg0);
   self_data->argv = malloc_err(sizeof(char*)*(argc+1), err);
   forwardError(*err, __LINE__, NULL);

   for (i=0;i<argc;i++){
      self_data->argv[i] = malloc_err(sizeof(char)*(strlen(argv[i])+1), err);
      forwardError(*err, __LINE__, NULL);
      strcpy(self_data->argv[i],argv[i]);
   }

   self_data->argv[argc] = NULL;
   res = init_forked(exec_loop,self_data,timout,free_exec,parent_init,name,ATTACH_IN | ATTACH_OUT,err);
   forwardError(*err,__LINE__,NULL);

   return res;
}

void exec_loop(void* data, error** err)
{
   exec_data *selfexec;
   forked_state *self;

   self      = (forked_state*) data;
   selfexec  = (exec_data*) (self->specific_data);

#ifdef COMM_DEBUG
   fprintf(stderr, "exec_loop start: (%s) (%s)\n", selfexec->arg0, selfexec->argv[0]);
#endif

   execvp(selfexec->arg0, selfexec->argv);

   /* If execvp, an error had occured */
   *err = addError(te_execvp, "Process (execvp) returned unexpectedly.\n"
		              "Make sure 'scamb' is executable and can be found\n."
			      "E.g., is '.' in the search path?", *err, __LINE__);
   perror("*** perror");
   return;
}
