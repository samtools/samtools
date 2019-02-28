#include <dlfcn.h>
#include "sam_hooks.h"

#define LOG_PREFIX "[hooks]"

typedef struct DLibData_t
	{
	SamHook base;
	/* stores the output of dlopen */
	void * handle;
	/* name of the plugin */
	char * name;
	}DLibData,*DLibDataPtr;


static char* hook_get_name(SamHookPtr h) {
return ((DLibDataPtr)h)->name;
}




struct sam_hook_t* hook_load_by_name(const char* name) {
	char libbame=NULL;
	char funname=NULL;
	/* allocate space for DLibDataPtr */
	DLibDataPtr *dlibdata = (DLibDataPtr*)calloc(sizeof(DLibData),0);
	if(dlibdata == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory. Cannot load hook \"%s\".",name);
		goto fail;	
		}
	/* copy name for this hook */
	dlibdata->name = strdup(name);
	if(dlibdata->name==NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory. Cannot load hook \"%s\".",name);
		goto fail;	
		}
        /* alloc string to build filename for the lib*.so file and build the filename */
	libname = malloc(strlen(name) +6);/* 'lib' '.so' */
	if(libname == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory for hook \"%s\".\n",name);
		goto fail;	
		}
	sprintf(libname,"lib%s.so",name);
	
	/* open dynamic library */
	dlibdata->handle = dlopen(libname, RTLD_LAZY);
	if(dlibdata->handle == 0 ) {
		const char* msg = dlerror();
		fprintf(stderr, LOG_PREFIX "Cannot open dynamic library %s \"%s\". Check your env variable $LD_LIBRARY_PATH ?\n",
			libname, (msg == NULL?"undefined error":msg));
		goto fail;	
		}
	free(libname);
	libname = NULL;
	
	 /* alloc string to build the functions names */
	funname = malloc(strlen(name) +12);
	if(funname == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory for hook \"%s\".\n",name);
		goto fail;	
		}
	sprintf(funname,"%s_initialize",name);
	dlibdata->base.initialize = dlsym(dlibdata->handle,funname);
	if(dlibdata->base.initialize==NULL)
		{
		fprintf(stderr,LOG_PREFIX "find function %s in %s.",funname,name);
		goto fail;
		}
	
	
	sprintf(funname,"%s_accept",name);
	dlibdata->base.accept = dlsym(dlibdata->handle,funname);
	if(dlibdata->base.accept==NULL)
		{
		fprintf(stderr, LOG_PREFIX "Find function %s in %s.",funname,name);
		goto fail;
		}

	sprintf(funname,"%s_dispose",name);
	dlibdata->base.dispose = dlsym(dlibdata->handle, funname);
	if(dlibdata->base.dispose==NULL)
		{
		fprintf(stderr, LOG_PREFIX "Find function %s in %s.",funname,name);
		goto fail;
		}

	free(funname);
	funname=NULL;

	dlibedata->base.get_name = hook_get_name;

	return (SamHookPtr)dlibdata;

	fail:
		free(libname);
		free(funname);
		DLibDataDispose(dlibdata);
		
	        return NULL;
	}


 int hook_accept_all(SamHookPtr hook,const bam_hdr_t* header, bam1_t* b) {
	while(hook!=NULL)
		{
		if(!hook->accept(hook,header,rec)) return 0;
		hook = hook->next;
		}
	return 1;
	}
void hook_dispose_all(SamHookPtr root)
	{
	int ret;
	if(hook == NULL) return;
	if(hook->next!=NULL) hook_dispose_all(hook->next);
	hook->dispose(hook);
	DLibDataPtr dlibdata = (DLibDataPtr)hook;
	ret = dlclose(dlibdata->base);
	if(ret != 0) {
		const char* msg = dlerror();
		fprintf(stderr, LOG_PREFIX "Cannot realease hoot %s : %s.\n",dlibdata->name,(msg == NULL?"undefined error":msg));
		}
	free(dlibdata->name);
	free(dlibdata);
	}

