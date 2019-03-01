/*
    Author: Pierre Lindenbaum PhD @yokofakun
            Institut du Thorax - Nantes - France

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <dlfcn.h>
#include "sam_dynreadfilter.h"

#define LOG_PREFIX "[dynreadfilter]"

/** concrete implementation of SamDynReadFilter */
typedef struct DLibData_t
	{
	SamDynReadFilter base;
	/* stores the output of dlopen */
	void * dl_handle;
	/* name of the plugin */
	char * name;
	}DLibData,*DLibDataPtr;


static char* dynreadfilter_get_name(SamDynReadFilterPtr h) {
return ((DLibDataPtr)h)->name;
}


SamDynReadFilterPtr dynreadfilter_load_by_name(const char* name) {
	char* libname=NULL;
	char* funname=NULL;
	/* allocate space for DLibDataPtr */
	DLibDataPtr dlibdata = (DLibDataPtr)calloc(1,sizeof(DLibData));
	if(dlibdata == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory. Cannot load hook \"%s\".",name);
		goto fail;	
		}
	/* copy name for this hook */
	dlibdata->name = strdup(name);
	if(dlibdata->name == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory. Cannot load hook \"%s\".",name);
		goto fail;	
		}
	/* set name callback */
	dlibdata->base.get_name = dynreadfilter_get_name;

        /* alloc string to build filename for the lib*.so file and build the filename */
	libname = malloc( strlen(name) + 7);/* 'lib' '.so' */
	if (libname == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory for hook \"%s\".\n",name);
		goto fail;	
		}
	sprintf(libname,"lib%s.so",name);
	
	/* open dynamic library */
	dlibdata->dl_handle = dlopen(libname, RTLD_LAZY);
	if (dlibdata->dl_handle == 0 ) {
		const char* msg = dlerror();
		fprintf(stderr, LOG_PREFIX "Cannot open dynamic library %s \"%s\".\n",
			libname, (msg == NULL?"undefined error":msg));
		goto fail;	
		}
	free(libname);
	libname = NULL;
	
	 /* alloc string to build the functions names */
	funname = malloc(strlen(name) + 15 );
	if (funname == NULL) {
		fprintf(stderr, LOG_PREFIX "Out of memory for hook \"%s\".\n",name);
		goto fail;	
		}
        /* set callback '_initialize' */
	sprintf(funname,"%s_initialize",name);
	dlibdata->base.initialize = dlsym(dlibdata->dl_handle,funname);
	if (dlibdata->base.initialize==NULL)
		{
		fprintf(stderr,LOG_PREFIX "Cannot find function %s in %s.\n",funname,name);
		goto fail;
		}
	
        /* set callback '_accept' */
	sprintf(funname,"%s_accept",name);
	dlibdata->base.accept = dlsym(dlibdata->dl_handle,funname);
	if (dlibdata->base.accept==NULL)
		{
		fprintf(stderr, LOG_PREFIX "Cannot find function %s in %s.\n",funname,name);
		goto fail;
		}
        
        /* set callback '_dispose' */
	sprintf(funname,"%s_dispose",name);
	dlibdata->base.dispose = dlsym(dlibdata->dl_handle, funname);
	if(dlibdata->base.dispose==NULL)
		{
		fprintf(stderr, LOG_PREFIX "Cannot find function %s in %s.\n",funname,name);
		goto fail;
		}

	free(funname);
	funname=NULL;

	return (SamDynReadFilterPtr)dlibdata;

        /* something is not right */
	fail:
		free(libname);
		free(funname);
		dynreadfilter_dispose_all((SamDynReadFilterPtr)dlibdata);
		
	        return NULL;
	}


int  dynreadfilter_accept_all(SamDynReadFilterPtr root,const bam_hdr_t* header, bam1_t* b) {
	while(root != NULL)
		{
		if(!root->accept(root,header,b)) return 0;
		root = root->next;
		}
	return 1;
	}

void  dynreadfilter_dispose_all(SamDynReadFilterPtr root)
	{
        DLibDataPtr dlibdata;
	if (root == NULL) return;
	if (root->next != NULL) dynreadfilter_dispose_all(root->next);
	if (root->dispose != NULL) root->dispose(root);
	dlibdata = (DLibDataPtr)root;
        if (dlibdata->dl_handle != NULL) {
		int ret = dlclose(dlibdata->dl_handle);
		if(ret != 0) {
			const char* msg = dlerror();
			fprintf(stderr, LOG_PREFIX "Cannot realease hoot %s : %s.\n",
                                  dlibdata->name,(msg == NULL?"undefined error":msg));
			}
		}
	free(dlibdata->name);
	free(dlibdata);
	}

SamDynReadFilterPtr dynreadfilter_append(SamDynReadFilterPtr root,SamDynReadFilterPtr child) {
if (root==NULL) {
  return child;
  }
else if (root->next==NULL) {
  root->next = child;
  }
else  {
  return dynreadfilter_append(root->next,child);
  }
 return root;
}
