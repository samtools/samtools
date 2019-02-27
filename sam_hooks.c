#include "htslib/sam.h"
#include "sam_hooks.h"

struct sam_hook_impl_t { 
	struct sam_hook_t base;
	void* handle;
	};

struct sam_hook_t* hook_load_by_name(const char* name) {
	struct sam_hook_impl_t *hook = (struct sam_hook_impl_t*)calloc(sizeof(struct sam_hook_impl_t),0);
	void * handle = dlopen("./libdog.so", RTLD_LAZY);
	
	*(void**)(&func_print_name) = dlsym(handle, "accept");
	
	return (struct sam_hook_t*)hook;
	}


