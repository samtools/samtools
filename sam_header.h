#ifndef __SAM_HEADER_H__
#define __SAM_HEADER_H__

#ifdef __cplusplus
extern "C" {
#endif

	void *sam_header_parse2(const char *headerText);
	void *sam_header_merge(int n, const void **dicts);
	void sam_header_free(void *header);
	char *sam_header_write(const void *headerDict);   // returns a newly allocated string

    /*
        // Usage example 
        const char *key, *val; 
        void *iter = sam_header_parse2(bam->header->text);
        while ( iter = sam_header_key_val(iter, "RG","ID","SM" &key,&val) ) printf("%s\t%s\n", key,val);
    */
    void *sam_header2key_val(void *iter, const char type[2], const char key_tag[2], const char value_tag[2], const char **key, const char **value);
	char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n);

	void *sam_header2tbl(const void *dict, char type[2], char key_tag[2], char value_tag[2]);
	const char *sam_tbl_get(void *h, const char *key);
	int sam_tbl_size(void *h);
	void sam_tbl_destroy(void *h);

#ifdef __cplusplus
}
#endif

#endif
