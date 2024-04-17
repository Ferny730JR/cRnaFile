#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "memory_utils.h"
#include "string_utils.h"

void append(char **s1, const char *s2) {
    const size_t len1 = *s1 ? strlen(*s1) : 0;
    const size_t len2 =  s2 ? strlen(s2)  : 0;

    if(len2 == 0) {
        return;     // s2 is NULL or empty, so dont modify contents of s1
    }

    *s1 = s_realloc(*s1,len1 + len2 + 1);

    if(len1 == 0) {
        strcpy(*s1, s2);    // s1 is NULL or empty, so copy contents of s2
    } else {
        strcat(*s1, s2);    // append contents of s2 onto s1
    }
}


void clean_seq(char *sequence, int do_substitute) {
	size_t ln = strlen(sequence)-1;

	if(sequence[ln] == '\n') {  // remove trailing new line character
		sequence[ln] = '\0';
	}

	for(unsigned int i = 0; sequence[i]; i++) {
		sequence[i] = toupper(sequence[i]);
		if(do_substitute && (sequence[i] == 'T' || sequence[i] == 't')) {
			sequence[i] = 'U';
		}
	}
}


void str_to_upper(char *str) {
    if(!str) {
        error_message("Unable to read string %s",str);
    }

    for(int i=0; str[i]; i++) {
        str[i] = toupper(str[i]);
    }
}


void seq_to_RNA(char *sequence) {
    unsigned int i;

    if(!sequence) {
        error_message("Unable to read string %s",sequence);
    }

    for (i = 0; sequence[i]; i++) {
        if (sequence[i] == 'T') {
            sequence[i] = 'U';
        } else if (sequence[i] == 't') {
            sequence[i] = 'u';
        }
    }
}


void remove_escapes(char *str) {

    if(!str) {  // str is NULL
        return;
    }

    size_t ln = strlen(str)-1;

    if(str[ln] == '\n') {  // remove trailing new line character
        str[ln] = '\0';
    }
    
    while(isspace(*str)) {  // move pointer past white space
        ++str;
    }
}
