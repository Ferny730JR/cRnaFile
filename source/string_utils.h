#ifndef STRING_UTILS_H
#define STRING_UTILS_H


/**
 *  @brief Appends the contents of the second string to the end of the first string.
 * 
 *  The function allocates memory for the combined result and updates the first string accordingly.
 *  If the first string is null or empty, it essentially duplicates the contents of `s2` into `s1`.
 *  Similarly, if `s2` is null or empty, then the contents of `s1` will remain the same.
 *
 *  @param s1 The pointer to the first string (modifiable).
 *  @param s2 The second string to append.
 */
void append(char **s1, const char *s2);


/**
 *  @brief Clean a sequence string.
 * 
 *  This function removes trailing newline character returned from rnaf_get(), capitalizes every
 *  nucleotide, and substitutes 'T' and 't' characters with 'U' if specified.
 * 
 *  @param sequence The null-terminated string to clean
 *  @param do_substitute Substitute 'T' and 't' characters with 'U'
*/
void clean_seq(char *sequence, int do_substitute);


/**
 *  @brief Capitalizes every lower case letter in string.
 * 
 *  @param str  String to capitalize
*/
void str_to_upper(char *str);


/**
 *  @brief Converts from DNA alphabet to RNA.
 * 
 *  The function will substitute all 'T' and 't' characters with 'U' and 'u' characters,
 *  respectively.
 * 
 *  @param sequence The sequence to substitute characters
*/
void seq_to_RNA(char *sequence);


void remove_escapes(char *sequence);


#endif
