#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define N	(100)


int main(void)
{
	unsigned long int	record;
	unsigned long int	ctr;
	int			c;
	unsigned long int	n;

	char* p = malloc(N);
	char* word = malloc(0);
	n = N;

	record = 0;
	ctr = 0;

	while ((c = getchar()) != EOF) {

		if (ctr == n) {
			n *= 2;
			p = realloc(p, n);
		}

		if (c == ' ' || c == '\n') {

			if (ctr > record) {
				record = ctr;
				free(word);
				word = realloc(p, record);
				n = N;
				p = malloc(n);
			}
			ctr = 0;
		} else {

			if (isalpha(c)) {
				p[ctr] = c;
				ctr++;
			} else {

				while (!(c == ' ' || c == '\n')) {
					c = getchar();
				}
				ctr = 0;
			}

		}
	}


	printf("%lu characters in longest word: %s\n", record, word);

	free(word);
	free(p);

}
