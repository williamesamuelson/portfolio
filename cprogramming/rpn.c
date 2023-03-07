#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include <string.h>

#define N       (10)

static void error(unsigned int line, int c, bool* err)
{
	char	buf[3];

	if (c == '\n')
		strcpy(buf, "\\n");
	else {
		buf[0] = c;
		buf[1] = 0;
	}

	printf("line %u: error at %s\n", line, buf);

	*err = true;
}

static bool isOp(int c) {
	return c == '+' || c == '-' || c == '*' || c == '/';
}

static bool checkErrors(int c, int d)
{
	return (!(isOp(c) || isdigit(c) || c == '\n' || c == ' ') ||
	(isOp(c) && d == 1));
}

int main(void)
{
        int             stack[N];
	int 		i;
        int             c;
        int             d;
        int             x;
        bool            num;
        bool            err;
        unsigned        line;

        x = 0;
        i = 0;
	d = 0;
        line = 1;
        num = false;
        err = false;

        while ((c = getchar()) != EOF) {

                if (err) {

                        if (c == '\n') {
                                line += 1;
                                err = 0;
                                i = 0;
				d = 0;
                        }

                        continue;
                } else if (isdigit(c)) {
                        x = x * 10 + c - '0';
                        num = true;
                        continue;
                } else if (num) {

                        if (i == N) {
                                error(line, '0' + x%10, &err);
				continue;
			}
                        else {
                                stack[d++] = x;
                                num = false;
                                x = 0;
				i++;
                        }
		}

		if (checkErrors(c, d)) {
			error(line, c, &err);
			continue;
		}

		if (isOp(c))
			d--;

		switch (c) {

			case '+' :
				stack[d-1] += stack[d];
				break;

			case '-' :
				stack[d-1] -= stack[d];
				break;

			case '*' :
				stack[d-1] *= stack[d];
				break;

			case '/' :
				if (stack[d] == 0) {
					error(line, c , &err);
					continue;
				} else {
					stack[d-1] /= stack[d];
					break;
				}

			case '\n' :
				if (d != 1) {
					printf("line %u: error at \\n\n", line);
				} else {
					printf("line %u: %d\n", line, stack[0]);
				}
				line++;
				d = 0;
				i = 0;
		}
        }
}
