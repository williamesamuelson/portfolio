#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "error.h"
#include "poly.h"

typedef struct term term;

struct term {
        int c;
        int deg;
        term* next;
};

struct poly_t {
        term* first;
        term* last;
};

term* new_term () {
        term* t = malloc(sizeof(term));
        t->c = 1;
        t->deg = 0;
        t->next = NULL;
        return t;
}

void insert_in_order(term* t, poly_t* poly)
{
        if (poly->first == NULL) {
                poly->first = t;
                poly->last = t;
        } else {
                term* p = poly->first;
                term* q = NULL;


                while(p != NULL) {

                        if (p->deg == t->deg) {
                                p->c += t->c;
                                free(t);
                                return;

                        } else if (t->deg > p->deg) {
                                t->next = p;
                                q->next = t;
                                return;
                        }
                        q = p;
                        p = p->next;
                }
                poly->last->next = t;
                poly->last = t;
        }
}

void print_term(int c, int deg, bool head, bool end)
{

        if (c < 0) {
                printf("- ");
        } else if (!head) {
                printf("+ ");
        }

        if (deg == 0) {
                printf("%d", abs(c));

        } else if (deg == 1) {

                if (c == 1){
                        printf("x");
                } else {
                        printf("%dx", abs(c));
                }

        } else {
                if (c == 1){
                printf("x^%d", deg);
                } else {
                        printf("%dx^%d", abs(c), deg);
                                }
                        }
        if (!end) {
                printf(" ");
        }
}

poly_t* new_poly_from_string(const char* p)
{

        poly_t* poly = malloc(sizeof(poly_t));
        poly->first = NULL;
        poly->last = NULL;
        term* t = new_term();


        int x = 0;
        bool num = false;
        bool exp = false;

        while (*p != 0) {

                if (isdigit(*p)) {
                        x = x * 10 + *p - '0';
                        num = true;

                } else if (num) {

                        if (exp) {
                                t->deg = x;
                                insert_in_order(t, poly);
                                exp = false;
                                t = new_term();

                        } else {
                                t->c *= x;
                                if (*p != 'x') {
                                        t->deg = 0;
                                        insert_in_order(t, poly);
                                        t = new_term();
                                } else if (*(++p) != '^') {
                                        t->deg = 1;
                                        insert_in_order(t, poly);
                                        t = new_term();
                                } else {
                                        //vi har ^
                                        exp = true;
                                }
                        }
                        num = false;
                        x = 0;

                } else if (*p == 'x') {
                        //ensamt x, koeff = 1
                        //denna del lite knas
                        if (*(++p) == '^') {
                                exp = true;
                        } else {
                                t->deg = 1;
                                insert_in_order(t, poly);
                                t = new_term();
                        }
                } else if (*p == '-') {
                        t->c *= -1;
                }

                p++;


        }

        //sista termen
        if(x != 0) {

                if (exp) {
                        t->deg = x;
                } else {
                        t->c *= x;
                }

                insert_in_order(t, poly);
        } else {
                free(t);
        }

        return poly;
}

void		free_poly(poly_t* poly)
{
        term* t = poly->first;
        term* s;

        free(poly);

        while(t != NULL) {
                s = t->next;
                free(t);
                t = s;
        }
}

poly_t*		mul(poly_t* p, poly_t* q)
{
        poly_t* res = malloc(sizeof(poly_t));
        res->first = NULL;
        res->last = NULL;
        term* t = p->first;
        term* s = q->first;
        term* n;

        while (t != NULL) {
                while (s != NULL) {
                        n = new_term();
                        n->c = t->c * s->c;
                        n->deg = t->deg + s->deg;
                        insert_in_order(n, res);
                        s = s->next;
                }
                s = q->first;
                t = t->next;
        }

        return res;

}

void print_poly(poly_t* poly)
{
        term* t = poly->first;

        print_term(t->c, t->deg, true, false);

        t = t->next;

        while(t->next != NULL) {

                print_term(t->c, t->deg, false, false);

                t = t->next;
        }

        print_term(t->c, t->deg, false, true);

        putchar('\n');

}
