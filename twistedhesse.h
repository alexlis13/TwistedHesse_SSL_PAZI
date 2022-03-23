#ifndef NULL
#define NULL (void*)0
#endif

#ifndef TWISTED_HESSE_H
#define TWISTED_HESSE_H
#include <openssl/bn.h>

#define p_str               "115792089237316195423570985008687907853269984665640564039457584007913111864739"
#define x_str               "44328971593885937857970623207174810055095945000614270339392047863929064377300"
#define y_str               "73987224069968535275377617159869580030126023743076722472100521420353122284142"
#define q_str               "115792089237316195423570985008687907853279740477758714817704293727807715164245"

// Значение, вычисленное на основе вышеперечисленных параметров с помощью Wolfram Mathematica
#define H_a_str           "8"
#define H_d_str           "48"

// структура для хранения параметров в виде больших чисел
struct par{
    BIGNUM* p;
    BIGNUM* u;
    BIGNUM* v;
    BIGNUM* q;
};

struct twisted_hesse{
    BIGNUM* X;
    BIGNUM* Y;
    BIGNUM* Z;
    BIGNUM* a;
    BIGNUM* d;
    BIGNUM* p;
};

// структура для хранения точек
struct point{
    BIGNUM* X;
    BIGNUM* Y;
    BIGNUM* Z;
};

void par_init(struct par* par);

void twisted_hesse_init(struct twisted_hesse* curve, struct par* par);

void print_in_affine(struct point* P);

void print_in_projective(struct point* P);

void point_init(struct point* point, char* X, char* Y, char* Z);

int is_point_equal(struct point* P1, struct point* P2, struct twisted_hesse* curve);

void rot_sum(struct point* P1, struct point* P2, struct point* P3, struct twisted_hesse* curve);

void std_sum(struct point* P1, struct point* P2, struct point* P3, struct twisted_hesse* curve);

void reverse_point (struct point *res, struct point *point, struct par* par);

int aff_point_check(struct point* Q, struct twisted_hesse* par);

void swap_to_affin(struct point* aff_point, struct point* P, struct twisted_hesse* curve);

void cra_find(struct point* kP, struct point* P, struct twisted_hesse* curve, BIGNUM* degree);

void FreePoint(struct point* P);

#endif // TWISTED_HESSE_H
