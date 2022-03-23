#include <openssl/bn.h>
#include "twistedhesse.h"


void par_init(struct par* par){
    BN_dec2bn(&par->p, p_str);
    BN_dec2bn(&par->u, x_str);
    BN_dec2bn(&par->v, y_str);
    BN_dec2bn(&par->q, q_str);
}

void point_init(struct point* point, char* X, char* Y, char* Z){
    BN_dec2bn(&point->X, X);
    BN_dec2bn(&point->Y, Y);
    BN_dec2bn(&point->Z, Z);
}

void twisted_hesse_init(struct twisted_hesse* curve, struct par* par){
    BN_CTX* tmp = BN_CTX_new ();
    BIGNUM* invert = BN_new ();
    BIGNUM* buf_x = BN_new ();
    BIGNUM* buf_y = BN_new ();
    BIGNUM* res = BN_new ();
    BIGNUM* res_x = BN_new();
    BIGNUM* res_y = BN_new ();
    BN_dec2bn(&curve->X, "0");
    BN_dec2bn(&curve->Y, "0");
    BN_dec2bn(&curve->Z, "0");
    BN_dec2bn(&curve->a, H_a_str);
    BN_dec2bn(&curve->d, H_d_str);
    BN_dec2bn(&curve->p, p_str);

    /* Формула перехода от Вейерштрасса к Скрученному Гессе
    th_x = (18*d^2 + 72*x)  /  (d^3-12*d*x-108*a+24*y)
    th_y = (1-(48*y)        /  (d^3-12*d*x-108*a+24*y)) */

    BIGNUM* buf1 = BN_new ();
    BN_dec2bn(&buf1, "12");
    BIGNUM* buf2 = BN_new ();
    BN_dec2bn(&buf2, "108");
    BIGNUM* buf3 = BN_new ();
    BN_dec2bn(&buf3, "24");
    BIGNUM* buf4 = BN_new ();
    BN_dec2bn(&buf4, "18");
    BIGNUM* buf5 = BN_new ();
    BN_dec2bn(&buf5, "72");
    BIGNUM* buf6 = BN_new ();
    BN_dec2bn(&buf6, "48");
    BIGNUM* buf0 = BN_new();
    BN_dec2bn(&buf0, "0");
    BIGNUM* buf11 = BN_new();
    BN_dec2bn(&buf11, "1");

    BIGNUM* buf_s1 = BN_new (); // Для промежуточных значений
    BIGNUM* buf_s2 = BN_new ();
    BIGNUM* buf_s3 = BN_new ();
    BIGNUM* buf_s4 = BN_new ();
    BIGNUM* buf_s5 = BN_new ();
    BIGNUM* buf_s6 = BN_new ();
    BIGNUM* buf_s7 = BN_new ();

    //d*d*d - 12 * d * par -> u - 108+a + 24*par->v
    BN_mod_sqr(buf_s1, curve->d, curve->p, tmp);
    BN_mod_mul(buf_s1, buf_s1, curve->d, curve->p, tmp);
    BN_mod_mul(buf_s2, curve->d, buf1, curve->p, tmp);

    BN_mod_mul(buf_s2, curve->d, par -> u, curve->p, tmp);
    BN_mod_mul(buf_s2, buf_s2, buf1, curve->p, tmp);

    BN_mod_mul(buf_s3, curve->a, buf2, curve->p, tmp);

    BN_mod_mul(buf_s4, par->v, buf3, curve->p, tmp);

    BN_mod_sub(res, buf_s1, buf_s2, curve->p, tmp);
    BN_mod_sub(res, res, buf_s3, curve->p, tmp);
    BN_mod_add(res, res, buf_s4, curve->p, tmp);

    //18*d*d +72 * par->u

    BN_mod_mul(buf_s5, curve->d, curve->d, curve->p, tmp);
    BN_mod_mul(buf_s5, buf_s5, buf4, curve->p, tmp);

    BN_mod_mul(buf_s6, par->u, buf5, curve->p, tmp);

    BN_mod_add(buf_x, buf_s6, buf_s5, curve->p, tmp);

    //48 * par->v
    BN_mod_mul(buf_y, par->v, buf6, curve->p, tmp);


    BN_mod_inverse(invert, res, curve->p, tmp);
    BN_mod_mul(res_x, buf_x, invert, curve->p, tmp);

    BN_mod_mul(res_y, buf_y, invert, curve->p, tmp);

    BN_mod_sub(res_y, buf11, res_y, curve->p, tmp);

    BN_copy(curve->X, res_x);
    BN_copy(curve->Y, res_y);
    BN_dec2bn(&curve->Z, "1");

    BN_free(buf1);
    BN_free(buf2);
    BN_free(buf3);
    BN_free(buf4);
    BN_free(buf5);
    BN_free(buf6);
    BN_free(buf0);
    BN_free(buf11);
    BN_free(buf_s1);
    BN_free(buf_s2);
    BN_free(buf_s3);
    BN_free(buf_s4);
    BN_free(buf_s5);
    BN_free(buf_s6);
    BN_free(buf_s7);
    BN_free(invert);
    BN_free(buf_x);
    BN_free(buf_y);
    BN_free(res);
    BN_free(res_x);
    BN_free(res_y);
    BN_CTX_free(tmp);
}


int is_point_equal(struct point* P1, struct point* P2, struct twisted_hesse* curve){
    struct point affine_point1 = {BN_new(), BN_new(),BN_new()};
    struct point affine_point2 = {BN_new(), BN_new(), BN_new()};
//    point_init(&affine_point1,"0","-1","1");
//    point_init(&affine_point2,"0","-1","1");
    swap_to_affin(&affine_point1,P1,curve);
    swap_to_affin(&affine_point2,P2,curve);
    int res;
    if(!(BN_cmp(affine_point1.X, affine_point2.X) && BN_cmp(affine_point1.Y, affine_point2.Y))){

        res = 0; // точки равны
    }
    else res = -1;// не равны

    FreePoint(&affine_point1);
    FreePoint(&affine_point2);

    return res;

}

int aff_point_check(struct point* Q, struct twisted_hesse* curve){
    BN_CTX* tmp = BN_CTX_new ();
    BIGNUM* left  = BN_new ();
    BIGNUM* right = BN_new ();
    BIGNUM* buf  = BN_new ();
    BIGNUM* num3  = BN_new ();
    BN_dec2bn(&num3, "3");

    //a * X^3 + Y^3 + Z^3 = 3 * d * X * Y * Z - подставить точки и сравнить левую правую часть
    BN_mod_exp(left, Q->X, num3, curve -> p, tmp);             // left = X^3
    BN_mod_mul(left, left, curve -> a, curve -> p, tmp);        // left = a * X^3
    BN_mod_exp(buf, Q->Y, num3, curve -> p, tmp);            // Y^3
    BN_mod_add(left, left, buf, curve -> p, tmp);           // left = a*X^3+Y^3
    BN_mod_exp(buf, Q->Z, num3, curve -> p, tmp);          // Z^3
    BN_mod_add(left, left, buf, curve -> p, tmp);         // left = a*X^3+Y^3+Z^3
    /*printf(BN_bn2dec(left));
    printf("\n");*/

    BN_mod_mul(right, Q->X, Q->Y, curve -> p, tmp);           // right = X * Y
    BN_mod_mul(right, right, curve -> d, curve -> p, tmp);     // right = d * X * Y
    BN_mod_mul(right, right, Q -> Z, curve -> p, tmp);      // right = d * X * Y * Z
    BN_mod_sub(buf, right, left, curve -> p, tmp);         // вычесть из правой части левую и сравнить с нулем
    /*printf(BN_bn2dec(right));
    printf("\n");*/

    int res = BN_is_zero(buf);
    BN_CTX_free(tmp);
    BN_free(left);
    BN_free(right);
    BN_free(buf);
    BN_free(num3);

    return res;
}

void rot_sum(struct point* P1, struct point* P2, struct point* P3, struct twisted_hesse* curve){
    BN_CTX *tmp = BN_CTX_new();
    BIGNUM* A = BN_new();
    BIGNUM* B = BN_new();
    BIGNUM* C = BN_new();
    BIGNUM* D = BN_new();
    BIGNUM* E = BN_new();
    BIGNUM* F = BN_new();
    BIGNUM* AB = BN_new();
    BIGNUM* CD = BN_new();
    BIGNUM* DE = BN_new();
    BIGNUM* FA = BN_new();
    BIGNUM* FC = BN_new();
    BIGNUM* BE = BN_new();
    BIGNUM* T1 = BN_new();
    BIGNUM* T2 = BN_new();
    BIGNUM* T3 = BN_new();

    //алгоритм сложения

    BN_mod_mul(A, P1->X, P2->Z, curve->p, tmp); //A = X1 * Z2
    BN_mod_mul(B, P1->Z, P2->Z, curve->p, tmp); //B = Z1 * Z2
    BN_mod_mul(C, P1->Y, P2->X, curve->p, tmp); //C = Y1 * X2
    BN_mod_mul(D, P1->Y, P2->Y, curve->p, tmp); //D = Y1 * Y2
    BN_mod_mul(E, P1->Z, P2 ->Y, curve->p, tmp); //E = Z1 * Y2
    BN_mod_mul(F, P1->X, P2->X, curve->p, tmp); //F = X1 * X2
    BN_mod_mul(F, curve->a, F, curve->p, tmp); //F = a * X1 * X2
    BN_mod_mul(AB, A, B, curve->p, tmp); //A * B
    BN_mod_mul(CD, C, D, curve->p, tmp); //C * D
    BN_mod_mul(DE, D, E, curve->p, tmp); //D * E
    BN_mod_mul(FA, F, A, curve->p, tmp); //F * A
    BN_mod_mul(FC, F, C, curve->p, tmp); //F * C
    BN_mod_mul(BE, B, E, curve->p, tmp); //B * E
    BN_mod_sub(T1, AB, CD, curve->p, tmp); //AB - CD
    BN_mod_sub(T2, DE, FA, curve->p, tmp); //DE - FA
    BN_mod_sub(T3, FC, BE, curve->p, tmp); //FC - BE
    BN_copy(P3->X, T1);
    BN_copy(P3->Y, T2);
    BN_copy(P3->Z, T3);


    BN_CTX_free(tmp);
    BN_free(T1);
    BN_free(T2);
    BN_free(T3);
    BN_free(A);
    BN_free(B);
    BN_free(C);
    BN_free(D);
    BN_free(E);
    BN_free(F);
    BN_free(AB);
    BN_free(CD);
    BN_free(DE);
    BN_free(FA);
    BN_free(FC);
    BN_free(BE);
}

void std_sum(struct point* P1, struct point* P2, struct point* P3, struct twisted_hesse* curve){
    BN_CTX *tmp = BN_CTX_new();
    BIGNUM* A = BN_new();
    BIGNUM* B = BN_new();
    BIGNUM* C = BN_new();
    BIGNUM* D = BN_new();
    BIGNUM* T1 = BN_new();
    BIGNUM* T2 = BN_new();
    BIGNUM* T3 = BN_new();

    BN_mod_sqr(A, P1->X, curve->p, tmp); //X1^2
    BN_mod_mul(B, P2->Y, P2->Z, curve->p, tmp); //Y2*Z2
    BN_mod_mul(C, A, B, curve->p, tmp); // X1^2*Y2*Z2
    BN_mod_sqr(A, P2->X, curve->p, tmp); //X2^2
    BN_mod_mul(B, P1->Y, P1->Z ,curve->p, tmp); //Y1*Z1
    BN_mod_mul(D, A, B, curve->p, tmp); //X2^2*Y1*Z1
    BN_mod_sub(T1, C, D, curve->p, tmp); //X1^2*Y2*Z2 - X2^2*Y1*Z1

    BN_mod_sqr(A, P1->Z, curve->p, tmp); //Z1^2
    BN_mod_mul(B, P2->X, P2->Y, curve->p, tmp); //X2*Y2
    BN_mod_mul(C, A, B, curve->p, tmp); //Z1^2*X2*Y2
    BN_mod_sqr(A, P2->Z, curve->p, tmp); //Z2^2
    BN_mod_mul(B, P1->X, P1->Y ,curve->p, tmp); //X1*Y1
    BN_mod_mul(D, A, B, curve->p, tmp); //Z2^2*X1*Y1
    BN_mod_sub(T2, C, D, curve->p, tmp); //Z1^2*X2*Y2 - Z2^2*X1*Y1

    BN_mod_sqr(A, P1->Y, curve->p, tmp); //Y1^2
    BN_mod_mul(B, P2->X, P2->Z, curve->p, tmp); //X2*Z2
    BN_mod_mul(C, A, B, curve->p, tmp); //Y1^2*X2*Z2
    BN_mod_sqr(A, P2->Y, curve->p, tmp); //Y2^2
    BN_mod_mul(B, P1->X, P1->Z ,curve->p, tmp); //X1*Z1
    BN_mod_mul(D, A, B, curve->p, tmp); //Y2^2*X1*Z1
    BN_mod_sub(T3, C, D, curve->p, tmp); //Y1^2*X2*Z2 - Y2^2*X1*Z1

    BN_copy(P3->X, T1);
    BN_copy(P3->Y, T2);
    BN_copy(P3->Z, T3);

    BN_CTX_free(tmp);
    BN_free(T1);
    BN_free(T2);
    BN_free(T3);
    BN_free(A);
    BN_free(B);
    BN_free(C);
    BN_free(D);
}


void print_in_affine(struct point* P){
    printf("x:%s\ny:%s\n\n",BN_bn2dec(P->X),BN_bn2dec(P->Y));
}

void print_in_projective(struct point* P){
    printf("X:%s\nY:%s\nZ:%s\n\n",BN_bn2dec(P->X),BN_bn2dec(P->Y),BN_bn2dec(P->Z));
}

void swap_to_affin(struct point* aff_point, struct point* P, struct twisted_hesse* curve){
    BN_CTX *tmp = BN_CTX_new ();
    BIGNUM *x = BN_new();
    BIGNUM *y = BN_new();
    BIGNUM *z = BN_new();


    BN_mod_inverse(x, P->Z, curve->p, tmp);
    BN_mod_mul(x, x, P->X, curve->p, tmp);
    BN_mod_inverse(y, P->Z, curve->p, tmp);
    BN_mod_mul(y, y, P->Y, curve->p, tmp);


    BN_dec2bn(&z, "1");

    BN_copy(aff_point->X, x);
    BN_copy(aff_point->Y, y);
    BN_copy(aff_point->Z, z);

    BN_CTX_free(tmp);
    BN_free(x);
    BN_free(y);
    BN_free(z);
}

void cra_find(struct point* kP, struct point* P, struct twisted_hesse* curve, BIGNUM* degree){

    int bits = BN_num_bits(degree);
    struct point R = {BN_new(),BN_new(),BN_new()};
    struct point Q = {NULL,NULL,NULL};
    point_init(&Q, "0", "-1", "1");

    BN_copy(R.X, P->X);
    BN_copy(R.Y, P->Y);
    BN_copy(R.Z, P->Z);

    for(int i = bits - 1; i >= 0; --i){
        if (BN_is_bit_set(degree, i)){
            std_sum(&Q, &R, &Q, curve);//Q = Q + R
            rot_sum(&R, &R, &R, curve); //R = R + R
        }
        else{
            std_sum(&R, &Q, &R, curve); //R = R + Q
            rot_sum(&Q, &Q, &Q, curve); // Q = Q + Q
        }
    }

    BN_copy(kP->X, Q.X);
    BN_copy(kP->Y, Q.Y);
    BN_copy(kP->Z, Q.Z);

    FreePoint(&Q);
    FreePoint(&R);
}

void reverse_point (struct point *res, struct point *point, struct par* par){
    BN_CTX* tmp = BN_CTX_new ();
    BIGNUM* ret = BN_new ();
    BIGNUM* r = BN_new();

    BN_mod_inverse(ret, point->Y, par->p, tmp); // ret = 1/Y
    BN_mod_mul(r, point->X, ret, par->p, tmp); // r = X * ret = X / Y
    BN_copy(res->X, r);
    BN_copy(res->Y, ret);
    BN_copy(res->Z, point->Z);
}

void FreePoint(struct point* P){
    BN_free(P->X);
    BN_free(P->Y);
    BN_free(P->Z);
}
