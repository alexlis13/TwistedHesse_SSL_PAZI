#include <openssl/bn.h>
#include "twistedhesse.h"
#include <stdio.h>
#include <string.h>

int main()
{
    struct par param = {NULL, NULL, NULL, NULL};
    par_init(&param);

    //вывожу наборы параметров
    printf("Набор параметров:\n");
    printf("p = %s\n", BN_bn2dec(param.p));
    printf("u = %s\n", BN_bn2dec(param.u));
    printf("v = %s\n\n", BN_bn2dec(param.v));

    //параметры Скрученного Гессе

    printf("Параметры Скрученного Гессе:\n");
    struct twisted_hesse curve={NULL, NULL, NULL, NULL, NULL, NULL};
    twisted_hesse_init(&curve, &param);
    printf("a = %s\n", BN_bn2dec(curve.a));
    printf("d = %s\n\n", BN_bn2dec(curve.d));

    printf("Точка P:\n");
    printf("X_base = %s\n", BN_bn2dec(curve.X));
    printf("Y_base = %s\n", BN_bn2dec(curve.Y));
    printf("Z_base = %s\n\n", BN_bn2dec(curve.Z));

    struct point P = {curve.X, curve.Y, curve.Z};
    //print_in_projective(&P);
    if(aff_point_check(&P, &curve)){
        printf("Точка P находится на кривой\n\n");
    }
    else{
        printf("Точка P не находится на кривой\n\n");
    }

    printf("Нейтральный элемент:\n");
    struct point O = {NULL, NULL, NULL};
    point_init(&O ,"0", "115792089237316195423570985008687907853269984665640564039457584007913111864738", "1");
    printf("В проективных координатах:\n");
    print_in_projective(&O);
    struct point Oa = {NULL, NULL, NULL};
    point_init(&Oa ,"0", "115792089237316195423570985008687907853269984665640564039457584007913111864738", "1");
    swap_to_affin(&Oa, &O, &curve);
    printf("В афинных:\n");
    print_in_affine(&O);

    if(aff_point_check(&O, &curve)){
        printf("Точка O находится на кривой\n\n");
    }
    else{
        printf("Точка O не находится на кривой\n\n");
    }
    /*********************************************************************/

    printf("Test 1: Лежит ли точка S=(2:3:4) на кривой\n");
    struct point S = {NULL, NULL, NULL};
    point_init(&S, "2", "3", "4");
//    printf("S в аффинных координатах:\n");
//    swap_to_affin(&P2, &P2, &curve);
//    print_in_affine(&P2);
    if(aff_point_check(&S,&curve)){
        printf("Точка S находится на кривой\n\n");
    }
    else{
        printf("Точка S не находится на кривой\n\n");
    }

    /*********************************************************************/

    printf("Test 2: Проверка нахождения [k]P на кривой\n");
    BIGNUM* k = BN_new();
    BN_dec2bn(&k, "184252");
    printf("k = %s\n", BN_bn2dec(k));
    struct point kP = {BN_new(), BN_new(), BN_new()};
    cra_find(&kP, &P, &curve, k);
    //print_in_projective(&kP);
    if(aff_point_check(&kP, &curve)){
        printf("Точка kP находится на кривой\n\n");
    }
    else{
        printf("Точка kP не находится на кривой\n\n");
    }

    swap_to_affin(&kP, &kP, &curve);
    //print_in_affine(&kP);


    /*********************************************************************/

    printf("Test 3: Проверка [q]P = O\n");
    printf("q = %s\n", BN_bn2dec(param.q));
    struct point qP = {BN_new(), BN_new(), BN_new()};
    cra_find(&qP, &P, &curve, param.q);
    if(!(is_point_equal(&qP, &O,&curve))){
        printf("[q]P равно O\n\n");
    }
    else{
        printf("[q]P не равно O\n\n");
    }
//    swap_to_affin(&qP, &qP, &curve);
//    print_in_affine(&qP);

    /*********************************************************************/

    printf("Test 4: P+O, O+O\n");

    struct point p_O = {BN_new(), BN_new(), BN_new()};
    struct point o_O = {BN_new(),BN_new(), BN_new()};
    rot_sum(&P, &O, &p_O, &curve);
    rot_sum(&O, &O, &o_O, &curve);
    print_in_projective(&p_O);
    print_in_projective(&o_O);


    /*********************************************************************/

    printf("Test 5: Проверка равенства [q+1]P = P и [q-1] = -P\n\n");
    BIGNUM* num = BN_new();
    BN_dec2bn(&num, "1");                     // num = 1;
    BN_add (num, num, param.q);               // num = q + 1
    struct point result = {BN_new(),BN_new(), BN_new()};
    cra_find(&result, &P, &curve, num);
    if(!(is_point_equal(&result, &P, &curve))){
              printf("[q+1]P равно P\n\n");
     }
    else{
        printf("[q+1]P не равно P\n\n");
    }

    BN_dec2bn (&num, "1");                     // num = 1
    BN_sub (num, param.q, num);                // num = q - 1
    cra_find(&result, &P, &curve, num);
    struct point negP = {BN_new(), BN_new(), BN_new()};
    reverse_point(&negP, &P, &param);
    if(!(is_point_equal(&result, &negP, &curve))){
              printf("[q-1]P равно -P\n\n");
    }
    else{
      printf("[q+1]P не равно P\n\n");
    }

    /*********************************************************************/

    printf("Test 6: Проверка равенства [k1]P + [k2]P = [k1 + k2]P\n");
    BIGNUM* k1 = BN_new();
    BIGNUM* k2 = BN_new();
    printf("Генерация k1 и k2\n");
    BIGNUM* maxrand = BN_new();
    BN_dec2bn(&maxrand, "100000000000000000");
    BN_rand_range(k1, maxrand);
    BN_rand_range(k2, maxrand);
    printf("k1 = %s\n", BN_bn2dec(k1));
    printf("k2 = %s\n", BN_bn2dec(k2));
    BN_add(k, k1, k2);                          //k = k1 + k2
    printf("k = k1 + k2 = %s\n\n", BN_bn2dec(k));
    struct point res1 = {BN_new(), BN_new(), BN_new()};
    struct point res2 = {BN_new(), BN_new(), BN_new()};
    struct point res3 = {BN_new(), BN_new(), BN_new()};
    struct point result2 = {BN_new(),BN_new(), BN_new()};
    cra_find(&res1, &P, &curve, k1);
    cra_find(&res2, &P, &curve, k2);
    cra_find(&res3, &P, &curve, k);
    rot_sum(&res1, &res2, &result2, &curve);

    if(!(is_point_equal(&result2, &res3, &curve))){
              printf("1) [k1]P + [k2]P равно [k1 + k2]P\n\n");
     }
    else{
        printf("1) [k1]P + [k2]P не равно [k1 + k2]P\n\n");
    }
    if(aff_point_check(&res3, &curve)){
        printf("2) Точка [k]P находится на кривой\n\n");
    }
    else{
        printf("2) Точка [k]P не находится на кривой\n\n");
    }
    FreePoint(&res1);
    FreePoint(&res2);
    FreePoint(&res3);
    BN_free(k1);
    BN_free(k2);
    BN_free(k);
    FreePoint(&negP);
    BN_free(num);
    FreePoint(&result2);
    FreePoint(&S);
    FreePoint(&p_O);
    FreePoint(&o_O);
    FreePoint(&O);
    FreePoint(&P);
}
