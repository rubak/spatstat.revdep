#ifndef AUCTIONBFNEW_H
#define AUCTIONBFNEW_H

// variant of auctionbf that also returns obj_to_pers
void auctionbf2(int *desirem, int *nn, int *pers_to_obj, int *obj_to_pers,
                double *price, double *profit, int *kk, double *eps);

// bf for back and forth: we alternate roles of bidders and objects
void auctionbf(int *desirem, int *nn, int *pers_to_obj,
               double *price, double *profit, int *kk, double *eps);

#endif

