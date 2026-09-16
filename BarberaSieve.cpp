/// This is an implementation of the bit wheel segmented sieve
/// with max modulus wheel choice 30, 210, 2310

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <stdint.h>
#include <time.h>

const int64_t PrimesBase[5]={2,3,5,7,11};
const int64_t n_PB_max = 5;

const int64_t del_bit[8] =
{
~(1 << 0),~(1 << 1),~(1 << 2),~(1 << 3),
~(1 << 4),~(1 << 5),~(1 << 6),~(1 << 7)
};

const int64_t bit_count[256] =
{
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

int64_t Euclidean_Diophantine(int64_t coeff_a, int64_t coeff_b)
{
    // return y in Diophantine equation coeff_a x + coeff_b y = 1
    int64_t k=1;
    std::vector<int64_t> div_t;
    std::vector<int64_t> rem_t;
    std::vector<int64_t> coeff_t;
    div_t.push_back(coeff_a);
    rem_t.push_back(coeff_b);
    coeff_t.push_back((int64_t)0);
    div_t.push_back((int64_t)div_t[0]/rem_t[0]);
    rem_t.push_back((int64_t)div_t[0]%rem_t[0]);
    coeff_t.push_back((int64_t)0);
    while (rem_t[k]>1)
    {
        k=k+1;
        div_t.push_back((int64_t)rem_t[k-2]/rem_t[k-1]);
        rem_t.push_back((int64_t)rem_t[k-2]%rem_t[k-1]);
        coeff_t.push_back((int64_t)0);
    }
    k=k-1;
    coeff_t[k]=-div_t[k+1];
    if (k>0)
        coeff_t[k-1]=(int64_t)1;
    while (k > 1)
    {
        k=k-1;
        coeff_t[k-1]=coeff_t[k+1];
        coeff_t[k]+=(int64_t)(coeff_t[k+1]*(-div_t[k+1]));
    }
    if (k==1)
        return (int64_t)(coeff_t[k-1]+coeff_t[k]*(-div_t[k]));
    else
        return (int64_t)(coeff_t[0]);
}

void segmented_bit_sieve_wheel(uint64_t n, int64_t max_bW)
{
    int64_t sqrt_n = (int64_t) std::sqrt(n);
    int64_t count_p=(int64_t)0;
    int64_t n_PB=(int64_t)3;
    int64_t bW=(int64_t)30;
    //get bW modulus equal to p1*p2*...*pn <=max_bW with n=n_PB
    while(n_PB<n_PB_max&&(bW*PrimesBase[n_PB]<=std::min(max_bW,sqrt_n)))
    {
        bW*=PrimesBase[n_PB];
        n_PB++;
    }
    for (int64_t i=0; i< n_PB; i++)
        if (n>PrimesBase[i])
            count_p++;
    if (n>1+PrimesBase[n_PB-1]){
        int64_t k_end = (n < bW) ? (int64_t)2 :(int64_t) (n/(uint64_t)bW+1);
        int64_t k_sqrt = (int64_t) std::sqrt(k_end/bW)+1;
        //find possible remainder of the congruence class
        std::vector<char> Remainder_i_t(bW+1,true);
        for (int64_t i=0; i< n_PB; i++)
            for (int64_t j=PrimesBase[i]*PrimesBase[i];j< bW+1;j+=PrimesBase[i])
                Remainder_i_t[j]=false;
        std::vector<int64_t> RW;
        for (int64_t j=PrimesBase[n_PB-1]+1;j< bW+1;j++)
            if (Remainder_i_t[j]==true)
                RW.push_back(-bW+j);
        RW.push_back(1);
        int64_t nR=RW.size();
        std::vector<int64_t> C1(nR*nR);
        std::vector<int64_t> C2(nR*nR);
        for (int64_t j=0; j<nR-2; j++)
        {
            int64_t rW_t,rW_t1;
            rW_t1=Euclidean_Diophantine(bW,-RW[j]);
            for (int64_t i=0; i<nR; i++)
            {
                if (i==j)
                {
                    C2[nR*i+j]=0;
                    C1[nR*i+j]=RW[j]+1;
                }
                else if(i==nR-3-j )
                {
                    C2[nR*i+j]=1;
                    C1[nR*i+j]=RW[j]-1;
                }
                else
                {
                    rW_t=(int64_t)(rW_t1*(-RW[i]))%bW;
                    if (rW_t>1)
                        rW_t-=bW;
                    C1[nR*i+j]=rW_t+RW[j];
                    C2[nR*i+j]=(int64_t)(rW_t*RW[j])/bW+1;
                    if (i==nR-1)
                        C2[nR*i+j]-=1;
                }
            }
            C2[nR*j+nR-2]=(int64_t)1;
            C1[nR*j+nR-2]=-(bW+RW[j])-1;
            C1[nR*j+nR-1]=RW[j]+1;
            C2[nR*j+nR-1]=(int64_t)0;
        }
        for (int64_t i=nR-2; i<nR; i++)
        {
            C2[nR*i+nR-2]=(int64_t)0;
            C1[nR*i+nR-2]=-RW[i]-1;
            C1[nR*i+nR-1]=RW[i]+1;
            C2[nR*i+nR-1]=(int64_t)0;
        }
        int64_t nB=nR/8;
        int64_t segment_size=1;
        int64_t p_mask_i=(int64_t)4;
        for (int64_t i=0; i<p_mask_i;i++)
            segment_size*=(bW+RW[i]);
        // if bW=30 =7*11*13*17
        while (segment_size<k_sqrt && p_mask_i<7)
        {
            segment_size*=(bW+RW[p_mask_i]);
            // if bW=30 max value =7*11*13*17*19*23*29
            p_mask_i++;
        }
        int64_t segment_size_b=nB*segment_size;
        std::vector<uint8_t> Primes(nB+segment_size_b, 0xff);
        std::vector<uint8_t> Segment_i(nB+segment_size_b, 0xff);
        int64_t pb,mb,mmin,ib,i,jb,j,k,kb;
        int64_t kmax = (int64_t) std::sqrt(segment_size/bW)+(int64_t)1;
        for (k =(int64_t)1; k <= kmax; k++)
        {
            kb=k*nB;
            for (jb = 0; jb<nB; jb++)
            {
                for (j = 0; j<8; j++)
                {
                    if(Primes[kb+jb] & (1 << j))
                    {
                        for (ib = 0; ib<nB; ib++)
                        {
                            for (i = 0; i<8; i++)
                            {
                                pb=nB*(bW*k+RW[j+jb*8]);
                                mmin=nB*(bW*k*k + k*C1[(i+ib*8)*nR+j+jb*8] + C2[(i+ib*8)*nR+j+jb*8]);
                                for (mb =mmin; mb <= segment_size_b && mb>=(int64_t)0; mb +=pb )
                                    Primes[mb+ib] &= del_bit[i];
                                if (pb<nB*(bW+RW[p_mask_i]) && k_end>segment_size)
                                {
                                    mb-=segment_size_b;
                                    while (mb<(int8_t)0)
                                        mb+=pb;
                                    for (; mb <= segment_size_b; mb +=pb )
                                        Segment_i[mb+ib] &= del_bit[i];
                                }
                            }
                        }
                    }
                }
            }
        }
        for (kb = nB; kb < std::min (nB+segment_size_b,nB*k_end); kb++)
            count_p+=bit_count[Primes[kb]];
        if (kb==nB*k_end && kb<=segment_size_b && kb>(int64_t)0)
            for (ib = 0; ib<nB; ib++)
                for (i = 0; i < 8; i++)
                    if(Primes[kb+ib]& (1 << i) && RW[i+ib*8]<(int64_t)(n%bW-bW))
                        count_p++;
        if (k_end>segment_size)
        {
            int64_t k_low, kb_low;
            std::vector<uint8_t> Segment_t(nB+segment_size_b);
            for (int64_t k_low = segment_size; k_low < k_end; k_low += segment_size)
            {
                kb_low=k_low*nB;
                for (kb = (int64_t)0; kb <(nB+segment_size_b); kb++)
                    Segment_t[kb]=Segment_i[kb];
                kmax=(std::min(segment_size,(int64_t)std::sqrt((k_low+segment_size)/bW)+2));
                j=p_mask_i;
                for(k=(int64_t)1; k<=kmax;k++)
                {
                    kb=k*nB;
                    for (jb = 0; jb<nB; jb++)
                    {
                        for (; j < 8; j++)
                        {
                            if (Primes[kb+jb]& (1 << j))
                            {
                                for (ib = 0; ib<nB; ib++)
                                {
                                    for (i = 0; i < 8; i++)
                                    {
                                        pb=bW*k+RW[j+jb*8];
                                        mmin=-k_low+bW*k*k+ k*C1[(i+ib*8)*nR+j+jb*8] + C2[(i+ib*8)*nR+j+jb*8];
                                        if (mmin<0)
                                            mmin=(mmin%pb+pb)%pb;
                                        mmin*=nB;
                                        pb*=nB;
                                        for (mb =mmin; mb <= segment_size_b; mb += pb)
                                            Segment_t[mb+ib] &= del_bit[i];
                                    }
                                }
                            }
                        }
                        j=(int64_t)0;
                    }
                }
                for ( kb =nB+kb_low; kb <std::min (kb_low+segment_size_b+nB,nB*k_end); kb++)
                    count_p+=bit_count[Segment_t[kb-kb_low]];
            }
            if (kb==nB*k_end && kb-kb_low<=segment_size_b && kb-kb_low>(int64_t)0)
                for (ib = 0; ib<nB; ib++)
                    for (i = 0; i < 8; i++)
                        if(Segment_t[kb-kb_low+ib]& (1 << i) && RW[i+ib*8]<(int64_t)(n%bW-bW))
                            count_p++;
        }
    }
    std::cout << " primes < " << n << ": "<< count_p<< std::endl;
}

int main()
{
    // segmented_bit_sieve_wheel(n, max_bW) with max_bW= 30 , 210 , 2310 for set modulus
    clock_t start = clock();
    segmented_bit_sieve_wheel(1000000000,210);
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;
    std::cout << " time: " << elapsed << " seconds" << std::endl;
    return 0;
}
