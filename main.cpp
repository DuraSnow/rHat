//
//  main.cpp
//  NoisyLRA
//
//  Created by Abelard on 16/3/14.
//  Copyright © 2016年 Abelard. All rights reserved.
//

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <map>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stack>
#include "ksw.h"
#include <unistd.h>
#include <sstream>
#include <unistd.h>
#include <queue>
#define PX(X) std::cout << X << std::endl
using namespace std;
int L=2048;
int length_K = 11;
uint32_t get_c[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    3};

uint32_t MASK = 0x3fffff;


typedef pair<uint32_t, int > window;
typedef pair<int, int> pointer;
typedef struct node
{
    int ref_index, read_index;
    int len;
}NODE;
NODE *node_list;
size_t **V;
int *dp;
size_t *path;
size_t p_index;
size_t node_i;
pointer plist[4194304];
window wlist[10000000];
uint64_t *read_wlist;
int countWindow[10000],Target_w[5];
int8_t  mat[25];
int8_t match=2, mism=-5, gape=2, gapo=1;
const char correspondTable[] = "MIDNSHP=X";
uint32_t *cigar;
int n_cigar = 0;
uint8_t  readqry[2048];
uint8_t  refqry[2048];
uint8_t *readqry_;
uint8_t *refqry_;
const uint8_t seq_nt4_tablet[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 2, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

bool cmp(window a,window b)
{
    return  a.first<b.first;
   
}
bool cmp2(uint64_t a, uint64_t b)
{
    return a<b;
}
void empty(stack<int> st){
    while (!st.empty()) {
        st.pop();
    }
}
typedef struct window_cnt {
    uint32_t index_of_W;
    size_t cnt;
    bool operator<(const window_cnt& w) const
    {
        return cnt > w.cnt;
    }
}WINDOW_CNT;
void maxget(int l){
    
    priority_queue<WINDOW_CNT> Hwin_cnt;
    for (int i = 0; i < 10000; i++)
    {
        WINDOW_CNT node;
        node.index_of_W = i;
        node.cnt = countWindow[i];
        Hwin_cnt.push(node);
        if (Hwin_cnt.size() > 5)
            Hwin_cnt.pop();
    }
 
    int num = 4;
    while (!Hwin_cnt.empty())
    {
        Target_w[num--] = Hwin_cnt.top().index_of_W;
        Hwin_cnt.pop();
    }
    
}

int getwindowlist(string s,int l,window list[])
{
    int i = 0;
    int wlen = 0;
    uint32_t next = 0,current;
    while(i< l)
    {
        next = (next << 2) | get_c[s[i]];
        current = next & MASK;
        if(i >= 10)
        {
            int k = i/1024;
            if (k == 0)
            {
                list[wlen].first = current;
                list[wlen].second = k;
           
                wlen++;
            }
            else
            {
                list[wlen].first = current;
                list[wlen].second = k;
                wlen++;
                list[wlen].first = current;
                list[wlen].second = k+1;
          
                wlen++;
            }
            
        }
        i++;
        
    }

    sort(list,list+wlen,cmp);
    
    return wlen;
}
void getpointlist(int wnum, window list1[], pointer list2[] )
{
    uint32_t p = list1[0].first;
    int  i = 0;
    list2[p].first=0;

    
    while (i < wnum)
    {
        if (list1[i].first != p)
        {
            list2[p].second = i-1 ;
            p = list1[i].first;
            list2[p].first = i ;
        }
        i++;
    }
   
    
}

void findwindow(string s, int wlen)
{
    int i = 0;
    uint32_t current, next;
    while (i < 1024)
    {
        next = (next << 2) | get_c[s[i]];
        current = next & MASK;
        if(i >= 10)
        {
            if (plist[current].first == plist[current].second)
            {
                
                if (plist[current].first != 0)
                {
                    countWindow[wlist[plist[current].second].second]++;
                }
            }
            else{
                for (int j = plist[current].first; j <= plist[current].second;j++)
                {
                    
                    countWindow[wlist[j].second]++;
                }
            }
        }
        i++;
    }
}

    
size_t getreadpointer(string s)
{
    read_wlist = new uint64_t[s.length()];
    uint64_t  current,next = 0,index;
    size_t i;
    for ( i=0; i < s.length(); i++)
    {
        next = (next << 2) | get_c[s[i]];
        current = next & MASK;
        if (i >= 10)
        {
            current = current<<32;
            index = i;
            current = current|index;
            read_wlist[i] = current;
        }
    }
    sort(read_wlist,read_wlist+i-10,cmp2);
    return i;
}

void add_node(uint32_t w, uint32_t r, size_t l)
{
    node_list[node_i].ref_index = w;
    node_list[node_i].read_index = r;
    node_list[node_i].len = l;
    ++node_i;
}
bool check_node(size_t i)
{
    if (node_i > 1 && i > node_list[node_i-1].ref_index && (i - node_list[node_i-1].ref_index) == node_list[node_i-1].len - 11 + 1)
    {
        ++node_list[node_i-1].len;
        return true;
    }
    return false;
}

uint32_t search(uint32_t pindex, size_t index)
{
    //binary search
    uint32_t l=0, r=index-1, m;
    while (l < r)
    {
        m = (l+r)/2;
        if (pindex <= (read_wlist[m]>>32)) r = m;
        else l = m + 1;
    }
    if ((read_wlist[r]>>32) == pindex) return r+1;
    return 0;
}
void create_matrix()
{
    V = new size_t*[node_i];
    for (size_t i=0; i<node_i; ++i) V[i] = new size_t[node_i];
    uint32_t iw, ir, jw, jr;
    size_t il;
    
    for (size_t i=0; i<node_i-1; ++i)
    {
        iw = node_list[i].ref_index;
        ir = node_list[i].read_index;
        il = node_list[i].len;
        for (size_t j=i+1; j<node_i; ++j)
        {
            V[i][j]=0;
            jw = node_list[j].ref_index;
            jr = node_list[j].read_index;
            if (jw >= il + iw && jr >= il + ir && jr <= ir + il + 1024)
            {
                V[i][j] = 1;
            }
        }
    }
}
void getpath()
{
    dp = new int[node_i];
    size_t thenode;
    for (size_t k=0; k<node_i; ++k){
        dp[k] = 0;
    }
    for (size_t i=0; i<node_i-1; ++i)
    {
        for (size_t j=i+1; j<node_i; ++j)
        {
            if (V[i][j]) {
                
                dp[j] = (dp[j] > node_list[j].len+dp[i] ? dp[j] : node_list[j].len+dp[i]);
                thenode = j;
                
            }
        }
    }
    
   
    size_t p_tmp;
    p_index = 0;
    path = new size_t[node_i];
    while (thenode != 0)
    {
        p_tmp = thenode;
        for (size_t i = 0; i< p_tmp; ++i)
        {
            if (i==thenode || V[i][thenode]==0) continue;
            if (dp[thenode] == node_list[thenode].len + dp[i])
            {
                thenode = i;
                //PX(i);
                path[p_index++] = i;
                break;
            }
        }
        if(thenode == p_tmp) break;
    }
}
inline void transIntoDec(uint8_t *transtr,char *str, int length)
{
    for (int i=0;i<length;++i) {
        transtr[i] = seq_nt4_tablet[str[i]];
    }
}

double do_alignment(char* ref, size_t window_up, size_t window_down, char*  read, size_t read_l, string& thecigar, FILE* outt)
{
    
    int k=0;
    for (int i=0; i<5; ++i)
    {
        for (int j=0; j<5; ++j)
        {
            if (i<4 && j<4) mat[k++] = i == j? match : mism;
            else mat[k++]=0;
        }
    }
    
    
    
    size_t PointerListLen =11;
    
    int read_len;
    int ref_len;
    int w_;
    
    uint32_t countM = 0;
    char     trans_cigar[500];
    int startPosCigar = 0;
    
    int qlen = 0;
    int tlen = 0;
    
    
    double score=0;
    
    size_t last_w=window_up + PointerListLen - 1, last_r=PointerListLen-1;
    size_t i, w, r, l;
    
    if (p_index < 3) return 0;
    if (p_index >= 3)
    {
        i = path[p_index-2];
        w = node_list[i].ref_index;
        r = node_list[i].read_index;
        l = node_list[i].len;
        ref_len = w-last_w;
        read_len = r-last_r;
        
        if (ref_len > PointerListLen / 2)
        {
            size_t tmp = ref_len - PointerListLen/2;
            ref_len = PointerListLen / 2;
            last_w += tmp;
        }
        if (read_len > PointerListLen / 2)
        {
            size_t tmp = read_len - PointerListLen/2;
            read_len = PointerListLen / 2;
            last_r += tmp;
        }
        if (w!=last_w && r!=last_r)
        {
            transIntoDec(readqry,read + last_r - PointerListLen + 1,read_len);
            transIntoDec(refqry,ref + last_w - PointerListLen + 1,ref_len);
            readqry_ = readqry;
            refqry_ = refqry;
            w_ = max(ref_len, read_len);
            n_cigar = 0;
            score += ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
            if (n_cigar) {
                for (int z=0;z<n_cigar-1;++z) {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
            
                    thecigar.append(trans_cigar);
                }
                
                if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                    countM = cigar[n_cigar-1] >> 4;
                } else {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                    thecigar.append(trans_cigar);
                }
                
            }
            free(cigar);
        }
        score+=(l)*match;
        countM+=l;
        last_w = w-PointerListLen+l+1;
        last_r = r-PointerListLen+l+1;
    }
    
    for (size_t o = p_index-2; o>0; --o)
    {
        if (o == p_index) continue;
        i = path[o-1];
        
        w = node_list[i].ref_index;
        r = node_list[i].read_index;
        l = node_list[i].len;
        
        
        
        ref_len = w-last_w-PointerListLen+1;
        read_len = r-last_r-PointerListLen+1;
        
             if (w>last_w+PointerListLen-1|| r>last_r+PointerListLen-1)
        {
            if (ref_len > 2048 || read_len > 2048)
            {
                return -INFINITY;
                cout << "hehe\t" << ref_len << " " << read_len << endl;
            }
            
            transIntoDec(readqry,read + last_r,read_len);
            transIntoDec(refqry,ref + last_w, ref_len);
            readqry_ = readqry;
            refqry_ = refqry;
            w_ = max(ref_len, read_len);
            n_cigar = 0;
            score += ksw_global(read_len,readqry_,ref_len,refqry_,5,mat,gapo,gape,w_,&n_cigar,&cigar);
            if (n_cigar - 1)
            {
                if (correspondTable[cigar[0]&0xf] == 'M') {
                    countM += (cigar[0] >> 4);
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                }
                else
                {
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                    thecigar.append(trans_cigar);
                }
                countM = 0;
                for (int z=1;z<n_cigar-1;++z) {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                    thecigar.append(trans_cigar);
                }
                
                if (correspondTable[cigar[n_cigar-1]&0xf] == 'M') {
                    countM = cigar[n_cigar-1] >> 4;
                } else {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[n_cigar-1]>>4,correspondTable[cigar[n_cigar-1]&0xf]);
                    thecigar.append(trans_cigar);
                }
                
            }
            else if (n_cigar>0)
            {
                if (correspondTable[cigar[0]&0xf] == 'M')
                    countM += (cigar[0] >> 4);
                else {
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                    thecigar.append(trans_cigar);
                    countM = 0;
                    
                }
            }
            free(cigar);
        }
        
        score+=(l)*match;

        countM+=l;
      
        last_w = w-PointerListLen+l+1;
        last_r = r-PointerListLen+l+1;
    }
    
  
    if (p_index >= 3)
    {
        w = window_down;
        r = read_l;
        
        ref_len = w-last_w;
        read_len = r-last_r;
        
        if (ref_len > PointerListLen/2) ref_len = PointerListLen / 2;
        if (read_len > PointerListLen/2) read_len = PointerListLen / 2;
        
        if (w>last_w && r>last_r)
        {
            transIntoDec(readqry,read + last_r,read_len);
            transIntoDec(refqry,ref + last_w, ref_len);
            readqry_ = readqry;
            refqry_ = refqry;
            w_ = max(ref_len, read_len);
            n_cigar = 0;
            qlen = 0;
            tlen = 0;
            score += ksw_extend_core(read_len, readqry_, ref_len, refqry_, 5, mat, gapo, gape, 40, read_len , &qlen, &tlen, &cigar, &n_cigar) - read_len;
            if (n_cigar) {
                if (correspondTable[cigar[0]&0xf] == 'M') {
                    countM += (cigar[0] >> 4);
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                }
                else
                {
                    startPosCigar += sprintf(trans_cigar,"%uM",countM);
                    thecigar.append(trans_cigar);
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[0]>>4,correspondTable[cigar[0]&0xf]);
                    thecigar.append(trans_cigar);
                }
                countM = 0;
                for (int z=1;z<n_cigar;++z) {
                    startPosCigar += sprintf(trans_cigar,"%u%c",cigar[z]>>4,correspondTable[cigar[z]&0xf]);
                    thecigar.append(trans_cigar);
                }
                
            }
            free(cigar);
        }
    }
        return score;
}
int main(int argc, const char * argv[])
{

   
    
    ifstream f;
    f.open("/Users/abelard/Desktop/File/E.coli.fa");
    string buffer,s="";
    getline(f,buffer);
    while(getline(f,buffer))
    {
        s +=buffer;
    }
    f.close();
    int len =s.length(),wnum;
     const char *c_dna_f = s.c_str();
    memset(wlist, 0, sizeof(wlist));
    memset(plist, 0, sizeof(plist));
    wnum = getwindowlist(s,len,wlist);
    getpointlist(wnum,wlist,plist);
    printf("Time used = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);
    
    char sss[10000];
    char thesss[10000];
    
    ifstream f2;
    f2.open("/Users/abelard/Desktop/File/E.coli-sim.fastq.coli-sim");
    string rs="";
    FILE *out;
    out = fopen("/Users/abelard/Desktop/File/out", "wb");
    while( getline(f2,buffer))
    {
        
        
        getline(f2, rs);
        getline(f2,buffer);
        
        getline(f2,buffer);
        
        
    
        
        
        const char *c_read = rs.c_str();
        int rs_len = rs.length()/2,len_up_stream,len_down_stream;
        int rs_hflen = rs_len - 512;
        len_up_stream = rs_hflen;
        len_down_stream =rs_len - len_up_stream - 1024;
        
        
        findwindow(rs.substr(rs_hflen,1024), wnum);
       
        memset(Target_w,0,sizeof(Target_w));
        
        maxget(10000);
        
        size_t index_rs = getreadpointer(rs);
        size_t read_len = rs.length();
        
        double t;
        
       
        double score = -INFINITY, score_enough=100;
        double score_too_low=-100;
        string thess(read_len+2048, 0);
       

        for(int  i = 0; i < 5; i++ )
        {
            size_t window_up, window_down;
            window_up = ( Target_w[i] * (L / 2) >= len_up_stream ) ? ( Target_w[i] * (L / 2) - len_up_stream ) : (0);
            window_down = ( Target_w[i] * (L / 2) + L + len_down_stream <= s.size() ) ? ( Target_w[i] * (L / 2) + L + len_down_stream ) : (s.size());
            size_t node_list_len  = window_down - window_up + 3;
            
            node_list = new NODE[node_list_len];
            node_i = 0;
            add_node(window_up+10, 10, 0);
            uint32_t tmp, tar;
            for (size_t i=window_up; i < window_down; ++i)
            {
                tmp = (tmp << 2) | get_c[s[i]];
                tar = tmp & MASK;
                if (i - window_up >= 10)
                {
                    uint32_t d;
                    d = search(tar,index_rs);
                    if (d)
                    {
                        
                        while (( read_wlist[d-1]>>32) == tar)
                        {
                            if (check_node(i))
                            {
                                ++d;
                                continue;
                            }
                            add_node(i,read_wlist[d-1], 11);
                            ++d;
                            
                        }
                    }
                }
            }
            add_node(window_down, read_len, 0);
            
            
            create_matrix();
            getpath();
            
            
    
            string ss;
            t = do_alignment(const_cast<char*>(c_dna_f), window_up, window_down, const_cast<char*>(c_read), read_len, ss, out);
            if ((t) > score)
            {
                score = t;
                thess = ss;
            }
      

            
        }
        
        
        
        fprintf(out, "%.0lf\n", score);
     
        fprintf(out, "%s\n", thess.c_str());
 
    }
    printf("Time used = %.2f s\n",  (double)clock() / CLOCKS_PER_SEC);
    fclose(out);



    return 0;
}



 