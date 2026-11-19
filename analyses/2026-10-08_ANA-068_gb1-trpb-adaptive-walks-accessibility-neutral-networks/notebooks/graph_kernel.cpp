#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <vector>
using namespace std;
struct UF { vector<int> p,s; UF(int n):p(n),s(n,1){iota(p.begin(),p.end(),0);} int find(int a){while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a;} void join(int a,int b){a=find(a);b=find(b);if(a==b)return;if(s[a]<s[b])swap(a,b);p[b]=a;s[a]+=s[b];}};
// Row-major output columns: degree, beneficial degree, neutral degree,
// uniform expected length, greedy expected length, uniform terminal rank,
// greedy terminal rank, uniform hit-top, greedy hit-top, uniform hit-best,
// greedy hit-best, any strict path to top, any strict path to best,
// any nondecreasing path to top, any nondecreasing path to best, rank,
// uniform terminal fitness, greedy terminal fitness.
extern "C" int analyse(int n,int q,int dim,const double* f,const uint8_t* valid,
 double eps,double top,double floor,double* out,int32_t* comp,int32_t* activecomp){
 vector<int> stride(dim,1), order;for(int j=dim-2;j>=0;--j)stride[j]=stride[j+1]*q;
 for(int i=0;i<n;i++)if(valid[i])order.push_back(i);
 sort(order.begin(),order.end(),[&](int a,int b){return f[a]>f[b]||(f[a]==f[b]&&a<b);});
 if(order.empty())return 1;double best=f[order[0]];int nv=order.size();
 for(int k=0;k<n*18;k++)out[k]=0;
 for(int b=0;b<nv;){int z=b+1;while(z<nv&&f[order[z]]==f[order[b]])z++;
  double rank=1.-(b+z)/(2.*nv);for(int k=b;k<z;k++)out[order[k]*18+15]=rank;b=z;}
 UF all(n),active(n);
 for(int u:order){double* a=out+u*18;double vmax=-INFINITY;int ng=0;
  for(int j=0;j<dim;j++){int base=u-((u/stride[j])%q)*stride[j];
   for(int k=0;k<q;k++){int v=base+k*stride[j];if(v==u||!valid[v])continue;
    double df=f[v]-f[u];a[0]++;
    if(abs(df)<=eps){a[2]++;if(u<v){all.join(u,v);if(f[u]>floor&&f[v]>floor)active.join(u,v);}}
    if(df>eps){const double* b=out+v*18;a[1]++;a[3]+=1+b[3];a[5]+=b[5];a[7]+=b[7];a[9]+=b[9];a[16]+=b[16];a[11]=max(a[11],b[11]);a[12]=max(a[12],b[12]);
      if(f[v]>vmax){vmax=f[v];ng=0;a[4]=a[6]=a[8]=a[10]=a[17]=0;}
      if(f[v]==vmax){ng++;a[4]+=1+b[4];a[6]+=b[6];a[8]+=b[8];a[10]+=b[10];a[17]+=b[17];}
    }
   }
  }
  if(a[1]){for(int k:{3,5,7,9,16})a[k]/=a[1];for(int k:{4,6,8,10,17})a[k]/=ng;}
  else{a[5]=a[6]=a[15];a[16]=a[17]=f[u];}
  if(f[u]>=top){a[7]=a[8]=a[11]=1;}
  if(f[u]==best){a[9]=a[10]=a[12]=1;}
 }
 // Reverse BFS on nondecreasing edges, including exact plateaus; epsilon does not enter.
 for(int which=0;which<2;which++){vector<int> queue;int col=13+which;
  for(int u:order)if(which?f[u]==best:f[u]>=top){out[u*18+col]=1;queue.push_back(u);}
  for(size_t idx=0;idx<queue.size();idx++){int v=queue[idx];
   for(int j=0;j<dim;j++){int base=v-((v/stride[j])%q)*stride[j];for(int k=0;k<q;k++){int u=base+k*stride[j];
    if(valid[u]&&!out[u*18+col]&&f[u]<=f[v]){out[u*18+col]=1;queue.push_back(u);}
   }}
  }
 }
 for(int u=0;u<n;u++){comp[u]=valid[u]?all.find(u):-1;activecomp[u]=(valid[u]&&f[u]>floor)?active.find(u):-1;}
 return 0;
}
extern "C" int walks(int n,int q,int dim,const double* f,const uint8_t* valid,
 const int32_t* starts,int nw,int steps,double eps,int greedy,uint64_t seed,int32_t* paths){
 vector<int> stride(dim,1);for(int j=dim-2;j>=0;--j)stride[j]=stride[j+1]*q;
 mt19937_64 rng(seed);
 for(int w=0;w<nw;w++){int u=starts[w];paths[w*(steps+1)]=u;
  for(int t=1;t<=steps;t++){vector<int> choices;double best=-INFINITY;
   for(int j=0;j<dim;j++){int base=u-((u/stride[j])%q)*stride[j];for(int k=0;k<q;k++){int v=base+k*stride[j];
    if(v==u||!valid[v]||f[v]-f[u]<=eps)continue;
    if(!greedy)choices.push_back(v);else if(f[v]>=best){if(f[v]>best){choices.clear();best=f[v];}choices.push_back(v);}
   }}
   if(choices.empty()){for(int k=t;k<=steps;k++)paths[w*(steps+1)+k]=u;break;}
   if(!choices.empty()){uniform_int_distribution<int> dist(0,(int)choices.size()-1);u=choices[dist(rng)];}
   paths[w*(steps+1)+t]=u;
  }
 }return 0;
}
