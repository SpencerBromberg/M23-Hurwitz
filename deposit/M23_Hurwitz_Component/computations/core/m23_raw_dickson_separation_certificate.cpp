#include <array>
#include <vector>
#include <queue>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <cstdint>
using namespace std;
constexpr int N=23,R=5;
using Perm=array<uint8_t,N>;
struct State{array<Perm,R> g; bool operator==(State const&o)const{return g==o.g;} bool operator<(State const&o)const{return g<o.g;}};
struct Hash{size_t operator()(State const&s)const noexcept{uint64_t h=1469598103934665603ULL;for(auto const&p:s.g)for(uint8_t x:p){h^=x;h*=1099511628211ULL;}return(size_t)h;}};
Perm idp(){Perm p{};for(int i=0;i<N;i++)p[i]=i;return p;} Perm mul(Perm const&p,Perm const&q){Perm r{};for(int i=0;i<N;i++)r[i]=p[q[i]];return r;} Perm invp(Perm const&p){Perm r{};for(int i=0;i<N;i++)r[p[i]]=i;return r;}
Perm cyc(initializer_list<vector<int>> cs){Perm p=idp();for(auto c:cs){vector<int>v;for(int x:c)v.push_back(x-1);for(size_t i=0;i<v.size();i++)p[v[i]]=v[(i+1)%v.size()];}return p;}
State canonical(State const&t){array<Perm,2*R> gens;for(int i=0;i<R;i++){gens[i]=t.g[i];gens[R+i]=invp(t.g[i]);}bool first=true;State best{};for(int root=0;root<N;root++){array<int,N>lab;lab.fill(-1);lab[root]=0;int nxt=1;array<int,N>q{};int qh=0,qt=0;q[qt++]=root;while(qh<qt){int x=q[qh++];for(auto const&gg:gens){int y=gg[x];if(lab[y]<0){lab[y]=nxt++;q[qt++]=y;}}}if(nxt!=N){cerr<<"nontransitive\n";exit(2);}State rt{};for(int k=0;k<R;k++)for(int old=0;old<N;old++)rt.g[k][lab[old]]=lab[t.g[k][old]];if(first||rt<best){best=rt;first=false;}}return best;}
State braid(State const&t,int i){State L=t;auto gi=L.g[i],gj=L.g[i+1];L.g[i]=mul(mul(gi,gj),invp(gi));L.g[i+1]=gi;return L;}
vector<int> invperm(vector<int> const&p){vector<int>r(p.size());for(int i=0;i<(int)p.size();i++)r[p[i]]=i;return r;} vector<int> comp(vector<int> const&p,vector<int> const&q){vector<int>r(p.size());for(int i=0;i<(int)p.size();i++)r[i]=p[q[i]];return r;}
vector<int> cyctype(vector<int> const&p){vector<char>s(p.size());vector<int>L;for(int i=0;i<(int)p.size();i++)if(!s[i]){int n=0,v=i;while(!s[v]){s[v]=1;n++;v=p[v];}L.push_back(n);}sort(L.begin(),L.end());return L;}
string ctstr(vector<int> const&p){map<int,int>h;for(int x:cyctype(p))h[x]++;string z;for(auto[k,v]:h){if(!z.empty())z+=" ";z+=to_string(k)+":"+to_string(v);}return z;}
int main(){
 State seed{};
 seed.g[0]=cyc({{2,6},{4,12},{5,8},{10,17},{13,22},{14,19},{15,16},{20,21}});
 seed.g[1]=cyc({{1,15},{2,19},{3,8},{4,17},{6,7},{11,23},{12,21},{14,18}});
 seed.g[2]=cyc({{2,7},{4,13},{6,8},{9,22},{10,11},{14,17},{15,20},{16,18}});
 seed.g[3]=cyc({{1,20},{3,7},{4,18},{6,19},{9,22},{10,23},{13,17},{14,16}});
 seed.g[4]=cyc({{3,19,14},{4,21,15},{5,7,8},{10,17,11},{12,16,20},{13,22,18}});
 seed=canonical(seed);
 // Colored orbit under Q0,Q1,Q2,Q3^2, retained to recover the exact D10 state keys.
 unordered_map<State,int,Hash> cidx;cidx.reserve(30000);vector<State> cs;queue<int>cq;cidx.emplace(seed,0);cs.push_back(seed);cq.push(0);
 auto cact=[&](State const&s,int gn){if(gn<3)return canonical(braid(s,gn));State t=braid(s,3);return canonical(braid(t,3));};
 while(!cq.empty()){int i=cq.front();cq.pop();for(int gn=0;gn<4;gn++){State ns=cact(cs[i],gn);auto[it,ins]=cidx.emplace(ns,(int)cs.size());if(ins){cs.push_back(ns);cq.push(it->second);}}}
 cout<<"COLORED_SIZE "<<cs.size()<<"\n";
 vector<vector<int>> cop(4,vector<int>(cs.size()));for(int gn=0;gn<4;gn++)for(int i=0;i<(int)cs.size();i++)cop[gn][i]=cidx.at(cact(cs[i],gn));
 // D10 = Fix(Q3^2) intersect width-5 cycles of Q0.
 vector<int>w(cop[0].size());vector<char>seen(cop[0].size());for(int i=0;i<(int)cop[0].size();i++)if(!seen[i]){vector<int>C;int v=i;while(!seen[v]){seen[v]=1;C.push_back(v);v=cop[0][v];}for(int u:C)w[u]=C.size();}
 vector<int>D10c;for(int i=0;i<(int)cs.size();i++)if(cop[3][i]==i&&w[i]==5)D10c.push_back(i);sort(D10c.begin(),D10c.end());cout<<"COLORED_D10";for(int v:D10c)cout<<" "<<v;cout<<"\n";
 // Full B5 orbit under four half twists Q0..Q3.
 unordered_map<State,int,Hash> idx;idx.reserve(150000);vector<State> st;st.reserve(110000);queue<int>q;idx.emplace(seed,0);st.push_back(seed);q.push(0);
 auto fact=[&](State const&s,int gn){return canonical(braid(s,gn));};
 while(!q.empty()){int i=q.front();q.pop();for(int gn=0;gn<4;gn++){State ns=fact(st[i],gn);auto[it,ins]=idx.emplace(ns,(int)st.size());if(ins){st.push_back(ns);q.push(it->second);}}}
 cout<<"FULL_SIZE "<<st.size()<<"\n";
 vector<vector<int>> Q(4,vector<int>(st.size()));for(int gn=0;gn<4;gn++)for(int i=0;i<(int)st.size();i++)Q[gn][i]=idx.at(fact(st[i],gn));
 vector<int>D10;for(int ci:D10c)D10.push_back(idx.at(cs[ci]));sort(D10.begin(),D10.end());cout<<"FULL_D10";for(int v:D10)cout<<" "<<v;cout<<"\n";
 vector<int>id(st.size());std::iota(id.begin(),id.end(),0);
 // Pure Dickson collision words A_{ij} in branch-slot frame.
 auto appperm=[&](vector<pair<int,int>> const&word){vector<int> p=id; for(auto [g,sg]:word){ if(sg==1)p=comp(Q[g],p); else p=comp(invperm(Q[g]),p);} return p;};
 auto Aword=[&](int i,int j){vector<pair<int,int>>w;for(int k=j-1;k>i;k--)w.push_back({k,1});w.push_back({i,1});w.push_back({i,1});for(int k=i+1;k<j;k++)w.push_back({k,-1});return w;};
 auto Wword=[&](int ss){vector<pair<int,int>>w;set<pair<int,int>>seenp;for(int i=0;i<5;i++){int j=(ss-i+5)%5;if(i==j)continue;auto pr=minmax(i,j);if(seenp.insert(pr).second){auto z=Aword(pr.first,pr.second);w.insert(w.end(),z.begin(),z.end());}}return w;};
 vector<vector<int>> W; for(int ss=0;ss<5;ss++){auto ww=Wword(ss);W.push_back(appperm(ww));cout<<"PURE_W "<<ss<<" WORDLEN "<<ww.size()<<" CTYPE "<<ctstr(W.back())<<"\n";}
 vector<vector<int>>gg=W;for(auto const&p:W)gg.push_back(invperm(p));vector<char>sv(st.size());map<int,int>hist;vector<vector<int>>ten;
 for(int ss=0;ss<(int)st.size();ss++)if(!sv[ss]){queue<int>b;vector<int>O;sv[ss]=1;b.push(ss);while(!b.empty()){int v=b.front();b.pop();O.push_back(v);for(auto const&p:gg){int u=p[v];if(!sv[u]){sv[u]=1;b.push(u);}}}hist[O.size()]++;if(O.size()==10){sort(O.begin(),O.end());ten.push_back(O);}}
 cout<<"PURE_H_ORBIT_HIST";for(auto[k,v]:hist)cout<<" "<<k<<":"<<v;cout<<" TEN_COUNT "<<ten.size()<<"\n";
 set<int>Dset(D10.begin(),D10.end());for(size_t k=0;k<ten.size()&&k<30;k++){cout<<"PURE_TEN_ORBIT "<<k<<" MEMBERS";for(int v:ten[k])cout<<" "<<v;int hit=0;for(int v:ten[k])hit+=Dset.count(v);cout<<" D10_HIT "<<hit<<"\n";}
 auto prod=id;for(int k=0;k<5;k++)prod=comp(W[k],prod);cout<<"PURE_PRODUCT_FORWARD_ID "<<(prod==id)<<"\n";prod=id;for(int k=4;k>=0;k--)prod=comp(W[k],prod);cout<<"PURE_PRODUCT_REVERSE_ID "<<(prod==id)<<"\n";
 return 0;
}
