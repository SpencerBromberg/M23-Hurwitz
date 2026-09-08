#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;

static constexpr int N=23;
using Perm=array<uint8_t,N>;
template<size_t K> using Tup=array<Perm,K>;

Perm idperm(){ Perm p{}; for(int i=0;i<N;i++) p[i]=i; return p; }
Perm mul(const Perm&p,const Perm&q){ Perm r{}; for(int i=0;i<N;i++) r[i]=p[q[i]]; return r; } // p after q
Perm invp(const Perm&p){ Perm r{}; for(int i=0;i<N;i++) r[p[i]]=i; return r; }
Perm cyc(initializer_list<vector<int>> cs){ Perm p=idperm(); for(auto c:cs){ for(auto &x:c)--x; for(size_t i=0;i<c.size();++i)p[c[i]]=c[(i+1)%c.size()]; } return p; }
int ord(const Perm&p){ vector<int> seen(N); int o=1; for(int i=0;i<N;i++) if(!seen[i]){ int x=i,l=0; while(!seen[x]){seen[x]=1;l++;x=p[x];} o=lcm(o,l);} return o; }
vector<int> lens(const Perm&p){ vector<int> seen(N),out; for(int i=0;i<N;i++) if(!seen[i]){int x=i,l=0;while(!seen[x]){seen[x]=1;l++;x=p[x];}out.push_back(l);} sort(out.begin(),out.end()); return out; }
int idx(const Perm&p){int s=0;for(int l:lens(p))s+=l-1;return s;}

template<size_t K> Perm product(const Tup<K>&t){Perm r=idperm();for(auto &g:t)r=mul(r,g);return r;}
template<size_t K> bool transitive(const Tup<K>&t){ vector<int>seen(N); queue<int>q; seen[0]=1;q.push(0); vector<Perm> gs; for(auto &g:t){gs.push_back(g);gs.push_back(invp(g));} while(!q.empty()){int x=q.front();q.pop();for(auto &g:gs){int y=g[x];if(!seen[y]){seen[y]=1;q.push(y);}}}return accumulate(seen.begin(),seen.end(),0)==N; }

template<size_t K> Tup<K> qmove(Tup<K> t,int i,int sign=1){ if(sign==1){Perm a=t[i],b=t[i+1];t[i]=mul(mul(a,b),invp(a));t[i+1]=a;} else {Perm a=t[i],b=t[i+1];t[i]=b;t[i+1]=mul(mul(invp(b),a),b);} return t; }

template<size_t K> struct THash{ size_t operator()(Tup<K> const&t)const noexcept{ uint64_t h=1469598103934665603ULL; for(auto &p:t)for(auto x:p){h^=(uint64_t)x+1;h*=1099511628211ULL;}return (size_t)h;} };

template<size_t K> Tup<K> relabel(const Tup<K>&t,const array<int,N>&lab){Tup<K>r{};for(int k=0;k<K;k++)for(int old=0;old<N;old++)r[k][lab[old]]=lab[t[k][old]];return r;}

template<size_t K> Tup<K> canonical(const Tup<K>&t){
 vector<Perm>gs;for(auto &g:t){gs.push_back(g);gs.push_back(invp(g));}
 bool have=false; Tup<K> best{};
 for(int root=0;root<N;root++){
   array<int,N> lab; lab.fill(-1); lab[root]=0; queue<int>q;q.push(root);int next=1;
   while(!q.empty()){int x=q.front();q.pop();for(auto &g:gs){int y=g[x];if(lab[y]<0){lab[y]=next++;q.push(y);}}}
   if(next!=N) throw runtime_error("nontransitive canonicalization");
   auto r=relabel(t,lab);
   if(!have || r<best){best=r;have=true;}
 }
 return best;
}

Tup<5> genact(Tup<5>t,int g){ if(g<3)return qmove(t,g,1); if(g==3){t=qmove(t,3,1);return qmove(t,3,1);} throw runtime_error("bad gen"); }

string histstr(const vector<int>&p){ vector<int>vis(p.size()); map<int,int>h;for(size_t i=0;i<p.size();i++)if(!vis[i]){int x=i,l=0;while(!vis[x]){vis[x]=1;l++;x=p[x];}h[l]++;} string s;for(auto [l,c]:h){if(!s.empty())s+=" ";s+=to_string(l)+"^"+to_string(c);}return s;}
vector<int> cycle_len_each(const vector<int>&p){ vector<int>res(p.size()),vis(p.size());for(size_t i=0;i<p.size();i++)if(!vis[i]){vector<int>c;int x=i;while(!vis[x]){vis[x]=1;c.push_back(x);x=p[x];}for(int y:c)res[y]=c.size();}return res;}

Tup<4> fuse_and_prepare(const Tup<5>&s,int i){
 Tup<4> t{}; int k=0; for(int j=0;j<5;j++){ if(j==i){t[k++]=mul(s[j],s[j+1]);j++;} else t[k++]=s[j]; }
 // Hurwitz bubble sort by type rank 2 < 3 < 5; preserves product one.
 auto rank=[](const Perm&p){int o=ord(p);return o==2?0:o==3?1:o==5?2:9;};
 for(int pass=0;pass<6;pass++) for(int j=0;j<3;j++) if(rank(t[j])>rank(t[j+1])) t=qmove(t,j,1);
 return t;
}

int main(){
 Tup<5> seed={
 cyc({{2,6},{4,12},{5,8},{10,17},{13,22},{14,19},{15,16},{20,21}}),
 cyc({{1,15},{2,19},{3,8},{4,17},{6,7},{11,23},{12,21},{14,18}}),
 cyc({{2,7},{4,13},{6,8},{9,22},{10,11},{14,17},{15,20},{16,18}}),
 cyc({{1,20},{3,7},{4,18},{6,19},{9,22},{10,23},{13,17},{14,16}}),
 cyc({{3,19,14},{4,21,15},{5,7,8},{10,17,11},{12,16,20},{13,22,18}})};
 if(product(seed)!=idperm()||!transitive(seed)) throw runtime_error("seed failure");

 unordered_map<Tup<5>,int,THash<5>> pos; vector<Tup<5>> states; queue<int> qq;
 auto start=canonical(seed);pos[start]=0;states.push_back(start);qq.push(0);
 while(!qq.empty()){int ii=qq.front();qq.pop();for(int g=0;g<4;g++){auto n=canonical(genact(states[ii],g));auto [it,ins]=pos.emplace(n,(int)states.size());if(ins){states.push_back(n);qq.push(it->second);}}}
 cout<<"FIVE_BRANCH_COMPONENT_SIZE "<<states.size()<<"\n";
 if(states.size()!=21456) return 2;
 vector<vector<int>> ops(4,vector<int>(states.size()));
 for(size_t i=0;i<states.size();i++)for(int g=0;g<4;g++){auto n=canonical(genact(states[i],g));ops[g][i]=pos.at(n);}
 for(int g=0;g<4;g++)cout<<"Q"<<g+1<<(g==3?"SQ":"")<<"_HIST "<<histstr(ops[g])<<"\n";
 vector<string> expected={"1^24 2^292 3^796 4^1248 5^980 6^1428","1^24 2^292 3^796 4^1248 5^980 6^1428","1^24 2^292 3^796 4^1248 5^980 6^1428","1^18 2^468 3^582 4^984 5^900 7^840 8^192 11^264"};
 for(int g=0;g<4;g++) if(histstr(ops[g])!=expected[g]) return 3;

 // Width-five clutching: each 5-cycle produces one ordered (2A,2A,3A,5A) inner class.
 for(int g=0;g<3;g++){
   vector<int>vis(states.size()); set<Tup<4>> keys; int cycles5=0,bad=0;
   for(size_t i=0;i<states.size();i++) if(!vis[i]){vector<int>c;int x=i;while(!vis[x]){vis[x]=1;c.push_back(x);x=ops[g][x];}if(c.size()==5){cycles5++;auto f=fuse_and_prepare(states[c[0]],g);vector<int>o;for(auto&p:f)o.push_back(ord(p));if(o!=vector<int>({2,2,3,5})||product(f)!=idperm()||!transitive(f)){bad++;continue;}keys.insert(canonical(f));}}
   cout<<"Q"<<g+1<<"_WIDTH5_CYCLES "<<cycles5<<"\n";
   cout<<"Q"<<g+1<<"_CLUTCHED_DISTINCT_INNER_KEYS "<<keys.size()<<"\n";
   cout<<"Q"<<g+1<<"_CLUTCHED_INVALID "<<bad<<"\n";
   if(cycles5!=980||keys.size()!=980||bad) return 4;
 }

 // Mixed boundary: q4^2 fixed set and restricted q1,q2 monodromy.
 vector<int> fixed; vector<int> local(states.size(),-1);for(size_t i=0;i<states.size();i++)if(ops[3][i]==(int)i){local[i]=fixed.size();fixed.push_back(i);} cout<<"MIXED_FIXED_SIZE "<<fixed.size()<<"\n"; if(fixed.size()!=18)return 5;
 vector<int>a(18),b(18);for(int i=0;i<18;i++){int x=fixed[i];if(local[ops[0][x]]<0||local[ops[1][x]]<0)return 6;a[i]=local[ops[0][x]];b[i]=local[ops[1][x]];}
 cout<<"MIXED_Q1_HIST "<<histstr(a)<<"\n";cout<<"MIXED_Q2_HIST "<<histstr(b)<<"\n";
 vector<int>ab(18);for(int i=0;i<18;i++)ab[i]=b[a[i]];cout<<"MIXED_Q1Q2_HIST "<<histstr(ab)<<"\n";
 if(histstr(a)!="2^1 3^2 5^2"||histstr(b)!="2^1 3^2 5^2"||histstr(ab)!="3^6")return 7;
 // transitivity
 vector<int>seen(18);queue<int>q;q.push(0);seen[0]=1;while(!q.empty()){int x=q.front();q.pop();for(auto &p:{a,b}){int y=p[x];if(!seen[y]){seen[y]=1;q.push(y);}}}int conn=accumulate(seen.begin(),seen.end(),0);cout<<"MIXED_MONODROMY_ORBIT_SIZE "<<conn<<"\n";if(conn!=18)return 8;
 // centralizer by possible image of 0
 vector<vector<int>> cents;
 for(int target=0;target<18;target++){vector<int>phi(18,-1);queue<int>z;phi[0]=target;z.push(0);bool ok=true;while(!z.empty()&&ok){int x=z.front();z.pop();for(auto &p:{a,b}){int y=p[x],py=p[phi[x]];if(phi[y]<0){phi[y]=py;z.push(y);}else if(phi[y]!=py){ok=false;break;}}}if(ok){set<int>s(phi.begin(),phi.end());if(s.size()==18)cents.push_back(phi);}}
 cout<<"MIXED_CENTRALIZER_ORDER "<<cents.size()<<"\n"; if(cents.size()!=2)return 9; int free_nontr=0, nontr_idx=-1;for(int ci=0;ci<(int)cents.size();ci++){auto &c=cents[ci];bool id=true,free=true;for(int i=0;i<18;i++){if(c[i]!=i)id=false;if(c[i]==i)free=false;}if(!id&&free){free_nontr++;nontr_idx=ci;}}cout<<"MIXED_FREE_NONTRIVIAL_CENTRALIZER_ELEMENTS "<<free_nontr<<"\n";if(free_nontr!=1)return 10;
 cout<<"MIXED_Q1_LOCAL";for(int x:a)cout<<" "<<x;cout<<"\n";
 cout<<"MIXED_Q2_LOCAL";for(int x:b)cout<<" "<<x;cout<<"\n";
 cout<<"MIXED_CENTRAL_INVOLUTION";for(int x:cents[nontr_idx])cout<<" "<<x;cout<<"\n";
 auto q1len=cycle_len_each(ops[0]); set<int>R10;for(int x:fixed)if(q1len[x]==5)R10.insert(x);
 cout<<"R10_LOCAL";for(int x:R10)cout<<" "<<local[x];cout<<"\n";
 cout<<"R10_GLOBAL";for(int x:R10)cout<<" "<<x;cout<<"\n";
 set<int>q2R;for(int x:R10)q2R.insert(ops[1][x]);cout<<"R10_SIZE "<<R10.size()<<"\n";cout<<"Q2_R10_EQUALS_R10 "<<(q2R==R10?"True":"False")<<"\n";if(R10.size()!=10||q2R==R10)return 11;
 cout<<"ALL_FIVE_BRANCH_CLUTCHING_MIXED_CHECKS_PASS True\n";
}
