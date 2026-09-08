#include <array>
#include <vector>
#include <queue>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <string>
#include <cstdint>
using namespace std;
constexpr int N=23, R=5;
using Perm=array<uint8_t,N>;
struct State { array<Perm,R> g; bool operator==(State const&o) const {return g==o.g;} bool operator<(State const&o) const {return g<o.g;} };
struct Hash { size_t operator()(State const&s) const noexcept { uint64_t h=1469598103934665603ULL; for(auto const&p:s.g) for(uint8_t x:p){h^=x;h*=1099511628211ULL;} return (size_t)h; } };
Perm idp(){Perm p{}; for(int i=0;i<N;i++) p[i]=i; return p;}
Perm mul(Perm const&p, Perm const&q){Perm r{}; for(int i=0;i<N;i++) r[i]=p[q[i]]; return r;}
Perm invp(Perm const&p){Perm r{}; for(int i=0;i<N;i++) r[p[i]]=i; return r;}
Perm cyc(initializer_list<vector<int>> cs){Perm p=idp(); for(auto c:cs){ vector<int> v; for(int x:c)v.push_back(x-1); for(size_t i=0;i<v.size();i++) p[v[i]]=v[(i+1)%v.size()]; } return p; }
State canonical(State const&t){
    array<Perm,2*R> gens; for(int i=0;i<R;i++){gens[i]=t.g[i]; gens[R+i]=invp(t.g[i]);}
    bool first=true; State best{};
    for(int root=0;root<N;root++){
        array<int,N> lab; lab.fill(-1); lab[root]=0; int nxt=1; array<int,N> q{}; int qh=0,qt=0; q[qt++]=root;
        while(qh<qt){int x=q[qh++]; for(auto const&gg:gens){int y=gg[x]; if(lab[y]<0){lab[y]=nxt++; q[qt++]=y;}}}
        if(nxt!=N){cerr<<"nontransitive\n"; exit(2);} State rt{};
        for(int k=0;k<R;k++) for(int old=0;old<N;old++) rt.g[k][lab[old]]=lab[t.g[k][old]];
        if(first||rt<best){best=rt; first=false;}
    }
    return best;
}
State braid(State const&t,int i){State L=t; auto gi=L.g[i], gj=L.g[i+1]; L.g[i]=mul(mul(gi,gj),invp(gi)); L.g[i+1]=gi; return L;}
State action(State const&s,int gn){ if(gn<3) return canonical(braid(s,gn)); State t=braid(s,3); t=braid(t,3); return canonical(t); }
State iota_raw(State const&s){ State t{}; t.g[0]=invp(s.g[3]); t.g[1]=invp(s.g[2]); t.g[2]=invp(s.g[1]); t.g[3]=invp(s.g[0]); t.g[4]=invp(s.g[4]); return t; }
vector<int> cycle_widths(vector<int> const&p, map<int,int>&hist){int n=p.size(); vector<int>w(n); vector<char> seen(n); for(int i=0;i<n;i++) if(!seen[i]){vector<int> c; int j=i; while(!seen[j]){seen[j]=1;c.push_back(j);j=p[j];} hist[(int)c.size()]++; for(int v:c) w[v]=c.size();} return w;}
vector<int> invperm(vector<int> const&p){vector<int>r(p.size()); for(int i=0;i<(int)p.size();i++)r[p[i]]=i; return r;}
vector<int> comp(vector<int> const&p,vector<int> const&q){vector<int>r(p.size());for(int i=0;i<(int)p.size();i++)r[i]=p[q[i]];return r;}
int main(){
 State seed{};
 seed.g[0]=cyc({{2,6},{4,12},{5,8},{10,17},{13,22},{14,19},{15,16},{20,21}});
 seed.g[1]=cyc({{1,15},{2,19},{3,8},{4,17},{6,7},{11,23},{12,21},{14,18}});
 seed.g[2]=cyc({{2,7},{4,13},{6,8},{9,22},{10,11},{14,17},{15,20},{16,18}});
 seed.g[3]=cyc({{1,20},{3,7},{4,18},{6,19},{9,22},{10,23},{13,17},{14,16}});
 seed.g[4]=cyc({{3,19,14},{4,21,15},{5,7,8},{10,17,11},{12,16,20},{13,22,18}});
 seed=canonical(seed);
 unordered_map<State,int,Hash> idx; idx.reserve(30000); vector<State> states; states.reserve(22000); queue<int> qq; idx.emplace(seed,0); states.push_back(seed); qq.push(0);
 while(!qq.empty()){int i=qq.front();qq.pop(); for(int gn=0;gn<4;gn++){State ns=action(states[i],gn); auto [it,ins]=idx.emplace(ns,(int)states.size()); if(ins){states.push_back(ns);qq.push(it->second);}}}
 cout<<"ORBIT_SIZE "<<states.size()<<"\n";
 vector<vector<int>> ops(4,vector<int>(states.size()));
 for(int gn=0;gn<4;gn++) for(int i=0;i<(int)states.size();i++) ops[gn][i]=idx.at(action(states[i],gn));
 vector<vector<int>> widths(4); vector<map<int,int>> hists(4);
 for(int gn=0;gn<4;gn++){widths[gn]=cycle_widths(ops[gn],hists[gn]); cout<<"HIST q"<<(gn<3?to_string(gn+1):string("4sq")); for(auto [k,v]:hists[gn])cout<<" "<<k<<":"<<v; cout<<"\n";}
 vector<int> chromfix; for(int i=0;i<(int)states.size();i++) if(ops[3][i]==i)chromfix.push_back(i); cout<<"CHROM_FIXED "<<chromfix.size()<<"\n";
 vector<int> iota(states.size(),-1); int miss=0; for(int i=0;i<(int)states.size();i++){State t=canonical(iota_raw(states[i])); auto it=idx.find(t); if(it==idx.end())miss++; else iota[i]=it->second;} cout<<"IOTA_MISSING "<<miss<<"\n";
 int invok=1,iftot=0; if(!miss){for(int i=0;i<(int)states.size();i++){if(iota[iota[i]]!=i)invok=0;if(iota[i]==i)iftot++;}} cout<<"IOTA_INVOLUTION "<<invok<<" IOTA_FIXED_TOTAL "<<iftot<<"\n";
 for(int mono=0;mono<3;mono++){
   map<int,int> dist; vector<int> diag; for(int i:chromfix){dist[widths[mono][i]]++; if(widths[mono][i]==5)diag.push_back(i);} int w5cy=hists[mono][5];
   cout<<"MONO q"<<mono+1<<" FIXED_WIDTH_DIST"; for(auto[k,v]:dist)cout<<" "<<k<<":"<<v; cout<<" W5_CYCLES "<<w5cy<<" W5_STATES "<<5*w5cy<<" DIAG "<<diag.size()<<"\n";
   set<int>D(diag.begin(),diag.end()); bool pres=true; for(int i:diag) if(!D.count(iota[i]))pres=false; cout<<"DIAG_IOTA_PRESERVES "<<pres<<"\n";
   set<int>used; int dfix=0; vector<pair<int,int>> orbs; for(int i:diag){if(iota[i]==i)dfix++; if(!used.count(i)){int j=iota[i]; used.insert(i);used.insert(j);orbs.push_back({i,j});}}
   cout<<"DIAG_IOTA_FIXED "<<dfix<<" ORBITS"; for(auto [a,b]:orbs)cout<<" ("<<a<<","<<b<<")"; cout<<"\n";
 }
 // conjugation of named operators
 if(!miss){vector<pair<string,vector<int>>> cand; for(int gn=0;gn<4;gn++){string n=gn<3?"q"+to_string(gn+1):"q4sq";cand.push_back({n,ops[gn]});cand.push_back({n+"^-1",invperm(ops[gn])});}
 for(int gn=0;gn<4;gn++){auto con=comp(iota,comp(ops[gn],iota)); cout<<"IOTA_CONJ "<<(gn<3?"q"+to_string(gn+1):"q4sq")<<" ->"; for(auto const&x:cand)if(x.second==con)cout<<" "<<x.first; cout<<"\n";}}
 // Six monochrome pair-collision operators, using minimal S4 braid representatives.
 auto word=[&](initializer_list<int> seq){vector<int> r(states.size());std::iota(r.begin(),r.end(),0);for(int k:seq)r=comp(ops[k],r);return r;};
 auto conj=[&](vector<int> const&w,vector<int> const&p){return comp(w,comp(p,invperm(w)));};
 vector<vector<int>> monos;
 monos.push_back(ops[0]); // 12
 monos.push_back(ops[1]); // 23
 monos.push_back(ops[2]); // 34
 auto w2=word({1}); monos.push_back(conj(w2,ops[0])); // 13
 auto w3=word({2}); monos.push_back(conj(w3,ops[1])); // 24
 auto w32=word({1,2}); monos.push_back(conj(w32,ops[0])); // 14
 // Four mixed collision operators: move the selected 2A slot to position 4.
 vector<vector<int>> chroms;
 chroms.push_back(ops[3]);
 auto a3=word({2}); chroms.push_back(conj(invperm(a3),ops[3]));
 auto a2=word({1,2}); chroms.push_back(conj(invperm(a2),ops[3]));
 auto a1=word({0,1,2}); chroms.push_back(conj(invperm(a1),ops[3]));
 // deduplicate defensively
 auto dedup=[](vector<vector<int>>&L){vector<vector<int>>u;for(auto const&x:L){bool seen=false;for(auto const&y:u)if(x==y){seen=true;break;}if(!seen)u.push_back(x);}L.swap(u);};
 dedup(monos);dedup(chroms);
 cout<<"MONO_PAIR_OPERATORS "<<monos.size()<<" CHROM_PAIR_OPERATORS "<<chroms.size()<<"\n";
 vector<vector<int>> mw; for(auto const&p:monos){map<int,int>hh;mw.push_back(cycle_widths(p,hh));}
 vector<vector<char>> cf; for(auto const&p:chroms){vector<char>f(states.size());for(int i=0;i<(int)states.size();i++)f[i]=(p[i]==i);cf.push_back(move(f));}
 set<int> wholeD; map<int,int> cellhist,multhist; map<int,int> mult;
 for(int a=0;a<(int)monos.size();a++)for(int b=0;b<(int)chroms.size();b++){int cnt=0;for(int i=0;i<(int)states.size();i++)if(mw[a][i]==5&&cf[b][i]){cnt++;wholeD.insert(i);mult[i]++;}cellhist[cnt]++;}
 for(auto const&kv:mult)multhist[kv.second]++;
 int wholefix=0,wholepres=1;for(int i:wholeD){if(iota[i]==i)wholefix++;if(!wholeD.count(iota[i]))wholepres=0;}
 cout<<"WHOLE_DIAGONAL_SIZE "<<wholeD.size()<<" IOTA_PRESERVED "<<wholepres<<" IOTA_FIXED "<<wholefix<<" CELL_COUNT_HIST";for(auto[k,v]:cellhist)cout<<" "<<k<<":"<<v;cout<<" MULTIPLICITY_HIST";for(auto[k,v]:multhist)cout<<" "<<k<<":"<<v;cout<<"\n";
 int if_mono5=0,if_chrom=0,if_both=0;for(int i=0;i<(int)states.size();i++)if(iota[i]==i){bool a=false,b=false;for(auto const&w:mw)a|=(w[i]==5);for(auto const&f:cf)b|=f[i];if(a)if_mono5++;if(b)if_chrom++;if(a&&b)if_both++;}
 cout<<"IOTA_FIXED_DIAGNOSTIC MONO5 "<<if_mono5<<" CHROM_FIXED "<<if_chrom<<" BOTH "<<if_both<<"\n";

}
