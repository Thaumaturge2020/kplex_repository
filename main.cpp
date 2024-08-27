#include <bits/stdc++.h>
#include <unordered_map>

// #include "Graph.hpp"

#define N 2000100
#define ui unsigned int
#define ull unsigned long long
#define puu std::pair<ull,ull>
#define fi first
#define se second

using namespace __gnu_cxx;
typedef std::unordered_map<ull,ui> ITIM;

inline int read(){
    int x=0;char c=getchar();bool y=1;
    for(;c<'0' || c>'9';c=getchar()) if(c=='-') y=0;
    for(;c>='0' && c<='9';c=getchar()) x=(x<<1)+(x<<3)+c-48;
    return y?x:-x;
}

int cnt;

struct MyBBMatrix{
    std::vector<ui> ver[N];
    std::vector<ui> point_array;
    ui n,m,K;
    std::vector<ui> best_solution_points;
    ui best_solution_size; 
    int degree[N];
    int degree_in_S[N];

    int prio_id[N];

    ITIM edge_map;

    inline ull mp(ull x,ull y){return x * (n + 2) + y;}
    
    inline void initialize(){
        for(int i=0;i<n;++i) point_array.push_back(i),prio_id[i]=i;
        return;
    }
    inline double getval(int x,int i){
        return 1.*degree[x] * (n-i)/n + 1.*degree_in_S[x]*i/n;
    }
    inline void initialize_degree(){
        for(int i=0;i<n;++i)
        degree[i] = ver[i].size();
        for(int i=0;i<n;++i)
        degree_in_S[i] = 0;
        return;
    }
    bool cmp(int x,int y){
        return degree[x] > degree[y];
    }
    inline void greedy_calc_kplex(){
        for(int i=0;i<n;++i){
            for(int j=i;j<n;++j){
                if(degree_in_S[prio_id[j]] > i - K && getval(prio_id[j],i) > getval(prio_id[i],i)){
                    std::swap(prio_id[i],prio_id[j]);
                }
            }
            best_solution_size = i+1;
            for(auto j:ver[prio_id[i]]){
                ++degree_in_S[j];
            }
            for(int j=0;j<i;++j){
                if(degree_in_S[prio_id[j]] <= i - K) break;
            }
        }
        for(int i=0;i<best_solution_size;++i) best_solution_points.push_back(prio_id[i]);
        return;
    }
    inline void add(int x,int y){
        ver[x].push_back(y);
        ver[y].push_back(x);

        edge_map[mp(x,y)] = 1;
        edge_map[mp(y,x)] = 1;
    }
    inline bool add_point_to_S(int &S_end,int &R_end,int i){
        if(S_end == R_end) return false;
        if(i<S_end) return false;
        std::swap(point_array[i],point_array[S_end]);
        for(auto i:ver[point_array[S_end]]){
            ++degree_in_S[i];
        }
        ++S_end;
        return true;
    }
    inline bool remove_point_from_S(int &S_end,int &R_end,int i){
        if(S_end == 0) return false;
        if(i>=S_end) return true;
        std::swap(point_array[i],point_array[S_end]);
        --S_end;
        for(auto i:ver[point_array[S_end]]){
            --degree_in_S[i];
        }
        return true;
    }
    inline bool add_point_to_R(int &R_end,int &S_end,int x){
        point_array.push_back(x);
        ++R_end;
        return true;
    }
    inline bool remove_point_from_R(int &R_end,int &S_end,int i){
        if(R_end==0) return false;
        if(i < S_end){
            remove_point_from_S(S_end,R_end,i);
            std::swap(point_array[S_end],point_array[R_end-1]);
        }
        else{
            std::swap(point_array[i],point_array[R_end-1]);
        }
        --R_end;
        point_array.pop_back();
        return true;
    }
    void read_graph(){
        n = read(),m=read();K=read();
        for(ui i=1;i<=m;++i){
            add(read(),read());
        }
    }
    void recalc_S(int S_end,int R_end){
        for(int i=0;i<R_end;++i)
        degree_in_S[point_array[i]] = 0;
        for(int i=0;i<S_end;++i){
            for(auto j:ver[i])
            ++degree_in_S[point_array[j]];
        }
        return;
    }
    bool check_S(){

        if(best_solution_size != best_solution_points.size()) return false;
        for(int i=0;i<best_solution_size ; ++i)
        std :: cerr << best_solution_points[i] << " ";
        std :: cerr << std :: endl;
        for(auto i : best_solution_points){
            int degree_check_i = 0;
            for(auto j : best_solution_points)
            if(edge_map.find(mp(i,j))!=edge_map.end())
            ++degree_check_i;

            std::cerr << i << " " << degree_check_i << std::endl;

            if(degree_check_i + K + 1 < best_solution_size)
            return false;
        }
        return true;
    }
    void BBMatrix(int S_end,int R_end){
        if(R_end <= best_solution_size)
        return;
    
        ++cnt;
        if(!(cnt & 65535)){
            std::cerr << S_end << " " << R_end << " " << cnt << " " << best_solution_size << std::endl;
            for(int i=0;i<S_end;++i){
                std::cerr << point_array[i] << " " << degree_in_S[point_array[i]] << std::endl;
            }
            if(!check_S())
            std::cerr << "Wrong Answer!!!" << std::endl;
        }
        if(R_end < best_solution_size) return;
        
        if(S_end != R_end)
        for(int i=0;i<R_end;++i){
            if(degree_in_S[point_array[i]] + K + 1 < R_end) goto S;
        }

        best_solution_size = R_end;
        best_solution_points = point_array;
        return;
        S:;

        for(int i=0;i<R_end;++i){
            if(degree_in_S[point_array[i]] + K + 1 < S_end){
            // std::cerr << degree_in_S[point_array[i]] + K << " " << S_end - 1 << std::endl;

            remove_point_from_R(R_end,S_end,i);--i;continue;            }
            if(degree_in_S[point_array[i]] + K + 1 >= R_end){
            if(add_point_to_S(S_end,R_end,i)) --i;continue;}
        }

        // std::cerr << " ? ?? " << " " << S_end << " " << R_end << std::endl;

        if(S_end != R_end && degree_in_S[point_array[S_end]] + K >= S_end){
            std::vector<ui> vec_ui = point_array;
            add_point_to_S(S_end,R_end,S_end);
            BBMatrix(S_end,R_end);
            point_array = vec_ui;
            recalc_S(S_end,R_end);
        }
        if(S_end != R_end){
            std::vector<ui> vec_ui = point_array;
            
            remove_point_from_R(R_end,S_end,S_end);

            BBMatrix(S_end,R_end);
            point_array = vec_ui;
            recalc_S(S_end,R_end);
        }         
        return;
    }

    void search(){
        BBMatrix(0,n);
        return;
    }

    void output(){
        printf("%u\n",best_solution_size);
        for(auto i:best_solution_points){
            printf("%u ",i);
        }
        std::cout << std::endl;
        return;
    }
};

int main(){
    freopen("data.in","r",stdin);
    freopen("data.out","w",stdout);
    MyBBMatrix myBBMatrix;
    myBBMatrix.read_graph();
    myBBMatrix.initialize_degree();
    myBBMatrix.greedy_calc_kplex();
    myBBMatrix.initialize_degree();
    myBBMatrix.search();
    myBBMatrix.output();
    return 0;
}
