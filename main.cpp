#include <bits/stdc++.h>
#include <unordered_map>

// #include "Graph.hpp"

#include "LinearHeap.h"

#define N 500100
#define ui unsigned int
#define ull unsigned long long
#define puu std::pair<ull,ull>
#define fi first
#define se second

#define DEBUG

using namespace __gnu_cxx;
typedef std::unordered_map<ull,ui> ITIM;

inline int read(){
    int x=0;char c=getchar();bool y=1;
    for(;c<'0' || c>'9';c=getchar()) if(c=='-') y=0;
    for(;c>='0' && c<='9';c=getchar()) x=(x<<1)+(x<<3)+c-48;
    return y?x:-x;
}

int cnt;

int debug_setting = 20;

struct MyBBMatrix{
    std::vector<ui> ver[N];
    std::vector<ui> point_array,point_id;
    int n,m,K;
    std::vector<ui> best_solution_points;
    ui best_solution_size; 
    int degree[N];
    int degree_in_S[N];

    std::map<int,int> mapper;
    std::map<int,int> mapper_inverse;

    int prio_id[N];

    int root_id;

    ITIM edge_map;

    ListLinearHeap *heap;

    MyBBMatrix(){
    }

    inline ull mp(ull x,ull y){return x * (n + 2) + y;}
    
    inline void initialize(){
        point_array.clear();
        point_id.clear();
        for(int i=0;i<n;++i)
        point_array.push_back(i),
        point_id.push_back(i),
        prio_id[i]=i;
        return;
    }
    inline double getval(int x,int i){
        return 1.*degree[x] * (n-i)/n + 1.*degree_in_S[x]*i/2;
    }
    inline void initialize_degree(){
        for(int i=0;i<n;++i)
        degree[i] = ver[i].size();
        for(int i=0;i<n;++i)
        degree_in_S[i] = 0;
        return;
    }
    inline void greedy_calc_kplex(){
        for(int i=0;i<n;++i) prio_id[i]=i;
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
            for(int j=0;j<=i;++j){
                // std::cerr << j << " ??? " << i << " " << degree_in_S[prio_id[j]] << " " << i-K << std::endl;
                if(degree_in_S[prio_id[j]] <= i - K) {
                    --best_solution_size;
                    goto S;
                }
            }
        }
        S:;
        for(int i=0;i<best_solution_size;++i) 
        best_solution_points.push_back(prio_id[i]);
        std::cerr << best_solution_size << " \n " << check_S() << std::endl;
        return;
    }
    inline void add(int x,int y){
        ver[x].push_back(y);
        ver[y].push_back(x);

        edge_map[mp(x,y)] = 1;
        edge_map[mp(y,x)] = 1;
    }

    struct swap_group{
        int fi,se,th,fo;
        swap_group(const int X=0,const int Y=0,const int Z=0,const int W=0):fi(X),se(Y),th(Z),fo(W){}
    };

    std::vector<swap_group> swap_vec;

    inline bool add_point_to_S(int &S_end,int &R_end,int i,int search_depth){
        if(S_end == R_end) return false;
        if(i<S_end) return false;
        if(i>=R_end) add_point_to_R(R_end,S_end,i,search_depth),i=R_end-1;
        std::swap(point_array[i],point_array[S_end]);
        std::swap(point_id[point_array[i]],point_id[point_array[S_end]]);
        for(auto i:ver[point_array[S_end]]){
            ++degree_in_S[i];
            heap->increment(i,1);
        }
        ++S_end;
        swap_vec.push_back(swap_group(i,S_end-1,0,search_depth));
        return true;
    }
    inline bool remove_point_from_S(int &S_end,int &R_end,int i,int search_depth){
        if(S_end == 0) return false;
        if(i>=S_end) return true;
        std::swap(point_array[i],point_array[S_end-1]);
        std::swap(point_id[point_array[i]],point_id[point_array[S_end-1]]);
        for(auto i:ver[point_array[S_end-1]]){
            --degree_in_S[i];
            heap->decrement(i,1);
        }
        swap_vec.push_back(swap_group(i,S_end-1,1,search_depth));
        --S_end;
        return true;
    }
    inline bool add_point_to_R(int &R_end,int &S_end,int i,int search_depth){
        if(i<R_end || i>n) return false;
        std::swap(point_array[i],point_array[R_end]);
        std::swap(point_id[point_array[i]],point_id[point_array[R_end]]);
        swap_vec.push_back(swap_group(i,R_end,2,search_depth));
        ++R_end;
        return true;
    }
    inline bool remove_point_from_R(int &R_end,int &S_end,int i,int search_depth){
        if(R_end==0) return false;
        if(i < S_end){
            remove_point_from_S(S_end,R_end,i,search_depth);
            std::swap(point_array[S_end],point_array[R_end-1]);
            std::swap(point_id[point_array[S_end]],point_id[point_array[R_end-1]]);
            swap_vec.push_back(swap_group(S_end,R_end-1,3,search_depth));
        }
        else{
            std::swap(point_array[i],point_array[R_end-1]);
            std::swap(point_id[point_array[i]],point_id[point_array[R_end-1]]);
            swap_vec.push_back(swap_group(i,R_end-1,3,search_depth));
        }
        --R_end;
        // point_array.pop_back();
        return true;
    }
    inline void invalid_operation(int depth,int &S_end,int &R_end){
        while(swap_vec.back().fo == depth && !swap_vec.empty()){
            int fi = swap_vec.back().fi,se = swap_vec.back().se,th = swap_vec.back().th;
            std::swap(point_array[fi],point_array[se]);
            std::swap(point_id[point_array[fi]],point_id[point_array[se]]);
            switch(th){
                case 0:--S_end;
                for(auto i:ver[point_array[fi]]){
                    --degree_in_S[i];
                }
                break;
                case 1:++S_end;
                for(auto i:ver[point_array[fi]]){
                    ++degree_in_S[i];
                }
                break;
                case 2:--R_end;break;
                case 3:++R_end;break;
            }
            swap_vec.pop_back();
        }
        return;
    }
    void read_graph(){
        n = read(),m=read();K=read();
        for(ui i=1;i<=m;++i){
            add(read(),read());
        }
    }

    
    inline void graph_reduce(){
        fprintf(stderr,"best solution size with K = %d: %d\n",K,best_solution_size);
        fprintf(stderr,"before graph_reduce: %d %d\n",n,m);
        fflush(stderr);
        edge_map.clear();
        std::set<int> s;
        for(int i=0;i<n;++i) s.insert(i),reduce_v[i] = 0;
        int min_degree = n;
        for(int i=0;i<n;++i){
            degree[i] = ver[i].size();
            min_degree = std::min(min_degree,degree[i]);
            if(degree[i] + K <= best_solution_size)
            reduce_v[i] = 1,my_queue.push(i);
        }
        while(!my_queue.empty()){
            int x = my_queue.front();my_queue.pop();
            degree[x] = 0;
            for(auto y:ver[x]){
                if(degree[y] == 0) continue;
                --degree[y];
                if(degree[y] + K <= best_solution_size && !reduce_v[y])
                reduce_v[y] = 1,my_queue.push(y);
            }
        }

        int las_n = n;

        point_array.clear();

        mapper.clear();

        for(int i=0;i<n;++i)if(degree[i] > 0){
            // std::cerr << degree[i] << std::endl;
            assert(degree[i] + K > best_solution_size);
            mapper[i] = point_array.size();
            mapper_inverse[point_array.size()] = i;
            point_array.push_back(i);
        }
        n = point_array.size();
        m = 0;

        edge_map.clear();
        
        for(int i=0;i<n;++i){
            int x = point_array[i];
            std::vector<unsigned int> new_ver;new_ver.clear();
            for(auto y:ver[x]){
                if(mapper.find(y)!=mapper.end())
                new_ver.push_back(mapper[y]),++m,edge_map[mp(i,mapper[y])] = 1,edge_map[mp(mapper[y],i)] = 1;
            }
            ver[i] = new_ver;
            assert(x >= i);
        }

        for(int i=0,lim = best_solution_points.size();i<lim;++i)
        best_solution_points[i] = mapper[best_solution_points[i]];

        for(int i=0;i<n;++i) degree[i] = ver[i].size(),assert(degree[i] + K > best_solution_size);

        m/=2;

        fprintf(stderr,"after graph_reduce: %d %d\n",n,m);
        fflush(stderr);

        for(int i=n;i<las_n;++i)
        ver[i].clear(),degree[i] = degree_in_S[i] = 0;
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

            if(degree_check_i + K < best_solution_size)
            return false;
        }
        return true;
    }
    void BBMatrix(int S_end,int R_end,int depth){
        if(R_end <= best_solution_size) return;++cnt;
        // std::cerr << R_end << " " << best_solution_size << std::endl;
        if(!(cnt & 1048575)){
            std::cerr << root_id << " root? " << S_end << " " << R_end << " " << cnt << " " << best_solution_size << std::endl;
            for(int i=0;i<S_end;++i){
                std::cerr << point_array[i] << " " << degree_in_S[point_array[i]] << std::endl;
                assert(degree_in_S[point_array[i]] >= 0);
            }
            if(!check_S())
            std::cerr << "Wrong Answer!!!" << std::endl;
        }
        if(R_end <= best_solution_size) return;
        for(int i=0;i<S_end;++i)    if(degree_in_S[point_array[i]] + K < S_end) return;
        
        // std::cerr << root_id << " root?kk " << S_end << " " << R_end << " " << cnt << " " << best_solution_size << std::endl;
        if(S_end > best_solution_size){
            best_solution_size = S_end;
            best_solution_points.clear();
            for(int i=0;i<best_solution_size;++i){
                std::cerr << degree_in_S[point_array[i]] << std::endl;
                best_solution_points.push_back(point_array[i]);
            }
        }
        if(S_end == R_end)
        return;
        // std::vector<int> record_vector;
        // for(int i=0;i<R_end;++i)
        // record_vector.push_back(point_array[i]);

        int heap_id,heap_key;
        while(heap->get_min(heap_id,heap_key)){
            if(heap_key + K - 1 < S_end){
                remove_point_from_R(R_end,S_end,point_id[heap_id],depth);
                heap->pop_min(heap_id,heap_key);
            }
        }
        while(heap->get_max(heap_id,heap_key)){
            if(heap_key + K >= R_end && point_id[heap_id] >= S_end){
                add_point_to_S(S_end,R_end,point_id[heap_id],depth);
                heap->pop_max(heap_id,heap_key);
            }
        }
        
        /*for(int i=0;i<R_end;++i){
            if(degree_in_S[point_array[i]] + K - 1 < S_end){
            remove_point_from_R(R_end,S_end,i,depth);
            --i;continue;}
            if(degree_in_S[point_array[i]] + K >= R_end && i>=S_end){
            if(add_point_to_S(S_end,R_end,i,depth)) --i;continue;}
        }*/

        if(S_end != R_end && degree_in_S[point_array[S_end]] + K > S_end){
            add_point_to_S(S_end,R_end,S_end,depth);
            BBMatrix(S_end,R_end,depth+1);
            remove_point_from_S(S_end,R_end,S_end-1,depth);
        }

        if(S_end != R_end){
            remove_point_from_R(R_end,S_end,S_end,depth);
            BBMatrix(S_end,R_end,depth+1);
            add_point_to_R(R_end,S_end,R_end,depth);
        }
        invalid_operation(depth,S_end,R_end);
        // for(int i=0;i<R_end;++i){
        //     assert(record_vector[i] == point_array[i]);
        // }
        return;
    }

    int initialize_certain_point(int x,int type = 0){
        std::set<int> point_set;
        std::set<ull> edge_set;
        root_id = x;
        point_array.clear();
        for(int i=0;i<n;++i)
        degree[i] = 0;
        for(int i=0;i<n;++i)
        degree_in_S[i] = 0;

        point_set.insert(x);
        for(auto y:ver[x]){
            if(edge_set.find(mp(x,y)) == edge_set.end()){
                edge_set.insert(mp(x,y));
                edge_set.insert(mp(y,x));
                ++degree[x];
                ++degree[y];
            }
            point_set.insert(y);
            for(auto z:ver[y]){
                if(edge_set.find(mp(y,z)) == edge_set.end()){
                    edge_set.insert(mp(y,z));
                    edge_set.insert(mp(z,y));
                    ++degree[y];
                    ++degree[z];
                }
                point_set.insert(z);
            }
        }

        if(type == 1){
            return point_set.size();
        }

        for(auto x:point_set) point_array.push_back(x);

        // for(auto element : edge_map){
        //     auto x = element.first;
        //     if(point_set.find(x/(n+2)) != point_set.end() && point_set.find(x%(n+2)) != point_set.end())
        //     ++degree[x%(n+2)];
        // }
        return 0;
    }

    bool reduce_v[N];
    std::queue<int> my_queue;

    int id[N];

    int point_search[N];

    int point_dis2_neighbor_number[N];

    #ifdef DEBUG
    int debug_point_v[N];
    #endif

    

    void search(){

        double t1=0,t2=0;

        clock_t clk1,clk2;
        
        for(int i=0;i<n;++i) point_search[i] = i,point_dis2_neighbor_number[i] = initialize_certain_point(i,1);

            std::cerr << "search_array:" << point_array.size()<< std::endl;
        std::sort(point_search,point_search+n,[&](int x,int y){
            return point_dis2_neighbor_number[x] > point_dis2_neighbor_number[y];
        });
        for(int i=0;i<n;++i) id[point_search[i]] = i;

        heap = new ListLinearHeap(n,n-1);
        heap->init(n,n-1,point_search,degree_in_S);

        for(int i=n-1;~i;--i){
            clk1 = clock();
            initialize_certain_point(point_search[i]);
            t1 += clock()-clk1;
            clk2 = clock();

            if(!(i&511))
            std::cerr << "search_id:" << point_search[i] << " ? " << point_array.size()<< " " << degree[point_search[i]] << " " << cnt << std::endl;
            
            BBMatrix(0,point_array.size(),0);

            degree[point_search[i]] = 0;
            for(auto y:ver[point_search[i]]){
                --degree[y];
                /*if(i && degree[y] < degree[point_search[i-1]]){
                    int Y = point_search[id[y]],X = point_search[i-1];
                    std::swap(point_search[id[y]],point_search[i-1]);
                    std::swap(id[Y],id[X]);
                }*/
                for(int j = 0,lim = ver[y].size();j<lim;++j){
                    if(ver[y][j] == point_search[i]){
                        std::swap(ver[y][j],ver[y][lim-1]);
                        ver[y].pop_back();
                        break;
                    }
                }
            }
            ver[point_search[i]].clear();

            #ifdef DEBUG
            debug_point_v[point_search[i]] = 1;
            #endif
            t2 += clock()-clk2;
            
            if(!(i&511)){
                fprintf(stderr,"t1:%.5f,t2:%.5f\n",t1,t2);
                fflush(stderr);
            }
        }
        return;
    }

    void output(){
        printf("%u\n",best_solution_size);
        for(auto i:best_solution_points){
            printf("%u ",mapper_inverse[i]);
        }
        std::cout << std::endl;
        return;
    }
};

int main(){

    clock_t time_begin = clock();

    freopen("data_4.in","r",stdin);
    freopen("data_4.out","w",stdout);
    MyBBMatrix myBBMatrix = MyBBMatrix();
    myBBMatrix.read_graph();
    
    
    myBBMatrix.initialize();
    myBBMatrix.initialize_degree();
    myBBMatrix.greedy_calc_kplex();
    myBBMatrix.initialize();
    myBBMatrix.initialize_degree();
    myBBMatrix.graph_reduce();


    // myBBMatrix.initialize();
    // myBBMatrix.initialize_degree();
    // myBBMatrix.greedy_calc_kplex();
    // myBBMatrix.initialize();
    // myBBMatrix.initialize_degree();
    // myBBMatrix.graph_reduce();
    
    myBBMatrix.search();
    myBBMatrix.output();

    fprintf(stderr,"time seconds:%.3fs",1.*(clock()-time_begin)/CLOCKS_PER_SEC);
    fflush(stderr);
    
    return 0;
}
