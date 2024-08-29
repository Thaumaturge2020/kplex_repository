#include <bits/stdc++.h>

using namespace std;

int main(){
    freopen("data4.in","r",stdin);
    freopen("data4.out","w",stdout);
    int n,m;
    scanf("%d %d",&n,&m);
        printf("%d %d\n",n,m);
    for(int i=1;i<=m;++i){
        int u,v;
        scanf("%d %d",&u,&v);
        --u,--v;
        printf("%d %d\n",u,v);
    }
    return 0;
}