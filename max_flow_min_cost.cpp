/*
 created by czr September，2015  for Mathematical Modeling
给出某个城市的45个点的经纬度（可以计算出点与点之间的距离）和其能提供和需要的流。
计算出最小费用最大流，最大流一定是将正数点上的流全部流到负数的点上
最终所有点变成零
 */

#include<iostream>
#include<cstring>
#include<cstdlib>
#include<cstdio>
#include<climits>
#include<algorithm>
#include<queue>
#include<math.h>
using namespace std;


int win;
int n,m;
const int N=17;
const int M=51;
const int MAX=0xffffff;
int pre[M];//存储前驱顶点
double dist[M];//存储到源点s的距离
double ff;
int inq[M];//每个顶点是否在队列中的标志
double min_c_f;//记录增广路径中的残留容量
int vertex;//顶点数
double sum;//保存最小费用
double out[16][4];
double inin[45]={120.094969230769,30.3410076923077,-45.933333334,
120.30225,30.4161166666667,-63.933333334,
120.21485,30.20213,-70.933333334,
120.163075,30.323625,76.066666666,
120.08296,30.167,1.066666666,
120.30606,30.30648,-47.933333334,
120.160536363636,30.2005636363636,-8.933333334,
120.298016666667,30.2230833333333,-17.933333334,
120.196511111111,30.26835,220.066666666,
120.0395,30.2657375,-70.933333334,
120.370492307692,30.3129615384615,-114.933333334,
120.188533333333,30.1183333333333,2.066666666,
120.1169,30.278675,126.066666666,
120.269253846154,30.1662846153846,-3.933333334,
120.225555555556,30.3572333333333,20.066666666};


struct element
{
    double c;//容量
    double f;//流
    double c_f;//残留容量
    double v;//价值
} G[N][N];

struct man//记录源点的坐标
{
    double x,y;
    int num;
} man[N];
struct house//记录汇点坐标
{
    double x,y;
    int num;
} house[N];


double EARTH_RADIUS = 6378.137;//地球半径
double rad(double d)
{
   return d * 3.14159265 / 180.0;
}

double GetDistance(double lng1, double lng2, double lat1, double lat2)
{
   double radLat1 = rad(lat1);
   double radLat2 = rad(lat2);
   double a = radLat1 - radLat2;
   double b = rad(lng1) - rad(lng2);

   double s = 2 * asin(sqrt(pow(sin(a/2),2) +
    cos(radLat1)* cos(radLat2)*pow(sin(b/2),2)));
   s = s * EARTH_RADIUS;
   s = round(s * 10000) / 10000;
   return s;
}

double absmin(double a,double b){
    if (a <-b){ return a;}
    else {return  -b;}
}


int mcase,hcase;//记录有多少个源点汇点
void init()
{
    sum=0;
    n=15;
    for (int i=0;i<45;i++){
        //cout<<"out"<<(i)/3+1<<" "<<(i)%3 +1<<endl;
        out[(i)/3+1][((i)%3) +1]=inin[i];
    }
    mcase=0;
    hcase=0;
    for(int i=1; i<=n; i++)
    {
            //cin>>coord[i][j];
            if(out[i][3]>0.001)
            {
                mcase++;
                man[mcase].num=i;
                man[mcase].x=out[i][1];
                man[mcase].y=out[i][2];
            }
            if(out[i][3]<0.001)
            {
                hcase++;
                house[hcase].num=i;
                house[hcase].x=out[i][1];
                house[hcase].y=out[i][2];

            }
    }
    
    vertex=mcase+hcase+1;//加入超源点0和超汇点，注意要+1，即抽象成网络流的结构
    for(int u=0; u<=vertex; u++)
    {
        for(int v=0; v<=vertex; v++)
        {
            G[u][v].c=G[v][u].c=0;
            G[u][v].c_f=G[v][u].c_f=0;
            G[u][v].f=G[v][u].f=0;
            G[u][v].v=G[v][u].v=MAX;
        }
    }
    
    for(int i=1; i<=mcase; i++)
    {
        G[0][i].v=0;//从超源点到各个源点之间的权值取为0
        G[0][i].c=G[0][i].c_f=out[man[i].num][3];
        for(int j=1; j<=hcase; j++)
        {
            double w=GetDistance(house[j].x,man[i].x,house[j].y,man[i].y);
            G[i][mcase+j].v=w;
            G[i][mcase+j].c=absmin(out[man[i].num][3],out[house[j].num][3]);
            G[i][mcase+j].c_f=G[i][mcase+j].c;
            G[mcase+j][vertex].v=0;
            double ttt=-out[house[j].num][3];
            G[mcase+j][vertex].c=ttt;
            G[mcase+j][vertex].c_f= G[mcase+j][vertex].c;
        }
    }

}

void SPFA(int s)//求最短路径的SPFA算法
{
    //cout<<"in===";
    queue<int> Q;
    int u;
    for(int i=0; i<=vertex; i++)//初始化
    {
        dist[i]=MAX;
        pre[i]=-1;
        inq[i]=0;
    }
    dist[s]=0;
    Q.push(s);
    inq[s] = 1;
    while(!Q.empty())
    {
        u=Q.front();
        Q.pop();
        inq[u]=0;
        for(int i=0; i<=vertex; i++)//更新u的邻接点的dist[], pre[], inq[]
        {
            int v=i;
            if(G[u][v].c_f<=1e-1)     // 表示(u,v)没有边
                continue;
            if(G[u][v].v==MAX){
                G[u][v].v=-G[v][u].v;
            }
            if(dist[v]-0.01>dist[u]+G[u][v].v)//松弛操作
            {

                dist[v]=dist[u]+G[u][v].v;

                pre[v]=u;
                if(pre[3]==14&&pre[14]==3){
                    ;
                }
                if(inq[v]==0)
                {
                    Q.push(v);
                    inq[v]=1;
                }
            }
        }
    }
    //cout<<"===out";
}

void ford_fulkerson(int s,int t)
{
    SPFA(s);
    while(pre[t]!=-1)//pre为-1表示没有找到从s到t的增广路径
    {
        if(pre[3]==14&&pre[14]==3){
            ;
        }

        min_c_f=MAX;
        int u=pre[t], v=t;//计算增广路径上的残留容量
        while(u!=-1)
        {
            if(min_c_f > G[u][v].c_f)
                min_c_f=G[u][v].c_f;
            v=u;
            u=pre[v];
            if(pre[3]==14&&pre[14]==3){
                ;
            }
        }
        ff+=min_c_f;
        cout<<ff<<endl;
        if(ff>435)win=1;
        u=pre[t], v=t;
        while(u!=-1)
        {
            G[u][v].f+=min_c_f; //修改流
            G[v][u].f=-G[u][v].f;
            G[u][v].c_f=G[u][v].c-G[u][v].f; //修改残留容量
            G[v][u].c_f=G[u][v].f;
            v=u;
            u=pre[v];
        }
        SPFA(s);
    }
}

int main()
{
    ff=0;//流
    init();
    ford_fulkerson(0,vertex);//计算从超源点0到超汇点vertex之间的最小费用最大流
    double sum2; //费用
    for(int i=1; i<=mcase; i++)
    {
        for(int j=1; j<=hcase; j++)
        {
            //cout<<G[i][j+mcase].c<<" ";
            sum2+=fabs(G[i][mcase+j].f)*G[i][mcase+j].v;
            
        }
        cout<<endl;
    }
    cout<<ff<<" "<<sum2<<endl; //输出最大流最小费用
    return 0;
}
