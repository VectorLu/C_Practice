/*
y = x^6 - 10x^5 - 26x^4 + 344x^3 + 193x^2 - 1846x - 1680
要求在 ( -8, +8 ) 间寻找使表达式达到最小值的 x，误差为 0.001
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SUM 20        // 定义群体的染色体数量
#define MAX_LOOP 1200 // 最大循环次数
#define ERROR 0.01 // 若两次最优值之差小于此数则认为结果没有改变
#define CROSS_P 0.7 // 交叉概率，所选中的双亲按此概率进行交叉
#define MP 0.04    // 变异概率

typedef struct gen // 定义染色体结构
{
	int info; // 染色体结构，用以整型数的后 14 位作为染色体编码
	float suitability; // 此染色体所对应的适应度函数值，即表达式的值
} T_GEN;

T_GEN gen_group[SUM]; // 定义有 20 个染色体的种群
T_GEN gen_new[SUM]; // 定义含有 20 个染色体的种群，记录交叉产生的自带染色体
T_GEN gen_result; // 记录上一轮循环中得到的最优染色体
// 当相邻两轮得到的最优值对应的适应度之差小于 ERROR 时，其值增 1，反之清零
int result_unchange_time; 			
                									
typedef struct log // 形成链表，记录每次循环所产生的最优的适应度
{
	float suitability;
	struct log *next;
} T_LOG;
T_LOG list_log;
T_LOG *head;
T_LOG *end;
int log_num; // 链表长度

/********* 下面是函数声明 **************/
void initiate(); // 初始化函数，主要负责产生初始种群
void evaluation(int flag); // 评估种群中各染色体的适应度，并据此进行排序
void cross();              // 交叉函数
void selection(); // 选择函数
int record(); // 记录每次循环所产生的最优解并判断循环是否终止
void mutation();       // 变异函数
void show_result(int); // 显示结果
/********* 以上函数由主函数直接调用 *****/
int rand_sign(float p); // 按概率 p 产生随机数 0, 1，其值为 1 的概率为 p
int rand_between(int i, int j); // 随机产生一个在 i, j 两个数之间的整数
int rand_gen(); // 随机产生一个由 14 个基因组成的染色体
int create_mask(int a); // 用于交叉操作
int d2b(float x);   // 对现实解空间的可能解 x 进行二进制编码（染色体形式）
float b2d(int x);   // 将二进制编码 x 转化为现实解空间的值

int main()
{
    int i, flag;
    flag = 0;
    initiate();     // 产生初始化种群
    evaluation(0);  // 对初始化种群进行评估、排序

    for (i = 0; i < MAX_LOOP; i++)
    { 
        // 进入进化循环，当循环次数超过 MAX_LOOP 所规定的值时终止循环
        // flag 的值保持为 0
        cross();    // 进行交叉操作
        evaluation(1);   // 对子种群进行评估、排序
        selection();    // 从父、子种群选择最优的 SUM 个作为新的父种群
        // 如果满足终止规则 1 时将 flag 置 1，停止循环
        if (record() == 1)  
        {
            flag = 1;
            break;
        }
        mutation();     // 进行变异操作
    }
    show_result(flag);  // 按 flag 值显示寻优结果
    return 0;
}


/*
1. 初始化随机序列
2. 建立规模为 SUM 的初始种群
3. 建立 T_LOG 链表头（以便更好地观察收敛情况，不是遗传算法的必要步骤）
*/
void initiate()
{
    int i, s_time;
    long l_time;
    l_time = time(NULL);
    s_time = (unsigned) l_time/2;
    srand(s_time);
    for (i = 0; i < SUM; i++)
    {
        // 建立初始种群
        gen_group[i].info = rand_gen();
    }
    gen_result.suitability = 1000;
    result_unchange_time = 0;

    //初始化链表
    head = end = (T_LOG *)malloc(sizeof(list_log));

    if (head == NULL)
    {
        printf("\n内存不够!\n");
        exit(0);
    }
    end->next = NULL;
    log_num = 1;
}


void evaluation(int flag)
{
    // 按 flag 的指示分别对父种群和子种群进行评估和排序
    int i, j;
    T_GEN *gen_p;

    // 临时变量用于交换
    int gen_info;
    float gen_suitability;
    float x;

    // 如果 flag 为 0，对父种群进行操作
    if (flag == 0)
    {
        gen_p = gen_group;
    }
    // 否则对子种群进行操作
    else
    {
        gen_p = gen_new;
    }

    // 计算各染色体对应的适应度
    for (i = 0; i < SUM; i++)
    {
        x = b2d(gen_p[i].info);
        gen_p[i].suitability = 
        x * ( x * ( x * ( x * ( x * (x-10)-26 ) + 344 ) + 193 ) - 1846) - 1680;
    }

    for (i = 0; i < SUM-1; i++)
    {
        // 按照适应度的大小进行排序
        for (j = i+1; j < SUM; j++)
        {
            if (gen_p[i].suitability > gen_p[j].suitability)
            {
                gen_info = gen_p[i].info;
                gen_p[i].info = gen_p[j].info;
                gen_p[j].info = gen_info;
                gen_suitability = gen_p[i].suitability;
                gen_p[i].suitability = gen_p[j].suitability;
                gen_p[j].suitability = gen_suitability;
            }
        }
    }
}

/*
1. 随机选择将要进行交叉的染色体对，
保证种群的每一个染色体都仅有唯一一次交叉的机会；
2. 按 CROSS_P 的概率对选择的染色体对进行交叉操作；
3. 未进行交叉的染色体直接复制作为子染色体。

在第 2 部分首先由 rand_between() 选择交叉位，
create_mask() 返回一个交叉位以右皆为 1 的二进制数，赋给 mask1;
mask2 为交叉位以左皆为 1 的二进制数，它们分别与父染色体 gen_group[i].info、
gen_group[j].info 进行与操作，所得的结果之和即为交叉的结果
*/
void cross()
{
    // 对父种群按概率 CROSS_P 进行交叉操作
    int i, j, k;
    k = 0;
    int mask1, mask2;
    int a[SUM];
    for ( i = 0; i < SUM; i++)
    {
        a[i] = 0;
    }
    
    for (i = 0; i < SUM; i++)
    {
        if (a[i] == 0)
        {
            for( ; ; )
            {
                // 随机找到一组未进行过交叉的染色体与 a[i] 交叉
                j = rand_between(i+1, SUM-1);
                if (a[j] == 0)
                {
                    break;  
                }
            }
            if (rand_sign(CROSS_P) == 1)
            {
                mask1 = create_mask(rand_between(0, 14));
                mask2 = -mask1;

                gen_new[k].info = 
                (gen_group[i].info)&mask1 + (gen_group[j].info)&mask2;

                gen_new[k+1].info = 
                (gen_group[i].info)&mask2 + (gen_group[j].info)&mask1;

                k += 2;
            }
            else
            {
                gen_new[k].info = gen_group[i].info;
                gen_new[k+1].info = gen_group[j].info;
                k += 2;
            }
            a[i] = a[j] = 1;
        }
    }
}


/*
selection() 函数从父代种群 gen_group 和子代种群 gen_new 中选择
最优的染色体组成新的父代种群。由于两个种群都已进行了排序，因此只需找到子代种群中
第 t 个个体使其适应度函数值满足下列条件：
大于父代种群中的第 NUM-t 个个体，并小于父代种群中的第 NUM-t+1 个个体，
然后用子代种群中的前 t 个个体代替父代种群中的后 t 个个体即可。
为了加快对 t 的寻找速度，从两个种群的中部开始比较。
*/
void selection()
{
    int i, j, k;
    j = 0;
    i = SUM/2 - 1;
    if (gen_group[i].suitability < gen_new[i].suitability)
    {
        for (j = 1; j < SUM/2; j++)
        {
            if (gen_group[i+j].suitability > gen_new[i-j].suitability)
            {break;}
        }
    }
    else
    {
        if (gen_group[i].suitability > gen_new[i].suitability)
        {
            for ( j = -1; j < -SUM/2; j--)
            {
                if (gen_group[i+j].suitability <= gen_new[i-j].suitability)
                {
                    break;
                }
            }
        }
    }
    
    for ( k = j; k < SUM/2 + 1; k++)
    {
        gen_group[i+k].info = gen_new[i-k].info;
        gen_group[i+k].suitability = gen_new[i-k].suitability;
    }
}

/*
record() 函数主要完成两个任务：
1. 用链表记录各轮循环中的最优值
2. 判断是否满足循环停止条件
*/
int record()
{
    float x;
    T_LOG *r;
    r = (T_LOG *) malloc(sizeof(list_log));
    if (r == NULL)
    {
        printf("\n内存不够！\n");
        exit(0);
    }
    r->next = NULL;
    // 记录各轮循环中的最优值
    end->suitability = gen_group[0].suitability;
    // TODO: 下面这行可能有印刷错误
    end->next = r;
    log_num++;

    // 下面部分判断是否满足循环停止条件 1
    x = gen_result.suitability - gen_group[0].suitability;
    if (x < 0) {x = -x;}
    if (x < ERROR)
    {
        result_unchange_time++;
        if (result_unchange_time >= 20)
        {
            return 1;
        }
    }
    else
    {
        gen_result.info = gen_group[0].info;
        gen_result.suitability = gen_group[0].suitability;
        result_unchange_time = 0;
    }
    return 0;
}