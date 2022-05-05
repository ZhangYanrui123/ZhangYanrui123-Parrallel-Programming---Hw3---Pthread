#include <iostream>
using namespace std;
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <arm_neon.h>
typedef struct {
	int k; //消去的轮次
	int t_id; // 线程 id
}threadParam_t;

#define N 4096
#define cnt 1
#define NUM_THREAD 4
#define RANDOM_ADD 10
#define OUTPUTN 6

float** M;
void init() {
	M = new float* [N];
	for (int i = 0; i < N; i++) {
		M[i] = new float[N];
		for (int j = 0; j < i; j++)
			M[i][j] = 0;

		for (int j = i; j < N; j++)
			M[i][j] = rand() % 50;
	}
	for (int k = 0; k < RANDOM_ADD; k++) {
		for (int i = 0; i < N; i++) {
			int temp = rand() % N;
			for (int j = 0; j < N; j++)
				M[temp][j] += M[i][j];
		}
	}
}
void m_reset() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < i; j++)
			M[i][j] = 0;

		for (int j = i; j < N; j++)
			M[i][j] = rand() % 50;
	}
	for (int k = 0; k < RANDOM_ADD; k++) {
		for (int i = 0; i < N; i++) {
			int temp = rand() % N;
			for (int j = 0; j < N; j++)
				M[temp][j] += M[i][j];
		}
	}
}
void output() {
	for (int i = 0; i < OUTPUTN; i++) {
		for (int j = 0; j < OUTPUTN; j++) {
			cout << setw(6) << fixed << setprecision(2)<< M[i][j];
		}
		cout << endl;
	}
}
void ori() {
	for (int k = 0; k < N; k++) {
		float tmp = M[k][k];
		for (int j = k; j < N; j++) {
			M[k][j] = M[k][j] / tmp;
		}
		for (int i = k + 1; i < N; i++) {
			float tmp2 = M[i][k];
			for (int j = k + 1; j < N; j++) {
				M[i][j] = M[i][j] - tmp2 * M[k][j];
			}
			M[i][k] = 0;
		}
	}
}
void* threadFunc_Act(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int k = p->k; //消去的轮次
	int t_id = p->t_id; //线程编号
	int i = k + t_id + 1; //获取自己的计算任务

	float32x4_t v0;
	float32x4_t v1;
	float32x4_t v2;
	for(;i<N;i+=NUM_THREAD){
		v1 = vmovq_n_f32(M[i][k]);
//        float32x4_t v2;
        int j;
		for (j = k + 1; j <= N - 4; j+= 4) {
			v2 = vld1q_f32(M[k] + j);
            v0 = vld1q_f32(M[i] + j);
            v2 = vmulq_f32(v1, v2);
            v0 = vsubq_f32(v0, v2);
            vst1q_f32(M[i] + j, v0);
        }
		for (j = j-4; j < N; j++)
            M[i][j] = M[i][j] - M[i][k] * M[k][j];
        M[i][k] = 0;		
	}
	pthread_exit(NULL);
	return nullptr;
}


int main() {
	cout<<NUM_THREAD<<" threads in total"<<endl;
	init();

	struct timespec sts,ets;
	time_t dsec;
	long dnsec;
	
	float dsec_record = 0;
	float dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		timespec_get(&sts, TIME_UTC);
		ori();
		timespec_get(&ets, TIME_UTC);
		dsec=ets.tv_sec-sts.tv_sec;
		dnsec=ets.tv_nsec-sts.tv_nsec;
		if (dnsec_record<0){
			dnsec_record--;
			dnsec_record+=1000000000ll;
		}
		dsec_record+=(float)dsec;
		dnsec_record+=(float)dnsec;
	}
	float ori = dsec_record * 1000 + dnsec_record/1000000;
	cout << "Ori =  " << ori << "ms" << endl;

	dnsec_record = 0;
	dsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		timespec_get(&sts, TIME_UTC);
		for (int k = 0; k < N; ++k) {
			float32x4_t v1 = vmovq_n_f32(M[k][k]);
        	float32x4_t v0;
        	int j;
			for (j = k + 1; j <= N - 4; j += 4) {
            	v0 = vld1q_f32(M[k] + i);
            	v0 = vdivq_f32(v0, v1);
            	vst1q_f32(M[k] + i, v0);
        	}
        	for (j = j-4; j < N; j++)
            	M[k][i] = M[k][i] / M[k][k];

			M[k][k] = 1.0;
			int work_count = NUM_THREAD;
			pthread_t* thread_handles;
			threadParam_t* param;
			//创建工作线程，进行消去操作
			if(N - k - 1>=NUM_THREAD){
				thread_handles = new pthread_t[NUM_THREAD];// 创建对应的 Handle
				param = new threadParam_t[NUM_THREAD];// 创建对应的线程数据结构
			}
			else{
				work_count = N - k - 1;
				thread_handles = new pthread_t[work_count];// 创建对应的 Handle
				param = new threadParam_t[work_count];// 创建对应的线程数据结构
			}
			//分配任务
			for (int t_id = 0; t_id < work_count; t_id++) {
				param[t_id].k = k;
				param[t_id].t_id = t_id;
			}
			//创建线程
			for (int t_id = 0; t_id < work_count; t_id++)
				pthread_create(&thread_handles[t_id], NULL, threadFunc_Act, (void*)&param[t_id]);

			//主线程挂起等待所有的工作线程完成此轮消去工作
			for (int t_id = 0; t_id < work_count; t_id++)
				pthread_join(thread_handles[t_id], NULL);

		}
		
		timespec_get(&ets, TIME_UTC);
		dsec=ets.tv_sec-sts.tv_sec;
		dnsec=ets.tv_nsec-sts.tv_nsec;
		if (dnsec_record<0){
			dnsec_record--;
			dnsec_record+=1000000000ll;
		}
		dsec_record+=(float)dsec;
		dnsec_record+=(float)dnsec;
	}

	float act = dsec_record * 1000 + dnsec_record/1000000;
	cout << "Act =  " << act << "ms" << endl;
	cout << "Ratio of Act = " << ori/act << endl;
	return 0;
}
