#include <iostream>
using namespace std;
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>
#include <arm_neon.h>
typedef struct {
	int t_id; // 线程 id
}threadParam_t;

#define N 996
#define cnt 1
#define NUM_THREAD 5
#define RANDOM_ADD 2

sem_t sem_leader;
sem_t sem_Division[NUM_THREAD-1];
sem_t sem_Elimination[NUM_THREAD-1];

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
			int temp = (k*2+i*3) % N;
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
			int temp = (k*2+i*3) % N;
			for (int j = 0; j < N; j++)
				M[temp][j] += M[i][j];
		}
	}
}
void output() {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cout << setw(12) << fixed << setprecision(2)<< M[i][j];
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
void* threadFunc_Sta2(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号
	for(int k = 0;k < N;k++){
		//t_id为0的线程做除法
		if(t_id == 0){
			for(int j = k+1;j<N;j++)
				M[k][j] = M[k][j]/M[k][k];
			M[k][k] = 1;
		}
		else
			sem_wait(&sem_Division[t_id-1]);
		//t_id为0的线程唤醒其他线程
		if(t_id == 0)
			for(int i = 0;i<NUM_THREAD - 1;++i)
				sem_post(&sem_Division[i]);
	
		//循环划分任务
		for(int i = k + 1 + t_id ; i<N;i+=NUM_THREAD){
			for(int j = k + 1;j<N;++j)
				M[i][j] = M[i][j] - M[i][k] * M[k][j];
			M[i][k] = 0;
		}
		if(t_id == 0){
			for(int i = 0;i<NUM_THREAD-1;i++)
				sem_wait(&sem_leader);
			for(int i = 0;i<NUM_THREAD-1;i++)
				sem_post(&sem_Elimination[i]);
		}
		else{
			sem_post(&sem_leader);
			sem_wait(&sem_Elimination[t_id-1]);
		}

	}
	pthread_exit(NULL);
}

void* threadFunc_Sta2_Multidiv(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号
	for(int k = 0;k < N;k++){
		//t_id为0的线程做除法
		if(t_id == 0){
			float32x4_t v0 = vmovq_n_f32(M[k][k]);
			float32x4_t v1;
			int j;
			for(j = k+1;j<=N-4;j+=4){
				v1 = vld1q_f32(M[k]+j);
				v0 = vdivq_f32(v1,v0);
				vst1q_f32(M[k] + j, v0);
			}
			for(j = j-4;j<N;j++){
				M[k][j] /=M[k][k];
			}
			M[k][k] = 1;
		}
		else
			sem_wait(&sem_Division[t_id-1]);
		//t_id为0的线程唤醒其他线程
		if(t_id == 0)
			for(int i = 0;i<NUM_THREAD - 1;++i)
				sem_post(&sem_Division[i]);
	
		//循环划分任务
		for(int i = k + 1 + t_id ; i<N;i+=NUM_THREAD){
			for(int j = k + 1;j<N;++j)
				M[i][j] = M[i][j] - M[i][k] * M[k][j];
			M[i][k] = 0;
		}
		if(t_id == 0){
			for(int i = 0;i<NUM_THREAD-1;i++)
				sem_wait(&sem_leader);
			for(int i = 0;i<NUM_THREAD-1;i++)
				sem_post(&sem_Elimination[i]);
		}
		else{
			sem_post(&sem_leader);
			sem_wait(&sem_Elimination[t_id-1]);
		}

	}
	pthread_exit(NULL);
}

int main() {
	cout<<NUM_THREAD<<" threads in total"<<endl;
	init();
	struct timespec sts,ets;
	time_t dsec;
	long dnsec;
	
	int dsec_record = 0;
	long dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		//get start time
		timespec_get(&sts, TIME_UTC);
		ori();
		//get end time
		timespec_get(&ets, TIME_UTC);
		dsec=ets.tv_sec-sts.tv_sec;
		dnsec=ets.tv_nsec-sts.tv_nsec;

		dsec_record+=dsec;
		dnsec_record+=dnsec;
	}
	if (dnsec_record<0){
		dnsec_record--;
		dnsec_record+=1000000000ll;
	}
	float ori = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Ori =  " << ori << "ms" << endl;
	
	pthread_t thread_handles[NUM_THREAD];
	threadParam_t param[NUM_THREAD];
	dsec_record = 0;
	dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_leader,0,0);
		for(int i = 0 ; i < NUM_THREAD - 1 ; ++i){
			sem_init(&sem_Division[i],0,0);
			sem_init(&sem_Elimination[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta2, (void*)&param[t_id]);
		}
		
		for(int t_id = 0;t_id<NUM_THREAD;++t_id)
			pthread_join(thread_handles[t_id], NULL);
		
		//get end time
		timespec_get(&ets, TIME_UTC);
		dsec=ets.tv_sec-sts.tv_sec;
		dnsec=ets.tv_nsec-sts.tv_nsec;
		dsec_record+=dsec;
		dnsec_record+=dnsec;
	}
	if (dnsec_record<0){
		dnsec_record--;
		dnsec_record+=1000000000ll;
	}
	float sta2 = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta2 =  " << sta2 << "ms" << endl;
	cout<<"Ratio of sta2 = "<<ori/sta2<<endl;


	dsec_record = 0;
	dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_leader,0,0);
		for(int i = 0 ; i < NUM_THREAD - 1 ; ++i){
			sem_init(&sem_Division[i],0,0);
			sem_init(&sem_Elimination[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta2_Multidiv, (void*)&param[t_id]);
		}
		
		for(int t_id = 0;t_id<NUM_THREAD;++t_id)
			pthread_join(thread_handles[t_id], NULL);
		
		//get end time
		timespec_get(&ets, TIME_UTC);
		dsec=ets.tv_sec-sts.tv_sec;
		dnsec=ets.tv_nsec-sts.tv_nsec;
		dsec_record+=dsec;
		dnsec_record+=dnsec;
	}
	if (dnsec_record<0){
		dnsec_record--;
		dnsec_record+=1000000000ll;
	}
	float sta2mul = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta2_Multidiv =  " << sta2mul << "ms" << endl;
	cout<<"Ratio of sta2_Multidiv = "<<ori/sta2mul<<endl;


	sem_destroy(&sem_leader);
	for(int i = 0 ; i < NUM_THREAD ; ++i){
		sem_destroy(&sem_Division[i]);
		sem_destroy(&sem_Elimination[i]);
	}
	
	return 0;
}
