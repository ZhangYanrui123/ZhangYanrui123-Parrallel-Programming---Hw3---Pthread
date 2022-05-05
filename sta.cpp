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

#define N 1024
#define cnt 5
#define NUM_THREAD 7
#define RANDOM_ADD 2

sem_t sem_main;
sem_t sem_workerstart[NUM_THREAD];
sem_t sem_workerend[NUM_THREAD];

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

void* threadFunc_Sta_NoSIMD(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号


	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//循环划分任务
		for(int i = k + 1 + t_id;i < N;i += NUM_THREAD){
			for(int j = k+1;j<N;j++)
				M[i][j] = M[i][j] - M[i][k] * M[k][j];
			
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
}

void* threadFunc_Sta(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号

	float32x4_t v0;
	float32x4_t v1;
	float32x4_t v2;

	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//循环划分任务
		for(int i = k + 1 + t_id;i < N;i += NUM_THREAD){
			v1 = vmovq_n_f32(M[i][k]);
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
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
}

void* threadFunc_Sta_Vertical_noSIMD(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号

	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//块划分任务
		for(int j = k + t_id + 1;j<N;j+=NUM_THREAD){
			M[k][j]/=M[k][k];
			for(int i = k+1;i<N;i++){
				M[i][j] = M[i][j] - M[i][k]*M[k][j]; 
			}
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
}

void* threadFunc_Sta_Vertical(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号

	float32x4_t v0;
	float32x4_t v1;
	float32x4_t v2;

	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//块划分任务
		int j;
		for(j = k+1;j<=N-4;j+=4){
			v1 = vld1q_f32(M[k+1]+j);
			for(int i = k+t_id+1;i<N;i++){
				v0 = vmovq_n_f32(M[i][k]);
				v2 = vld1q_f32(M[i]+j);
				v0 = vmulq_f32(v1,v0);
				v2 = vsubq_f32(v2,v0);
				vst1q_f32(M[i]+j,v2);
			}
		}
		for(j = j-4;j<N;j++){
			for(int i = k + t_id + 1;i<N;i+=NUM_THREAD){
				M[i][j] = M[i][j] - M[i][k] * M[k][j];
				M[i][k] = 0;
			}
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
}

void* threadFunc_Sta_Block(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号

	float32x4_t v0;
	float32x4_t v1;
	float32x4_t v2;

	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//块划分任务
		int my_task = (N-k-1)/NUM_THREAD;
		//如果任务太少不能块划分就串行
		
		int my_end = k + (t_id + 1) * my_task + 1;

		for(int i = k + t_id * my_task + 1;i < my_end;i++){
			v1 = vmovq_n_f32(M[i][k]);
			int j;
			for (j = k + 1; j <= N - 4; j+= 4) {
				v2 = vld1q_f32(M[k] + j);
	            v2 = vmulq_f32(v1, v2);
            	v0 = vld1q_f32(M[i] + j);
	            v0 = vsubq_f32(v0, v2);
            	vst1q_f32(M[i] + j, v0);
        	}
			for (j = j-4; j < N; j++)
            	M[i][j] = M[i][j] - M[i][k] * M[k][j];
			M[i][k] = 0;
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
}

int main() {
	cout<<"N = "<<N<<endl;
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
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta_NoSIMD, (void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){

			for(int j = k+1;j<N;j++)
				M[k][j]/=M[k][k];

			//唤醒工作线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerstart[t_id]);
				
			//主线程睡眠
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_wait(&sem_main);
			//主线程再次唤醒下一轮消去线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerend[t_id]);
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
	float staNoSIMD = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta_NoSIMD =  " << staNoSIMD << "ms" << endl;
	cout<<"Ratio of staNoSIMD = "<<ori/staNoSIMD<<endl;

	

	dsec_record = 0;
	dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta_Vertical_noSIMD ,(void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){


			//唤醒工作线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerstart[t_id]);
			for(int i = k+1;i<N;i++)
				M[i][k] = 0;
			//主线程睡眠
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_wait(&sem_main);
			//主线程再次唤醒下一轮消去线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerend[t_id]);
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
	float stavNoSIMD = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta_v_NoSIMD =  " << stavNoSIMD << "ms" << endl;
	cout<<"Ratio of Sta_v_NoSIMD = "<<ori/stavNoSIMD<<endl;


	dsec_record = 0;
	dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta, (void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){
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
			//唤醒工作线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerstart[t_id]);
				
			//主线程睡眠
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_wait(&sem_main);
			//主线程再次唤醒下一轮消去线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerend[t_id]);
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
	float sta = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta =  " << sta << "ms" << endl;
	cout<<"Ratio of sta = "<<ori/sta<<endl;


	dsec_record = 0;
	dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta_Block, (void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){
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

            for(int i = N-(N-k-1)%NUM_THREAD+1;i<N;i++)
            	for(int j = k+1;j<N;j++)
            		M[i][j] = M[i][j] - M[i][k] * M[j][k];
			//唤醒工作线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerstart[t_id]);
				
			//主线程睡眠
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_wait(&sem_main);
			//主线程再次唤醒下一轮消去线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerend[t_id]);
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
	float stablo = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta_Block =  " << stablo << "ms" << endl;
	cout<<"Ratio of stablo = "<<ori/stablo<<endl;


	dsec_record = 0;
	dnsec_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		timespec_get(&sts, TIME_UTC);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta_Vertical, (void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){
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
			//唤醒工作线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerstart[t_id]);
				
			//主线程睡眠
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_wait(&sem_main);
			//主线程再次唤醒下一轮消去线程
			for(int t_id = 0;t_id<NUM_THREAD;++t_id)
				sem_post(&sem_workerend[t_id]);
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
	float staver = (float)dsec_record * 1000 + (float)dnsec_record / 1000000;
	cout << "Sta_Vertical =  " << staver << "ms" << endl;
	cout<<"Ratio of staver = "<<ori/staver<<endl;


	sem_destroy(&sem_main);
	for(int i = 0 ; i < NUM_THREAD ; ++i){
		sem_destroy(&sem_workerstart[i]);
		sem_destroy(&sem_workerend[i]);
	}
	
	return 0;
}
