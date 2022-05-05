#include <iostream>
#include <stdio.h>
using namespace std;
#include <iomanip>
#include <time.h>
#include <pthread.h>
#include <semaphore.h>
#include <windows.h>
#include <pmmintrin.h>
#include <xmmintrin.h>
#include <immintrin.h>
typedef struct {
	int t_id; // 线程 id
}threadParam_t;

#define N 1555
#define cnt 1
#define NUM_THREAD 7
#define RANDOM_ADD 2

sem_t sem_main;
sem_t sem_workerstart[NUM_THREAD];
sem_t sem_workerend[NUM_THREAD];
__m128 t1, t2, t3, t4;

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
void* threadFunc_Sta(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号
	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//循环划分任务
		for(int i = k + 1 + t_id;i < N;i += NUM_THREAD){
			for(int j = k + 1;j < N;++j){
				M[i][j] = M[i][j] - M[i][k] * M[k][j];
			}
			M[i][k] = 0;
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
	return nullptr;
}

void* threadFunc_Sta_SIMD(void* param) {
	threadParam_t* p = (threadParam_t*)param;
	int t_id = p->t_id; //线程编号
	for(int k = 0;k < N;k++){
		//阻塞，等待主线程完成除法操作
		sem_wait(&sem_workerstart[t_id]);
		//循环划分任务
		for(int i = k + 1 + t_id;i < N;i += NUM_THREAD){
			float tmp[4] = { M[i][k], M[i][k], M[i][k], M[i][k]};
            t1 = _mm_loadu_ps(tmp);
			for (int j = N - 4; j >k;j -= 4){
                t2 = _mm_loadu_ps(M[i] + j);
                t3 = _mm_loadu_ps(M[k] + j);
                t4 = _mm_sub_ps(t2,_mm_mul_ps(t1, t3));
                _mm_storeu_ps(M[i] + j, t4);
            }
            for (int j = k + 1; j % 4 !=(N % 4); j++){
                M[i][j] = M[i][j] - M[i][k] * M[k][j];
            }
			M[i][k] = 0;
		}
		sem_post(&sem_main);
		sem_wait(&sem_workerend[t_id]);
	}
	pthread_exit(NULL);
	return nullptr;
}

int main() {
	cout<<NUM_THREAD<<" threads in total"<<endl;
	init();
	long long head, tail, freq; 
	QueryPerformanceFrequency((LARGE_INTEGER *)&freq);

	long time_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		//get start time
		QueryPerformanceCounter((LARGE_INTEGER *)&head);

		ori();
		//get end time
		QueryPerformanceCounter((LARGE_INTEGER *)&tail );

		time_record += (tail - head) * 1000.0 / freq;
	}
	float ori = (float)time_record;
	cout << "Ori =  " << ori << "ms" << endl;
	
	pthread_t thread_handles[NUM_THREAD];
	threadParam_t param[NUM_THREAD];
	time_record = 0;
	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		QueryPerformanceCounter((LARGE_INTEGER *)&head);

		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta, (void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){
			for(int j = k+1;j<N;j++)
				M[k][j] = M[k][j]/M[k][k];
			M[k][k] = 1;
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
		QueryPerformanceCounter((LARGE_INTEGER *)&tail );
		time_record += (tail - head) * 1000.0 / freq;
	}

	float sta = (float)time_record;
	cout << "Sta =  " << sta << "ms" << endl;
	cout<<"Ratio of sta = "<<ori/sta<<endl;


	for (int i = 0; i < cnt; ++i){
		m_reset();
		sem_init(&sem_main,0,0);
		for(int i = 0 ; i < NUM_THREAD ; ++i){
			sem_init(&sem_workerstart[i],0,0);
			sem_init(&sem_workerend[i],0,0);
		}
		//get start time
		QueryPerformanceCounter((LARGE_INTEGER *)&head);
		for(int t_id=0;t_id<NUM_THREAD;++t_id){
			param[t_id].t_id = t_id;
			pthread_create(&thread_handles[t_id], NULL, threadFunc_Sta_SIMD, (void*)&param[t_id]);
		}
		for(int k = 0 ; k < N ;++k){
			float tmp[4] = { M[k][k], M[k][k], M[k][k], M[k][k] };
        	t1 = _mm_loadu_ps(tmp);
			int j;
			for (j = k + 1; j <= N - 4; j += 4) {
            	t2 = _mm_loadu_ps(M[k] + j);
            	t3 = _mm_div_ps(t2, t1);
            	_mm_storeu_ps(M[k] + j, t3);
        	}
        	if (k % 4 != (N % 4)){
            	for (int j = k; j % 4 != ( N% 4); j++){
                	M[k][j] = M[k][j] / tmp[0];
            	}
       		}
        	for (int j = (N % 4) - 1; j>= 0; j--){
            	M[k][j] = M[k][j] / tmp[0];
        	}
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
		QueryPerformanceCounter((LARGE_INTEGER *)&tail );
		time_record += (tail - head) * 1000.0 / freq;
	}
	float sta_SIMID = (float)time_record;
	cout << "Sta_SIMD =  " << sta_SIMID << "ms" << endl;
	cout<<"Ratio of sta_SIMD = "<<ori/sta_SIMID<<endl;


	sem_destroy(&sem_main);
	for(int i = 0 ; i < NUM_THREAD ; ++i){
		sem_destroy(&sem_workerstart[i]);
		sem_destroy(&sem_workerend[i]);
	}
	
	return 0;
}
