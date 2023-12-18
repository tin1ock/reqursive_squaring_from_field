#include "cuda_runtime.h"
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <cmath>

// if N >= 65 535 then error


#include <vector>
#include <stdio.h>



__device__ int division_f(int a, int b, int f){
        int result = 0;
        int count = 1;
        while (count <= f){
            if ((b*count)%f == a){
                return count;
            }
            count += 1;
        }
        return -1;
    }

 __device__ int getsize(int *ptr, int field){ //works fine
         //std::cout<<ptr[0]<<ptr[1]<<ptr[2]<<ptr[3]<<std::endl;
        int result = 0;

        while (ptr[result] != field+1 ) 
        {
        result++;
        }
        return result+1;
    }

__device__ int diff_f(int a, int b, int f){
        if ((a-b)<0){
            return a-b+f;
        }
        else{
            return a-b;
        }
    }

__device__ int sum_f(int a, int b, int f){
        if ((a+b)>=f){
            return (a+b)-f;
        }
        else{
            return a+b;
        }
    }




__device__ int division(int* A_in, int* B_in, int f, int* answer){
        const int a_len = getsize(A_in, f);
        const int b_len = getsize(B_in, f);
        int* B = new int[b_len];
        for (int i=0;i<b_len;i++) B[i] = B_in[i];

        int result_flag = 0;
        int* result = new int[a_len]; 
        for (int i=0; i<a_len; i++) {
            result[result_flag] = A_in[i]; 
            result_flag++;
            }
        int code = 0;
        
        while (true) {
            int cursorIndex = 1;
            if (getsize(result, f) < getsize(B, f)) {
                break;
            }
            else{

                int res_size = getsize(result,f);
                int* A = new int[res_size]; 
                int* tmp; 
                tmp=result; 
                for (int i = 0; i<res_size; i++)  
                {
                    A[i] = *tmp; 
                    tmp++;
                }  

                result_flag = 0;
                for (int r=0; r<res_size;r++){
                    result[r] = f+1;
                }
                int c = division_f(A[0], B[0], f);

                if (c==-1){
                    code = -1;
                    break;
                }
                int flag = 0;
                for (int i = 1; i<getsize(B, f)-1; i++){ 
                    cursorIndex += 1;
                    if (flag == 0) {
                        if (A[i] == c*B[i]){
                            continue;
                        }
                        else {
                            flag = 1;
                            result[result_flag] = diff_f(A[i]%f, c*B[i]%f, f); 
                            result_flag++;
                        }
                    }
                    else {
                        result[result_flag] = diff_f(A[i]%f, c*B[i]%f, f); 
                        result_flag++;
                    }
                    
                }
                if (cursorIndex < getsize(A, f)-1){
                    for (cursorIndex; cursorIndex<getsize(A, f); cursorIndex++){
                        result[result_flag] = A[cursorIndex]; 
                        result_flag++;
                    }
                }
            }

        }
        int i =0;

        if ((getsize(result,f) == 2) && result[0] == 0){
            
                answer[0] = f+1;
            
        }
        else{
                for (int an=0; an<getsize(result,f); an++){
                    answer[an] = result[an];
                }
            }
        return 0;
    }
    


int strInt(std::string inputString){
    int result = 0;
    int razr = 1;
    for (int i = inputString.length()-1; i>=0; i--){
        int f  = inputString[i]-'0';
        result = result + f*razr;
        razr = razr*10;
    }
    return result;
}

__device__ bool bigger(int* A, int* B, int f){

    if (getsize(A,f) > getsize(B,f)) return true;
    if (getsize(A,f) < getsize(B,f)) return false;
    for (int i=0; i<getsize(A,f); i++){
        if (A[i] > B[i]) return true;
        if (A[i] < B[i]) return false;
    }
}


class Polinom {     
    public:
    std::vector<int> koefs;
   
    Polinom(std::vector<int> input_koefs){
        for (int i: input_koefs)
            this->koefs.push_back(i);
    }
    void printPolinom(){
        for (int i: koefs){
            std::cout<<i;
        }
        std::cout<<std::endl;
    }

    std::vector<int> getPolinom(){
        return koefs;
    }

    int getLength(){
        return koefs.size();
    }

    int division_f(int a, int b, int f){
        int result = 0;
        int count = 1;
        while (count <= f){
            if ((b*count)%f == a){
                return count;
            }
            count += 1;
        }
        return -1;
    }

    int diff_f(int a, int b, int f){
        if ((a-b)<0){
            return a-b+f;
        }
        else{
            return a-b;
        }
    }

    int sum_f(int a, int b, int f){
        if ((a+b)>=f){
            return (a+b)-f;
        }
        else{
            return a+b;
        }
    }

    Polinom sum (Polinom B, int f){
        std::vector<int> result = koefs;
        if (koefs.size() > B.koefs.size()){
           int dif = koefs.size() - B.koefs.size();
           for (int i =0; i<B.koefs.size(); i++){
                result[dif+i] = result[dif+i] + B.koefs[i];
           }
        }
        else{
        result = B.koefs;
        int dif = B.koefs.size() - koefs.size();
           for (int i =0; i<koefs.size(); i++){
                result[dif+i] = result[dif+i] + koefs[i];
           }
        }
        return Polinom(result);
    }

    bool bigger(Polinom B, int f) {
        if (koefs.size() == B.getLength()){
            for (int i=0; i<B.getLength(); i++){
                if (koefs[i] == B.koefs[i]) {
                    continue;
                }
                else{
                    if (koefs[i]%f > B.koefs[i]%f){
                        return true;
                    }
                    else{
                        return false;
                    }
                }
            }
        }
        else {
            if(koefs.size() > B.getLength()){
                return true;
            }
            else {
                return false;
            }
        }
    return true;
    }

    Polinom division(Polinom B, int f){
        Polinom A(koefs);
        std::vector<int> result=koefs;
        int code = 0;
        while (true) {
            int cursorIndex = 1;
            if (result.size() < B.koefs.size()) {
                break;
            }
            else{
                A.koefs = result;
                result={};
                int c = division_f(A.koefs[0], B.koefs[0], f);
                if (c==-1){
                    code = -1;
                    break;
                }
                int flag = 0;
                for (int i = 1; i<B.koefs.size(); i++){
                    cursorIndex += 1;
                    if (flag == 0) {
                        if (A.koefs[i] == c*B.koefs[i]){
                            continue;
                        }
                        else {
                            flag = 1;
                            result.push_back(diff_f(A.koefs[i]%f, c*B.koefs[i]%f, f));
                        }
                    }
                    else {
                        result.push_back(diff_f(A.koefs[i]%f, c*B.koefs[i]%f, f));
                    }
                    
                }
                if (cursorIndex < A.getLength()){
                      for (cursorIndex; cursorIndex<A.getLength(); cursorIndex++){
                        result.push_back(A.koefs[cursorIndex]);
                    }
                }
            }

        }
        return Polinom(result);
    }
    

  
};



Polinom recAlg(Polinom A, Polinom B, int f){
    if ( A.getLength() != 0 && B.getLength() != 0){
        if ((A.koefs[0] == 0) || (B.koefs[0] == 0)){
            return A.sum(B, f);
        }
        if (A.bigger(B, f)) {
            return recAlg(A.division(B,f), B, f);
        }
        else{
            return recAlg(B.division(A,f), A, f);
        }
    }
    return A.sum(B, f);
}



__device__ void recAlg_norm(int* A_in, int* B_in, int f, int* nod){
        int a_len = getsize(A_in, f);
        int b_len = getsize(B_in, f);
        int* B = new int[b_len];
        int* A = new int[a_len];
        for (int i=0;i<b_len;i++) B[i] = B_in[i];
        for (int i=0;i<a_len;i++) A[i] = A_in[i];
    
    while (true){
    if (getsize(A, f) != 1 && getsize(B, f) != 1){
            if (A[0] == f+1) {
                for (int an=0; an<getsize(B,f); an++){
                    nod[an] = B[an];
                }
                break;
                }
                if (B[0] == f+1) {
                for (int an=0; an<getsize(A,f); an++){
                    nod[an] = A[an];
                }
                break;
                }
        
        if (bigger(A, B, f)) {

            int* C = new int[getsize(A, f)];
            division(A,B,f, C);
            int c_res = getsize(C,f);
            for (int c_size=0;c_size<c_res; c_size++) A[c_size] = C[c_size];


        }
        else{
            int* C = new int[getsize(B, f)];
            division(B,A,f, C);
            for (int c_size=0;c_size<getsize(C,f); c_size++) B[c_size] = C[c_size];
        }
    }
    else{
        if (A[0] == f+1) {
            for (int an=0; an<getsize(B,f); an++){
                    nod[an] = B[an];
            }
            }
        if (B[0] == f+1) {
            for (int an=0; an<getsize(A,f); an++){
                    nod[an] = A[an];
            }
            }
            break;
    }
    }
    }


__global__ void kernel_counter(int* ptr, int* input_array, int field, int count_polinoms, int size_array, int point_of_answer){

    int tid = threadIdx.x+blockIdx.x*250;
    if(tid<size_array) {
        input_array[tid]=ptr[tid];

    }
    if (tid < count_polinoms){
    while (input_array[point_of_answer-1] == field)
    {
    int i = 0;
    int counter_f = 0;
    int start_point = 0;
    int first_pol_size=0;
    int middle_point = 0;
    int end_point = 0;
    int second_pol_size = 0;
    int answer_point = 0;
    if (counter_f < tid+count_polinoms){
            while (counter_f < 2*tid){
                if (ptr[i] == field + 1){
                    counter_f++;
                }
                i++;
            }
            start_point = i;
        
        
            while (counter_f < 2*tid+1){
                if (ptr[i] == field + 1){
                    counter_f++;
                }
                i++;
                first_pol_size++;
            }
            middle_point = i;
        
            while (counter_f < 2*tid+2){
                if (ptr[i] == field + 1){
                    counter_f++;
                }
                i++;
                second_pol_size++;
            }
            end_point = i;


            while (counter_f < tid+count_polinoms){
                if (ptr[i] == field + 1){
                    counter_f++;
                }
                i++;
            }
            answer_point = i;
    }


        if (input_array[start_point] == field){
          
        continue;
        }

        if (input_array[middle_point] == field){
          
            continue;
        }

    int* A = new int[first_pol_size];
    int* B = new int[second_pol_size];
    int build_counter = start_point;
    int build_counter_a = 0;
    int build_counter_b = 0;
    while (build_counter < end_point){
        while (build_counter < middle_point){
            A[build_counter_a] = input_array[build_counter];
            build_counter++;
            build_counter_a++;
        }
        B[build_counter_b] = input_array[build_counter];
        build_counter++;
        build_counter_b++;
    }

    int result_value = getsize(A,field);
    if (result_value<getsize(B,field))
        result_value = getsize(B,field);
    int* nod = new int[result_value];
    recAlg_norm(A, B, field, nod);
    if (getsize(nod,field) == 2 && nod[0]==1){
      input_array[point_of_answer-1] =1;
    }
    for (int num =0; num<getsize(nod, field); num++){
        input_array[answer_point] = nod[num];
        answer_point++;
    }
      }
    }
    }
  



int main(){
std::vector<Polinom> unsortedPolinoms;
    std::vector<Polinom> polinoms_1;
    std::vector<Polinom> polinoms;

    int field = 0;
    std::string fieldCin;

    std::cout << "enter field\n";
    std::cin >> fieldCin;

    field = pow(2,strInt(fieldCin));

    int polinomCount;
    std::cout << "enter polinoms, press 'e' when u want to stop\n";

    int flag_field = 0;

    while (true) {
        std::string g;
        std::cin >> g;

        if (g == "e") {
            break;
        }

        std::vector<int> test_vec;
        for (int letter = 0; letter < g.length(); letter++) {
            int digit;
            digit = g[letter] - '0';
            test_vec.push_back(digit);
        }
        unsortedPolinoms.push_back(Polinom(test_vec));
    }
    
    

    int maxLength = 0; 
    
    for (Polinom i: unsortedPolinoms){
        if (maxLength < i.getLength()){
            maxLength = i.getLength();
        }
    }
    
    for (int j=0; j<=maxLength; j++){
        for (Polinom p: unsortedPolinoms){
            if (p.getLength() == j){
                polinoms_1.push_back(p);
       //         p.printPolinom();
            }
        }
    }

    int pol_count = polinoms_1.size()-1;

    int cor_pol = 0;

    while (cor_pol < pol_count-cor_pol){
            if (cor_pol == pol_count - cor_pol) polinoms.push_back(polinoms_1[cor_pol]);
            polinoms.push_back(polinoms_1[cor_pol]);
            polinoms.push_back(polinoms_1[pol_count-cor_pol]);
            cor_pol++;
        
    }

    std::cout << "sort by division result: " << "\n";
    for (Polinom p:polinoms) p.printPolinom();

    int count_symbols = 0;


    int count_polinoms = polinoms.size();

    for (int p_s=0; p_s<polinoms.size();p_s++){
        for (int pol_s=0; pol_s<polinoms[p_s].koefs.size();pol_s++){
            count_symbols ++;
        }
    }

    int all_symbols_in_result_array = (count_polinoms*2-1) + count_symbols + maxLength*(count_polinoms-1);

    int input_array[1000];


    int great_counter = 0;

    for (int p_s=0; p_s<polinoms.size();p_s++){
        for (int pol_s=0; pol_s<polinoms[p_s].koefs.size();pol_s++){
            input_array[great_counter] = polinoms[p_s].koefs[pol_s];
            great_counter++;
        }
        input_array[great_counter] = field+1;
        great_counter++;
    }

    for (int extra_pol =0; extra_pol<count_polinoms-1;extra_pol++){
        for (int field_number=0; field_number<maxLength; field_number++){
            input_array[great_counter] = field;
            great_counter++;
        }
        input_array[great_counter] = field+1;
        great_counter++;
    }


    std::cout << "max length: "<< maxLength << " polinom count: " << count_polinoms << "\n";

    int *input_for_kernel;
    int *result_after_sorting;

    int sort_result[1000];

  cudaMalloc((void**)&input_for_kernel, great_counter * sizeof(int));
  cudaMalloc((void**)&result_after_sorting, great_counter * sizeof(int));


cudaMemcpy(input_for_kernel, input_array, great_counter * sizeof(int), cudaMemcpyHostToDevice);

cudaEvent_t start, stop;
    float gpuTime = 0.0;

    cudaEventCreate( &start );
    cudaEventCreate( &stop );
    cudaEventRecord( start, 0 );

    

  kernel_counter <<<30, 250>>> (input_for_kernel, result_after_sorting, field, count_polinoms, great_counter, great_counter-maxLength);

cudaEventRecord( stop, 0 );
    cudaEventSynchronize( stop );
    cudaEventElapsedTime( &gpuTime, start, stop );
    printf("viduxa time\n", gpuTime);

    cudaEventDestroy( start );
    cudaEventDestroy( stop );


  cudaMemcpy(sort_result, result_after_sorting, great_counter * sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(input_for_kernel);
  cudaFree(result_after_sorting);

  std::cout<<"NOD = ";
  for (int i =great_counter-(maxLength+1); i<great_counter;i++){
    if (sort_result[i] == field+1) break;
    if (sort_result[i] == field) break;
    std::cout<<sort_result[i];
  }

  std::cout << "\n";

  system("pause");

  return 1;
}
