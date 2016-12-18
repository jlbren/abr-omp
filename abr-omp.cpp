#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

using namespace std;

#ifdef _OPENMP
#include <omp.h>
#endif

//unsigned 64 bit or longer int (platform dependent)
typedef unsigned long long UINT64;

typedef chrono::high_resolution_clock Clock;
/*Default Values
 * kmer = word size
 * mm = num of missmatches in index
*/
int kmer = 8;
int mm = 0;
char * index_file;
char * query_file;
vector<string> query_kmer_array;
vector<string> index_kmer_array;

ifstream in;
/*
 * count = num occurrences in index
 * pos = position in index_kmer_array
 * missmatches = pos within index related by mm missmatches
*/
struct Node {
    int count = 0;
    //vector<string*> pos ;
    //Node ** missmatches;
};
UINT64 index_size = 4;
Node * Index;

void help(){
    fprintf(stderr,"ab-omp --help|-h --kmer|-k --missmatch|-m --index|-i --query|-q\n");
    return;
}
void build_kmer_array_v(char * file, vector <string> &kmer_array){

    ifstream in;
    in.open(file);
    if(!in.is_open())
    {
        cout << "The read file could not be opened. Check the location.\n";
        throw;
    }
    int i,j,len;
    char * word=new char [kmer+1];

    cout << "Reading: " << file << endl;
    string str_word,line,header;
    getline(in,header);       //gets header

    while(in.peek()!=EOF)
    {
        len=0;
        str_word="";
        while( !(in.peek()=='>' || in.peek()==EOF))
        {
            getline(in,line);
            len+=line.size();
            str_word+=line;
        }
        


        for(i=0; i<kmer; i++)
        {
            word[i]=str_word[i];
            if(word[i]<97) word[i]+=32;                         //makes it lowercase
        }
        word[kmer]='\0';
        kmer_array.push_back(word);

        for(i=1; i<(len-kmer+1); i++)                      //read until the end of the file
        {
            //shift
            for(j=0; j<(kmer-1); j++) word[j]=word[j+1];
            word[(kmer-1)]=str_word[kmer+i-1];
            if(word[kmer-1]<97) word[kmer-1]+=32;     //makes it lowercase
            word[kmer]='\0';
           // fprintf(stderr,"Pushing back kmer: %s\n",word);
            kmer_array.push_back(word);
        }



        if(in.peek()!=EOF)
        {
            getline(in,header);       //gets header
        }
    }

    in.clear();
    in.close();

    return;


}

void compare_query(){
    UINT64 n = query_kmer_array.size();
    UINT64 num_hits = 0;
    UINT64 val;
    unsigned char temp;
    char c;

#pragma omp parallel for private(val) private(c) reduction(+:num_hits)
    for(UINT64 j = 0; j < n; ++j){
        val = 0;
        c=query_kmer_array[j].at(0);
       // fprintf(stderr,"query_kmer_array[%d]=%s\n",j,query_kmer_array[j].c_str());
       // cout<<"query_kmer_array["<<j<<"]="<<query_kmer_array[j]<<endl;

        switch(c)
        {
            case 'a':
                temp=0;
                break; //shifts 00 for a
            case 't':
                temp=1;
                break; //shifts 10 for t
            case 'c':
                temp=2;
                break;//shifts 01 for c
            case 'g':
                temp=3;
                break;//shifts 11 for g
            default: break;
        }
        val=temp;

        for (int i=1; i<kmer; ++i)
        {
            val=val<<2;
            c=query_kmer_array[j].at(i);
            switch(c)
            {
                case 'a':
                    temp=0;
                    break; //shifts 00 for a
                case 't':
                    temp=1;
                    break; //shifts 10 for t
                case 'c':
                    temp=2;
                    break;//shifts 01 for c
                case 'g':
                    temp=3;
                    break;//shifts 11 for g
                default: break;
            }
            val+=(int)temp;
        }


        if (val>index_size)
            cout <<"Out of range: " + query_kmer_array[j] <<endl;

        //printf("Hash val for %s : %u\n", word.c_str(), val);
        num_hits ++;
        //num_hits += Index[val].count;
        //printf("Num hits: %d\n", num_hits);
    }




        /* #pragma omp critical
             index[val].pos.push_back(&index_kmer_array[i]);*/

    printf("Total number of hits: %u\n",num_hits);

}
void build_kmer_array(vector<string> (&kmer_array), char * file){
    string chunk;
    char * temp = new char [kmer+1];
    char * old;
    char * in_buffer = new char[501];
    int chunk_size, count, genome_size, pos;
    genome_size = 0;
    count = 1;

    in.open(file);
    if(in.peek()=='>'){
        in.getline(in_buffer, 500);
        printf("Reading: %s\n",in_buffer);
    }
    in.read(in_buffer,kmer);
    for (int i = 0 ; i < kmer; ++i)
        temp[i] = in_buffer[i];
    temp[kmer] = '\0';
    //printf("Reading: %s\n", temp);
    kmer_array[0] = string(temp);
    while(in.peek()!=EOF){
        in_buffer = new char[501];
        if(in.peek()=='>'){
            in.getline(in_buffer, 500);
            printf("Reading: %s",in_buffer);
        }
        kmer_array.resize((500*count)+1);
        in.get(in_buffer,500,'\n');
        chunk = in_buffer;
        chunk_size = chunk.size();
        genome_size += chunk_size;
        //#pragma omp parallel for // this will not work sadly
        for (int i =0; i < chunk_size; ++i){
            old = temp;
            for(int k =0; k < kmer- 1; ++k)
                temp[k] = old[k+1];
            temp[kmer-1] = chunk[i];
            temp[kmer] ='\0';
            //printf("Reading %s\n",temp);
            pos = i + (500 * (count-1));
            kmer_array[pos] = temp;
            //printf("Kmer: %s\n",kmer_array[pos].c_str());
        }
        //implicit barrier
        count ++ ;
        delete(in_buffer);
    }
    kmer_array.resize(genome_size);
    in.clear();
    in.close();
    return;
}

void build_index(){ UINT64 num_hits = 0;
    UINT64 val;
    unsigned char temp;
    char c;
    UINT64 n = index_kmer_array.size();


    for (int i = 0; i <kmer; ++i)
        index_size *= 4;
    Index = new Node [index_size];
    printf("Index size: %u\n", index_size);
#pragma omp parallel for private(val) private(c)
    for(UINT64 j = 0; j < n; ++j){
        val = 0;
       // fprintf(stderr,"index_kmer_array[%d]=%s\n",j,index_kmer_array[j].c_str());
        c= index_kmer_array[j].at(0);
        //fprintf(stderr, "char at index_kmer_array[%d].at(0)=%c\n", j,c);

        switch(c)
        {
            case 'a':
                temp=0;
                break; //shifts 00 for a
            case 't':
                temp=1;
                break; //shifts 10 for t
            case 'c':
                temp=2;
                break;//shifts 01 for c
            case 'g':
                temp=3;
                break;//shifts 11 for g
            default: break;
        }
        val=temp;

        for (int i=1; i<kmer; ++i)
        {
            val=val<<2;
            c=index_kmer_array[j].at(i);
            switch(c)
            {
                case 'a':
                    temp=0;
                    break; //shifts 00 for a
                case 't':
                    temp=1;
                    break; //shifts 10 for t
                case 'c':
                    temp=2;
                    break;//shifts 01 for c
                case 'g':
                    temp=3;
                    break;//shifts 11 for g
                default: break;
            }
            val+=(int)temp;
        }


        if (val>index_size)
            cout <<"Out of range: " + query_kmer_array[j] <<endl;

        //printf("Hash val for %s : %u\n", word.c_str(), val);
#pragma omp atomic
        Index[val].count++;

    }
    printf("\nIndex successfully built\n");
    return;
}

int main(int argc, char* argv[]) {

    int max_threads = 1;
    #ifdef _OPENMP
        max_threads = omp_get_max_threads();
        fprintf(stderr, "OpenMP enabled\n");
        #pragma omp parallel
        {
            int thread_num = omp_get_thread_num();
            fprintf(stderr,"Thread num %d checking in\n", thread_num);
        }
    #endif

    if (argc < 3){
        fprintf(stderr, "Min number of arguments not met. Index and query files must be specified.\n");
        help();
        return 1;
    }
    if (argc ==3){
        index_file = argv[1];
        query_file = argv[2];
    }
    else {
        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                help();
                return 1;
            }
            else if (strcmp(argv[i], "--kmer") == 0 || strcmp(argv[i], "-k") == 0) {
                i++;
                if (isdigit(*argv[i]))
                    kmer = atoi(argv[i]);
            }
            else if (strcmp(argv[i], "--missmatch") == 0 || strcmp(argv[i], "-m") == 0) {
                i++;
                if (isdigit(*argv[i])) if (atoi(argv[i]) <= kmer)
                        mm = atoi(argv[i]);
            }
            else if (strcmp(argv[i], "--index") == 0 || strcmp(argv[i], "-i") == 0) {
                index_file = argv[i];
            }
            else if (strcmp(argv[i], "--query") == 0 || strcmp(argv[i], "-q") == 0) {
                query_file = argv[i];
            }
        }
    }

    in.open(index_file);
    if(!in){
        fprintf(stderr, "Invalid first argument, cannot open index file: %s", index_file);
        return 1;
    }
    in.close();
    in.open(query_file);
    if(!in){
        fprintf(stderr, "Invalid second argument, cannot open query file: %s", query_file);
        return 1;
    }
    in.close();

    fprintf(stderr, "Max threads: %d\nIndex: %s\nQuery: %s\nKmer: %d Missmatches: %d\n",
            max_threads, index_file, query_file, kmer, mm);

    auto t1 = Clock :: now();
    build_kmer_array_v(index_file,index_kmer_array);
    auto t2 = Clock :: now();
    auto time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    printf("Index kmer array size: %d\n", index_kmer_array.size());
    printf("Time (ms) to build kmer array for %s: %f\n", index_file,
            time_span.count()*1000);
   /* for(int i = 0; i < index_kmer_array.size(); ++i)
        printf("Index %d = %s\n",i,index_kmer_array[i].c_str());*/

    t1 = Clock :: now();
    build_index();
    t2 = Clock :: now();
    time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    printf("Time (ms) to build index: %f\n",
           time_span.count()*1000);

    t1 = Clock :: now();
    build_kmer_array_v(query_file,query_kmer_array);
    t2 = Clock :: now();
    time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    printf("Query size: %d\n", query_kmer_array.size());
    printf("Time (ms) to build kmer array for %s: %f\n", query_file,
           time_span.count()*1000);

    t1 = Clock :: now();
    compare_query();
    t2 = Clock :: now();
    time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
    printf("Time (ms) to map hits: %f\n", time_span.count()*1000);


}

