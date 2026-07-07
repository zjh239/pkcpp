//types_para.cxx
// Here is the partition to set thread number.
// I put it here so that it is low level enough.
module;

#include <thread>

export module types:para;

export {

size_t thread_number = 1;

void set_thread(size_t k = 0){
    if (k == 0){
        thread_number = std::thread::hardware_concurrency();
    } else {
        thread_number = k;
    }
}

}
